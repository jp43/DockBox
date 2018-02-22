import os
import stat
import shutil
import subprocess
import argparse

import ambertools
from MOBPred.tools import mol2

def do_clustering(files_r, files_l, mode='clustering', cutoff=None, nclusters=None, cleanup=True, mask=None, maskfit=None):
    """
    do_clustering(files_r, files_l=None)

    Performs Amber's cpptraj clustering

    Parameters
    ----------
    files_r: filenames for receptor (.pdb)
    files_l: list of filenames (.mol2) for ligand, when ligand-protein complex

"""
    # get current directory
    curdir = os.getcwd()

    # create directory where minimization will be performed
    workdir = mode
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of ligand files
    if len(files_l) > 2:
        files_l = [os.path.abspath(file_l) for file_l in files_l]
    else:
        raise ValueError('At least 2 ligand files are required for clustering')

    # get full path of receptor files
    if len(files_r) == 1:
        nfiles_l = len(files_l)
        files_r = [os.path.abspath(files_r[0]) for idx in range(nfiles_l)]
    else:
        files_r = [os.path.abspath(file_r) for file_r in files_r]

    # Check if same number of receptors and ligands 
    if len(files_r) != len(files_l):
        raise ValueError('Number of receptors and ligands should be the same!')

    # change working directory
    os.chdir(workdir)

    new_files_r = []
    # prepare receptors
    for idx, file_r in enumerate(files_r):
        new_file_r = 'protein-%s.pdb'%idx
        # prepare receptor
        ambertools.prepare_receptor(new_file_r, file_r)
        new_files_r.append(new_file_r)

    # amber clustering
    do_amber_clustering(new_files_r, files_l, mode, cutoff=cutoff, nclusters=nclusters, cleanup=cleanup, mask=mask, maskfit=maskfit)
    os.chdir(curdir)

def prepare_leap_config_file(filename, files_r, files_l, files_rl, forcefield='leaprc.ff14SB'):

    linespdb = ""
    for idx, file_rl in enumerate(files_rl):
        if idx == 0:
            linespdb += """p = loadPdb %s
saveAmberParm p protein-ligand.prmtop protein-ligand.inpcrd
savepdb p %s\n"""%(file_rl,file_rl)
        else:
            linespdb += """p = loadPdb %s
savepdb p %s\n"""%(file_rl,file_rl)

    name = ambertools.get_ligand_name(files_l[0])

    linespdb = linespdb[:-1]
    with open(filename, 'w') as ff:
        contents ="""source %(forcefield)s
source leaprc.gaff
loadamberparams frcmod.ionsjc_tip3p
loadamberparams frcmod.ionslm_1264_tip3p
%(name)s = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
%(linespdb)s
quit"""% locals()
        ff.write(contents)

def prepare_cpptraj_config_file(filename, files_rl, cutoff=None, nclusters=None, mode='clustering', mask=None, maskfit=None, prmtop='protein-ligand.prmtop'):

    lines_trajin = ""
    for file_rl in files_rl:
        lines_trajin += "trajin %s\n"%(file_rl)

    # remove last \n
    lines_trajin = lines_trajin[:-1]

    # write cpptraj config file to cluster frames
    with open(filename, 'w') as file:
        if mode == 'clustering':
            if cutoff and nclusters:
                ValueError('Both cutoff value and nclusters provided. Only one of those parameters should be given!')
            elif cutoff:
                option = " epsilon %s "%cutoff
            elif nclusters:
                option = " clusters %s"%nclusters
            else:
                option = ""

            contents = """parm %(prmtop)s
%(lines_trajin)s
rms first %(maskfit)s
cluster %(mask)s nofit%(option)s summary summary.dat info info.dat\n"""% locals()
            file.write(contents)
        elif mode == 'pca':
            contents = """parm %(prmtop)s
%(lines_trajin)s
rms first %(maskfit)s
createcrd md-trajectories
run
crdaction md-trajectories matrix covar name covar %(mask)s
runanalysis diagmatrix covar out evecs.dat vecs 2 name myEvecs
crdaction md-trajectories projection md-pca modes myEvecs %(mask)s out pca.out\n"""% locals()
            file.write(contents)
        elif mode == 'fit':
            contents = """parm %(prmtop)s
%(lines_trajin)s
rms first %(maskfit)s
trajout ref.rst restart onlyframes 1
trajout struct.pdb multi\n"""% locals()
            file.write(contents)
        elif mode == 'rmsd':
            lines_rms = ""
            for idx, file_rl in enumerate(files_rl):
                jdx = idx + 1
                lines_rms += """reference %(file_rl)s [ref%(jdx)s]
rms ref [ref%(jdx)s] %(maskfit)s 
rms ref [ref%(jdx)s] %(mask)s nofit out rmsd_%(jdx)s.txt\n""" %locals()
            contents = """parm protein-ligand.prmtop
%(lines_trajin)s
%(lines_rms)s"""% locals()
            file.write(contents)

def do_amber_clustering(files_r, files_l, mode, cutoff=None, nclusters=None, cleanup=False, mask=None, maskfit=None):

    # (A) Prepare ligand and PDB files
    os.mkdir('PDB')
    files_rl = []
    for idx, file_l in enumerate(files_l):
        file_r = files_r[idx]
        file_rl = 'PDB/protein-ligand-%s.pdb'%(idx+1)
        shutil.copyfile(file_l, 'ligand.mol2')
        ambertools.prepare_ligand(file_r, 'ligand.mol2', file_rl)
        files_rl.append(file_rl)
        os.remove(file_r)
        if idx == 0:
            shutil.copyfile('ligand.mol2', 'ligand_ref.mol2')

    if not mask or not maskfit:
        ligname = ambertools.get_ligand_name('ligand_ref.mol2')
        if not mask:
            mask = ':%s&!@H='%ligname
        if not maskfit:
            maskfit = '@CA,C,N,O&!:%s'%ligname 

    # (B) Run tleap
    prepare_leap_config_file('leap.in', files_r, files_l, files_rl)
    subprocess.check_output('tleap -f leap.in > leap.log', shell=True, executable='/bin/bash')

    # (C) Run cpptraj
    prepare_cpptraj_config_file('cpptraj.in', files_rl, cutoff=cutoff, nclusters=nclusters, mode=mode, mask=mask, maskfit=maskfit)
    subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True, executable='/bin/bash')

    if cleanup:
        # (D) remove PDB folder and other large files
        shutil.rmtree('PDB', ignore_errors=True)
        os.remove('leap.log')
        os.remove('protein-ligand.prmtop')

def create_arg_parser():

    parser = argparse.ArgumentParser(description="Run Amber Clustering")

    parser.add_argument('-l',
        type=str,
        dest='files_l',
        nargs='+',
        default=None,
        help = 'Ligand coordinate file(s): .mol2')

    parser.add_argument('-r',
        type=str,
        dest='files_r',
        nargs='+',
        required=True,
        help = 'Receptor coordinate file(s): .pdb')

    parser.add_argument('-rmsd',
        type=str,
        dest='cutoff',
        default=None,
        help = 'RMSD cutoff for clustering analysis')

    parser.add_argument('-n',
        type=str,
        dest='nclusters',
        default=None,
        help = 'Number of clusters for clustering analysis')

    parser.add_argument('-mode',
        type=str,
        dest='mode',
        default='clustering',
        help = 'Cpptraj mode (clustering, pca, fit, rmsd)')

    parser.add_argument('-cleanup',
        dest='cleanup',
        action='store_true',
        default=False,
        help = 'Remove intermediate files')

    parser.add_argument('-m',
        dest='mask',
        default='',
        help = 'Mask used for clustering or pca')

    parser.add_argument('-mf',
        dest='maskfit',
        default='',
        help = 'Mask used for fitting prior to clustering or pca')

    return parser

def run():

    parser = create_arg_parser()
    args = parser.parse_args()

    do_clustering(args.files_r, args.files_l, mode=args.mode, cutoff=args.cutoff, nclusters=args.nclusters, cleanup=args.cleanup, mask=args.mask, maskfit=args.maskfit)
