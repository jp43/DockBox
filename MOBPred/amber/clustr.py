import os
import stat
import shutil
import subprocess
import argparse

import minimz as mn
from MOBPred.tools import mol2

def do_clustering(files_r, files_l, cutoff=2.0, mode='clustering', cleanup=True):
    """
    do_clustering(files_r, files_l=None)

    Performs Amber's cpptraj clustering

    Parameters
    ----------
    files_r: filenames for receptor (.pdb)
    files_l: list of filenames (.mol2) for ligand, when ligand-protein complex

    Steps
    -----
    antechamber, parmchk (ligand-protein complex) 
"""
    # get current directory
    curdir = os.getcwd()

    # create directory where minimization will be performed
    workdir = mode
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of receptor files
    if len(files_r) == 1:
        nfiles_l = len(files_l)
        files_receptor = [os.path.abspath(files_r[0]) for idx in range(nfiles_l)]
    else:
        files_receptor = []
        for file_r in files_r:
            files_receptor.append(os.path.abspath(file_r))

    # get full path of ligand files
    if len(files_l) > 2:
        files_ligand = []
        for file_l in files_l:
            files_ligand.append(os.path.abspath(file_l))
    else:
        raise ValueError('At least 2 ligand files are required for clustering')

    # Check if same number of receptors and ligands 
    if len(files_receptor) != len(files_ligand):
        raise ValueError('Number of receptors and ligands should be the same!')

    # change working directory
    os.chdir(workdir)

    new_files_receptor = []
    # prepare receptors
    for idx, file_r in enumerate(files_receptor):
        new_file_r = 'protein-%s.pdb'%idx
        # prepare receptor
        mn.prepare_receptor(new_file_r, file_r, False)
        new_files_receptor.append(new_file_r)

    # amber clustering
    do_amber_clustering(new_files_receptor, files_ligand, cutoff, mode, cleanup=cleanup)
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

    linespdb = linespdb[:-1]
    with open(filename, 'w') as ff:
        contents ="""source %(forcefield)s
source leaprc.gaff
loadamberparams frcmod.ionsjc_tip3p
loadamberparams frcmod.ionslm_1264_tip3p
LIG = loadmol2 ligand_ref.mol2
loadamberparams ligand.frcmod
%(linespdb)s
quit"""% locals()
        ff.write(contents)

def prepare_cpptraj_config_file(filename, files_rl, cutoff, mode='clustering'):

    lines_trajin = ""
    for file_rl in files_rl:
        lines_trajin += "trajin %s\n"%(file_rl)

    # remove last \n
    lines_trajin = lines_trajin[:-1]

    # write cpptraj config file to cluster frames
    with open(filename, 'w') as file:
        if mode == 'clustering':
            contents ="""parm protein-ligand.prmtop
%(lines_trajin)s
rms first "@CA,C,N,O & !:LIG"
cluster ":LIG & !@/H" nofit mass epsilon %(cutoff)s summary summary.dat info info.dat\n"""% locals()
            file.write(contents)
        elif mode == 'pca':
            contents ="""parm protein-ligand.prmtop
%(lines_trajin)s
createcrd md-trajectories
run
crdaction md-trajectories matrix covar name covar :LIG
runanalysis diagmatrix covar out evecs.dat vecs 2 name myEvecs
crdaction md-trajectories projection md-pca modes myEvecs :LIG out pca.out\n"""% locals()
            file.write(contents)

#crdaction md-trajectories rms @CA,C,N,O&!:LIG

def do_amber_clustering(files_r, files_l, cutoff, mode, cleanup=False):

    # (A) Prepare ligand and PDB files
    os.mkdir('PDB')
    files_rl = []
    for idx, file_l in enumerate(files_l):
        mol2.update_mol2file(file_l, 'ligand.mol2', ligname='LIG')
        if len(files_r) != 1:
            file_r = files_r[idx]
        else:
            file_r = files_r[0]
        file_rl = 'PDB/protein-ligand-%s.pdb'%(idx+1)
        mn.prepare_ligand(file_r, 'ligand.mol2', file_rl)
        files_rl.append(file_rl)
        os.remove(file_r)
        if idx == 0:
            shutil.copyfile('ligand.mol2','ligand_ref.mol2')

    # (B) Run tleap
    prepare_leap_config_file('leap.in', files_r, files_l, files_rl)
    subprocess.check_output('tleap -f leap.in > leap.log', shell=True, executable='/bin/bash')

    # (C) Run cpptraj
    prepare_cpptraj_config_file('cpptraj.in', files_rl, cutoff, mode=mode)
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
        default=2.0,
        help = 'RMSD cutoff for clustering analysis')

    parser.add_argument('-mode',
        type=str,
        dest='mode',
        default='clustering',
        help = 'Cpptraj mode (clustering, pca)')

    parser.add_argument('-cleanup',
        dest='cleanup',
        action='store_true',
        default=False,
        help = 'Remove intermediate files')

    return parser

def run():

    parser = create_arg_parser()
    args = parser.parse_args()

    do_clustering(args.files_r, args.files_l, cutoff=args.cutoff, mode=args.mode, cleanup=args.cleanup)
