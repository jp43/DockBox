import os
import stat
import shutil
import subprocess

import minimz as mn
from MOBPred.tools import mol2

def do_clustering(files_r, files_l, cutoff=2.0, cleanup=True):
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
    workdir = 'clustr'
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of receptor files
    if isinstance(files_r, str):
        files_r = [files_r]
    r_tmp = []
    for file_r in files_r:
        r_tmp.append(os.path.abspath(file_r))
    files_r = r_tmp

    # get full path of ligand files
    if len(files_l) > 2:
        l_tmp = []
        for file_l in files_l:
            l_tmp.append(os.path.abspath(file_l))
        files_l = l_tmp
    else:
        raise ValueError('At least 2 ligand files are required for clustering')

    # Check if same number of receptors and ligands 
    if len(files_r) != 1:
        if len(files_r) != len(files_l):
            raise ValueError('Number of receptors and ligands should be the same!')

    # change working directory
    os.chdir(workdir)

    new_files_r = []
    # prepare receptors
    for idx, file_r in enumerate(files_r):
        new_file_r = 'rec-%s.pdb'%idx
        # prepare receptor
        mn.prepare_receptor(new_file_r, file_r, False)
        new_files_r.append(new_file_r)

    # amber clustering
    do_amber_clustering(new_files_r, files_l, cutoff, cleanup=cleanup)
    os.chdir(curdir)

def prepare_ligand(file_r, file_l, file_rl):

    # Use antechamber to prepare the ligand
    script_name = 'prepare_ligand.sh'
    with open(script_name, 'w') as ff:
        script ="""#!/bin/bash
# prepare mol2 file with antechamber
antechamber -fi mol2 -i %(file_l)s -fo mol2 -o tmp.mol2 -c gas > antchmb.log
mv tmp.mol2 %(file_l)s
if [ ! -f lig.frcmod ]; then
  parmchk -i %(file_l)s -f mol2 -o lig.frcmod # run parmchk
fi 

# prepare complex.pdb
antechamber -fi mol2 -i %(file_l)s -fo pdb -o lig.pdb > /dev/null # create pdb
cp %(file_r)s %(file_rl)s\n"""%locals()
        ff.write(script)

    os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)
    subprocess.check_output('./' + script_name, shell=True, executable='/bin/bash')

    with open('lig.pdb', 'r') as ligf:
        with open(file_rl, 'a') as cmplxf:
            for line in ligf:
                if line.startswith(('ATOM','HETATM','TER')):
                    cmplxf.write(line)

    os.remove(script_name)

def prepare_leap_config_file(filename, files_r, files_l, files_rl):

    linespdb = ""
    for idx, file_rl in enumerate(files_rl):
        if idx == 0:
            linespdb += """p = loadPdb %s
saveAmberParm p rec-lig.prmtop rec-lig.inpcrd
savepdb p %s\n"""%(file_rl,file_rl)
        else:
            linespdb += """p = loadPdb %s
savepdb p %s\n"""%(file_rl,file_rl)

    linespdb = linespdb[:-1]
    with open(filename, 'w') as ff:
        contents ="""source leaprc.ff12SB
source leaprc.gaff
loadamberparams frcmod.ionsjc_tip3p
loadamberparams frcmod.ionslm_1264_tip3p
LIG = loadmol2 ligref.mol2
loadamberparams lig.frcmod
%(linespdb)s
quit"""% locals()
        ff.write(contents)

def prepare_cpptraj_config_file(filename, files_rl, cutoff):

    lines_trajin = ""
    for file_rl in files_rl:
        lines_trajin += "trajin %s\n"%(file_rl)

    # remove last \n
    lines_trajin = lines_trajin[:-1]

    # write cpptraj config file to cluster frames
    with open(filename, 'w') as file:
        contents ="""parm rec-lig.prmtop
%(lines_trajin)s
rms first "@CA,C,N & !:LIG"
cluster ":LIG & !@/H" nofit mass epsilon %(cutoff)s summary summary.dat info info.dat
"""% locals()
        file.write(contents)

def do_amber_clustering(files_r, files_l, cutoff, cleanup=True):

    # (A) Prepare ligand and PDB files
    os.mkdir('PDB')
    files_rl = []
    for idx, file_l in enumerate(files_l):
        mol2.update_mol2file(file_l, 'lig.mol2', ligname='LIG')
        if len(files_r) != 1:
            file_r = files_r[idx]
        else:
            file_r = files_r[0]
        file_rl = 'PDB/rec-lig-%s.pdb'%(idx+1)
        prepare_ligand(file_r, 'lig.mol2', file_rl)
        files_rl.append(file_rl)
        os.remove(file_r)
        if idx == 0:
            shutil.copyfile('lig.mol2','ligref.mol2')
        
    # (B) Run tleap
    prepare_leap_config_file('leap.in', files_r, files_l, files_rl)
    subprocess.check_output('tleap -f leap.in > leap.log', shell=True, executable='/bin/bash')

    # (C) Run cpptraj
    prepare_cpptraj_config_file('cpptraj.in', files_rl, cutoff)
    subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True, executable='/bin/bash')

    if cleanup:
        # (D) remove PDB folder and other large files
        shutil.rmtree('PDB', ignore_errors=True)
        os.remove('leap.log')
        os.remove('rec-lig.prmtop')
