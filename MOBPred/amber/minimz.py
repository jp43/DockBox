import os
import sys
import stat
import shutil
import subprocess
import argparse

from MOBPred.tools import mol2
from MOBPred.tools import reader
import ambertools as ambt

leap_default_settings = {'solvate': False, 'PBRadii': None, 'forcefield': 'leaprc.ff14SB'}

def do_minimization(file_r, files_l=None, restraints=None, keep_hydrogens=False):
    """
    do_minimization(file_r, files_l=None, keep_hydrogens=False)

    Performs Amber minimization

    Parameters
    ----------
    file_r: filename for receptor (.pdb)
    files_l: list of filenames (.mol2) for ligand, when ligand-protein complex

    Steps
    -----
    antechamber, parmchk (ligand-protein complex) 
"""
    # get current directory
    curdir = os.getcwd()

    # create directory where minimization will be performed
    workdir = 'minimz'
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of ligand and receptor files
    file_r = os.path.abspath(file_r)
    if files_l:
        if isinstance(files_l, str):
            files_l = [files_l] # make it a list
        l_tmp = []
        for file_l in files_l:
            l_tmp.append(os.path.abspath(file_l))
        files_l = l_tmp

    # change working directory
    os.chdir(workdir)

    # only keep atom lines
    with open(file_r, 'r') as tmpf:
        with open('rec.pdb', 'w') as recf:
            for line in tmpf:
                # check if non-hydrogen atom line
                if line.startswith(('ATOM','TER')):
                    recf.write(line)
            # if last line not TER, write it
            if not line.startswith('TER'):
                recf.write('TER\n')

    # prepare receptor
    ambt.prepare_receptor('rec.pdb', file_r, keep_hydrogens=keep_hydrogens)

    if not isinstance(restraints, list):
        restraints = [restraints]

    # amber minimization
    do_amber_minimization('rec.pdb', files_l, restraints=restraints, keep_hydrogens=keep_hydrogens)
    os.chdir(curdir)

def prepare_minimization_config_file(script_name, restraint=None):

    if restraint:
        restraint_lines = """\nibelly=1,
 bellymask='%(restraint)s',"""%locals()
    else:
        restraint_lines = ""

    with open(script_name, 'w') as minf:
        script ="""In-Vacuo minimization with restraints
&cntrl
 imin=1, maxcyc=5000,
 ntb=0,
 ncyc=1000,
 ntmin=1,
 ntpr=5,
 cut=10.0,%(restraint_lines)s
&end\n"""%locals()
        minf.write(script)

def prepare_and_minimize(restraints, keep_hydrogens=False):

    # run tleap
    subprocess.check_output('tleap -f leap.in > /dev/null', shell=True, executable='/bin/bash')
    shutil.copyfile('start.inpcrd', 'start.rst')

    for idx, restraint in enumerate(restraints):
        infile = 'min%i.in'%idx
        outfile = 'min%i.out'%idx 
        prepare_minimization_config_file(infile, restraint)
        try:
            # run minimization
            subprocess.check_output('sander -O -i %s -o %s -c start.inpcrd -p start.prmtop -ref start.rst -r end.inpcrd > /dev/null'%(infile, outfile), shell=True, executable='/bin/bash')
            shutil.move('end.inpcrd', 'start.inpcrd')
            status = 0 # minimization finished normaly
        except subprocess.CalledProcessError as e:
            # the minimization failed
            return e.returncode

    # get output configuration
    subprocess.check_output('cpptraj -p start.prmtop -y start.inpcrd -x complex_out.pdb > /dev/null', shell=True, executable='/bin/bash')
    return status

def do_amber_minimization(file_r, files_l, restraints=None, keep_hydrogens=False):

    if files_l:
        for idx, file_l in enumerate(files_l):
            shutil.copyfile(file_l, 'lig.mol2')

            # change ligand name to LIG
            ligname = reader.open('lig.mol2').ligname
            mol2.update_mol2file('lig.mol2', 'lig.mol2', ligname='LIG')

            if idx == 0:
                ambt.prepare_leap_config_file('leap.in', file_r, 'lig.mol2', 'complex.pdb')

            # since the charges of the original .mol2 file are eventually kept using pdb2mol2
            # do not regenerate the charges using antechamber when preparing the ligand
            ambt.prepare_ligand(file_r, 'lig.mol2', 'complex.pdb', charge_method='None')
            status = prepare_and_minimize(restraints, keep_hydrogens=keep_hydrogens)

            if status == 0:
                is_ligand = False
                # get ligand atom positions from complex file
                with open('complex_out.pdb', 'r') as cf:
                    with open('rec-%s.out.pdb'%(idx+1), 'w') as recf:
                        with open('lig.out.pdb', 'w') as tmpf:
                            for line in cf:
                                if line.startswith(('ATOM', 'HETATM')):
                                    if line[17:20] == 'LIG':
                                        is_ligand = True
                                if is_ligand:
                                    tmpf.write(line)
                                else:
                                    recf.write(line)
                if not is_ligand:
                    raise ValueError('Ligand not found in complex_out.pdb')

                mol2file = 'lig-%s.out.mol2'%(idx+1)
                mol2.pdb2mol2('lig.out.pdb', mol2file, file_l)
                mol2.update_mol2file(mol2file, mol2file, ligname=ligname)
    else:
        ambt.prepare_leap_config_file('leap.in', file_r, None, 'complex.pdb')
        prepare_and_minimize(restraints, keep_hydrogens=keep_hydrogens)
        shutil.move('complex_out.pdb', 'rec.out.pdb')

def create_arg_parser():

    parser = argparse.ArgumentParser(description="Run Amber Minimization")

    parser.add_argument('-l',
        type=str,
        dest='input_file_l',
        nargs='+',
        default=None,
        help = 'Ligand coordinate file(s): .mol2')

    parser.add_argument('-r',
        type=str,
        dest='input_file_r',
        required=True,
        help = 'Receptor coordinate file(s): .pdb')

    parser.add_argument('-restraints',
        type=str,
        dest='restraints',
        default=None,
        help = 'Restraints')

    parser.add_argument('-removeh',
        dest='removeh',
        action='store_true',
        help = 'Remove Hydrogens')

    return parser

def run():

    parser = create_arg_parser()
    args = parser.parse_args()

    do_minimization(args.input_file_r, files_l=args.input_file_l, restraints=args.restraints, keep_hydrogens=not args.removeh)
