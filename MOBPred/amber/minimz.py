import os
import sys
import stat
import shutil
import subprocess
import argparse

from MOBPred.tools import mol2
from MOBPred.tools import reader

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
    prepare_receptor('rec.pdb', file_r, keep_hydrogens=keep_hydrogens)

    if not isinstance(restraints, list):
        restraints = [restraints]

    # amber minimization
    do_amber_minimization('rec.pdb', files_l, restraints=restraints, keep_hydrogens=keep_hydrogens)
    os.chdir(curdir)

def prepare_receptor(file_r_out, file_r, keep_hydrogens=False):

    # only keep atom lines
    with open(file_r, 'r') as tmpf:
        with open(file_r_out, 'w') as recf:
            for line in tmpf:
                # check if non-hydrogen atom line
                if line.startswith(('ATOM', 'HETATM', 'TER')):
                    recf.write(line)
            # if last line not TER, write it
            if not line.startswith('TER'):
                recf.write('TER\n')

    # remove atoms and hydrogen with no name recognized by AMBER
    correct_hydrogen_names(file_r_out, keep_hydrogens=keep_hydrogens)

def run_antechamber(infile, outfile, at='gaff', c='gas'):
    """ use H++ idea of running antechamber multiple times with bcc's 
charge method to estimate the appropriate net charge!!"""

    logfile = 'antchmb.log'
    max_net_charge = 30
    net_charge = [0]

    if c.lower() != 'none':
        for nc in range(max_net_charge):
            net_charge.extend([nc+1,-(nc+1)])

        for nc in net_charge:
            iserror = False
            cmd = 'antechamber -i %(infile)s -fi mol2 -o %(outfile)s -fo mol2 -at %(at)s -c %(c)s -nc %(nc)s -du y -pf y > %(logfile)s'%locals()
            subprocess.call(cmd, shell=True)
            with open(logfile, 'r') as lf:
                for line in lf:
                    line_st = line.strip()
                    line_sp = line_st.split()
                    if 'Warning' in line or 'Error' in line:
                        iserror = True
                    if line_st.startswith("does not equal"):
                        nc_suggested = int(float(line_sp[8][1:-1]))
                        if nc_suggested == nc:
                            return
                if not iserror:
                    break
        if iserror:
            raise ValueError("No appropriate net charge was found to run antechamber's %s charge method"%c)
    else: # do not regenerate charges
       cmd = 'antechamber -i %(infile)s -fi mol2 -o %(outfile)s -fo mol2 -at %(at)s -du y -pf y > %(logfile)s'%locals()
       subprocess.call(cmd, shell=True)

def correct_hydrogen_names(file_r, keep_hydrogens=False):

    chainIDs = []
    atoms_info = load_PROTON_INFO()

    nremoved = 0
    removed_lines = []

    chainID = None
    is_first_residue = True
    first_residues = []

    # determine which residues are first residues
    with open(file_r, 'r') as rf:
        for line in rf:
            if line.startswith('ATOM'): # atom line
                resnum = line[22:26].strip()

                if is_first_residue:
                    first_residues.append(resnum)
                    is_first_residue = False

            elif line.startswith('TER'):
                is_first_residue = True

    resnum = ''
    with open(file_r, 'r') as rf:
        with open('tmp.pdb', 'w') as wf:
            for line in rf:
                remove_line = False
                if line.startswith('ATOM'): # atom line
                    resname = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    chainID = line[21:22].strip()
                    resnum = line[22:26].strip()

                    if resname in atoms_info:
                       # atom (if atom name starts with a digit, correct it)
                        if atom_name[0].isdigit():
                            atom_name = atom_name[1:] + atom_name[0]
 
                        # check if hydrogen should be removed
                        if atom_name[0] == 'H':
                            is_hydrogen_from_nterminal = resnum in first_residues and atom_name == 'H'
                            is_hydrogen_known = atom_name in atoms_info[resname] and not is_hydrogen_from_nterminal
                            if keep_hydrogens and not is_hydrogen_known:
                                #print line
                                remove_line = True
                                removed_lines.append(line)
                                #print hydrogens_info[resname], atom_name
                                nremoved += 1
                            elif not keep_hydrogens:
                                remove_line = True
                                nremoved += 1
                        # check if non-hydrogen atom should be removed
                        else:
                            is_atom_known = atom_name in atoms_info[resname]
                            if not is_atom_known:
                                remove_line = True
                                removed_lines.append(line)
                                nremoved += 1

                if not remove_line:
                    wf.write(line)
    #print '\n'.join(removed_lines)
    shutil.move('tmp.pdb', file_r)
    #print "Number of atom lines removed: %s" %nremoved

def load_PROTON_INFO():

    filename = os.path.dirname(os.path.abspath(__file__)) + '/PROTON_INFO'
    info = {}

    with open(filename) as ff:
        ff.next() # skip first line
        for line in ff:
            line_s = line.split()
            is_residue_line = len(line_s) == 2 and line_s[1].isdigit()
            is_hydrogen_line = len(line_s) >= 4 and \
                all([c.isdigit() for c in line_s[:4]])
            is_heavy_atom_line = not is_residue_line and \
                not line_s[0].isdigit()

            if is_residue_line:
                resname = line_s[0]
                info[resname] = []
            elif is_hydrogen_line:
                info[resname].extend(line[15:].split())
            elif is_heavy_atom_line:
                info[resname].extend(line_s)

    no_h_residues = ['PRO']
    for resname in info:
        if resname not in no_h_residues:
            info[resname].append('H')

    info['NME'] = []
    return info

def load_atomic_ions():
    """Load formal charge libraries of monoatomic ions"""

    filename = os.path.dirname(os.path.abspath(__file__)) + '/atomic_ions.cmd'
    info = {}
    with open(filename) as ff:
        for line in ff:
            if line.startswith('i = createAtom'):
                charge = float(line.split()[5])
                is_new_atom = True
            elif line.startswith('r = createResidue') and is_new_atom:
                resname = line.split()[3]
                info[resname] = charge
            elif not line.strip():
                is_new_atom = False
    return info

def prepare_ligand(file_r, file_l, file_rl, charge_method='gas'):

    file_l_prefix, ext = os.path.splitext(file_l)
    file_l_prefix = os.path.basename(file_l_prefix)

    run_antechamber(file_l, 'tmp.mol2', at='gaff', c=charge_method)
    shutil.move('tmp.mol2', file_l)
    subprocess.check_output('parmchk -i %s -f mol2 -o %s.frcmod'%(file_l,file_l_prefix), shell=True, executable='/bin/bash')
    subprocess.check_output('antechamber -fi mol2 -i %s -fo pdb -o %s.pdb > /dev/null'%(file_l,file_l_prefix), shell=True, executable='/bin/bash')

    shutil.copyfile(file_r, file_rl)
    subprocess.check_output('cat %s.pdb >> %s'%(file_l_prefix, file_rl), shell=True, executable='/bin/bash')


def get_ions_number(logfile, concentration=0.15):

    with open(logfile, 'r') as lf:
        for line in lf:
            line_s = line.strip().split()
            if len(line_s) > 2:
                if line_s[0] == 'Added' and line_s[2] == 'residues.': 
                    nresidues = int(line_s[1])
                    ncl = int(round(nresidues * concentration * 0.0187))
                    nna = ncl
            if line.startswith("WARNING: The unperturbed charge"):
                net_charge = int(round(float(line_s[7])))
                if net_charge > 0:
                    ncl += abs(net_charge)
                elif net_charge < 0:
                    nna += abs(net_charge)
    return nna, ncl

def prepare_leap_config_file(script_name, file_r, file_l, file_rl, solvate=False, PBRadii=None, forcefield='leaprc.ff14SB', nna=0, ncl=0, distance=10.0, closeness=1.0, remove=None):

    lines_solvate = ""
    if solvate:
        lines_solvate = "\nsolvateBox complex TIP3PBOX %.2f %.2f"%(distance,closeness)

    lines_pbradii = ""
    if PBRadii:
        lines_pbradii = "\nset default PBRadii %s"%PBRadii

    remove_lines = ""
    if remove:
        for idx in remove:
            remove_lines += "\nremove complex complex.%i"%idx 

    lines_ions = ""
    if nna > 0:
        lines_ions += "\naddions complex Na+ %i"%nna
    if ncl > 0:
        lines_ions += "\naddions complex Cl- %i"%ncl

    if file_l:
        file_l_prefix, ext = os.path.splitext(file_l)
        file_l_prefix = os.path.basename(file_l_prefix)
        with open(script_name, 'w') as leapf:
                script ="""source %(forcefield)s
source leaprc.gaff
loadoff atomic_ions.lib
loadamberparams frcmod.ionsjc_tip3p
loadamberparams frcmod.ionslm_1264_tip3p
LIG = loadmol2 %(file_l)s
loadamberparams %(file_l_prefix)s.frcmod%(lines_pbradii)s
complex = loadPdb %(file_rl)s%(lines_solvate)s%(lines_ions)s%(remove_lines)s
saveAmberParm complex start.prmtop start.inpcrd
savePdb complex start.pdb
quit\n"""%locals()
                leapf.write(script)
    else:
        with open(script_name, 'w') as leapf:
                script ="""source %(forcefield)s
loadoff atomic_ions.lib
loadamberparams frcmod.ionsjc_tip3p
loadamberparams frcmod.ionslm_1264_tip3p
complex = loadPdb %(file_r)s%(lines_solvate)s%(lines_ions)s
saveAmberParm complex start.prmtop start.inpcrd
savePdb complex start.pdb
quit\n"""%locals()
                leapf.write(script)

def get_restraints_with_kept_hydrogens(pdb_before_leap, pdb_after_leap):
    """allow the atoms added by tleap to be minimized"""

    lines_res = [[],[]]
    lines = [[],[]]
    for jdx, filename in enumerate([pdb_before_leap, pdb_after_leap]):
        resnum = 0
        with open(filename, 'r') as pdbf:
            for idx, line in enumerate(pdbf):
                if line.startswith('ATOM'):
                    if int(line[22:26]) != resnum: # new residue
                        if resnum != 0:
                            lines_res[jdx].append(curres)
                        resnum = int(line[22:26])
                        curres = [idx]
                    else:
                        curres.append(idx)
                lines[jdx].append(line)
            # last residue is the ligand

    lines_added = [] 
    for idx in range(len(lines_res[1])):
        for kdx in lines_res[1][idx]:
            is_added_after_tleap = True
            for jdx in lines_res[0][idx]:
                #print lines[0][jdx], lines[1][kdx]
                if lines[0][jdx][30:54] == lines[1][kdx][30:54]:
                    is_added_after_tleap = False
            #print is_added_after_tleap
            if is_added_after_tleap:
                lines_added.append(lines[1][kdx])
  
    print "Number of atoms added by tleap: %s"%len(lines_added)
    constraints_line = ' | '.join([':%s@%s'%(line[22:26].strip(),line[12:16].strip()) for line in lines_added])
    #return constraints_line

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

def create_pdbfile_with_restraints(file_rl, file_rst, force=50.0):

    with open(file_rl, 'r') as startfile:
         with open(file_rst, 'w') as rstf:
            for line in startfile:
                if line.startswith(('ATOM', 'HETATM')):
                    atomname = line[12:16].strip()
                    resname = line[17:20].strip()
                    if resname not in ['WAT', 'LIG'] and atomname in ['C', 'CA', 'N', 'O']:
                        newline = line[0:30] + '%8.3f'%force + line[38:]
                    else:
                        newline = line[0:30] + '%8.3f'%0.0 + line[38:]
                else:
                    newline = line
                rstf.write(newline)

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
        prepare_leap_config_file('leap.in', file_r, 'lig.mol2', 'complex.pdb')
        for idx, file_l in enumerate(files_l):
            shutil.copyfile(file_l, 'lig.mol2')
            # change ligand name to LIG
            ligname = reader.open('lig.mol2').ligname
            mol2.update_mol2file('lig.mol2', 'lig.mol2', ligname='LIG')

            # prepare ligand
            # since the charges of the original .mol2 file are eventually kept using pdb2mol2
            # do not regenerate the charges using antechamber when preparing the ligand
            prepare_ligand(file_r, 'lig.mol2', 'complex.pdb', charge_method='None')
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
        prepare_leap_config_file('leap.in', file_r, None, 'complex.pdb')
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
