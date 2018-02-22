import os
import sys
import shutil
import subprocess

from MOBPred.tools import reader

def get_nwaters(logfile):

    with open(logfile, 'r') as logf:
        for line in logf:
            line_s = line.split()
            if len(line_s) == 3 and line_s[0] == 'Added' and line_s[-1] == 'residues.':
                return int(line_s[1])

def get_removed_waters(file_r, files_l, file_c, nwaters_tgt, boxsize, step=0.01, ntries=5):

    if not files_l:
        file_c = file_r

    # determine the number of solute residues
    nsolutes = int(subprocess.check_output("echo `cpptraj -p %(file_c)s -mr '*' | awk '{print NF - 1;}'`"%locals(), shell=True))

    print "Targeted number of water residues:", nwaters_tgt
    print "Number of residues found for solute:", nsolutes

    distance = [boxsize]
    for idx in range(ntries):
        distance.append(boxsize + (idx+1)*step)
        distance.append(boxsize - (idx+1)*step)

    diff_best = 1e10
    for d in distance:
        c = 1.0
        lastdiff = 1e10
        while True:
            prepare_leap_config_file('leap.in', file_r, files_l, file_c, solvate=True, distance=d, closeness=c)
            subprocess.check_output('tleap -f leap.in > leap.log', shell=True)
            nwaters = get_nwaters('leap.log')
            diff = nwaters - nwaters_tgt
            #print d, c, diff, nwaters
            if diff > 0:
                if lastdiff < 0:
                    if diff < diff_best:
                        diff_best = diff
                        dbest = d
                        cbest = c
                        nwaters_best = nwaters
                    break
                else:
                   lastdiff = diff
                c += 0.01
            elif diff < 0:
                c -= 0.01
                lastdiff = diff
            elif diff == 0:
                diff_best = diff
                dbest = d
                cbest = c
                nwaters_best = nwaters
                break
        if diff == 0:
            break

    print "Closest number of water residues found:", nwaters_best
    print "Removing %i water residues..."%diff_best

    removed_waters = [nsolutes + nwaters_best - diff_best + idx + 1 for idx in range(diff_best)]
    return removed_waters, dbest, cbest

def get_solvent_mask(pdbfile, residues='WAT'):

    solvent_residues = map(str.strip, residues.split(','))

    # set residue numbers
    resnum = None
    with open(pdbfile) as pdbf:
        for line in pdbf:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                if resname in solvent_residues:
                    if not resnum:
                        resnum_in = line[22:27].strip()
                    resnum = line[22:27].strip()

    resnum_fin = resnum
    solvent_mask = ':%s-%s'%(resnum_in, resnum_fin)

    return solvent_mask

def get_ions_number(logfile, concentration=0.15):

    # Nions = Cions * Nwater * 1/55.5 where 55.5 M is the concentration of pure water

    with open(logfile, 'r') as lf:
        for line in lf:
            line_s = line.strip().split()
            if len(line_s) > 2:
                if line_s[0] == 'Added' and line_s[2] == 'residues.':
                    nresidues = int(line_s[1])
                    ncl = int(round(nresidues * concentration / 55.5))
                    nna = ncl
            if line.startswith("WARNING: The unperturbed charge"):
                net_charge = int(round(float(line_s[7])))
                if net_charge > 0:
                    ncl += abs(net_charge)
                elif net_charge < 0:
                    nna += abs(net_charge)
    return nna, ncl

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

    info['ACE'] = ['CH3', 'C', 'O', 'HT1', 'HT2', 'HT3']
    info['NME'] = ['CH3', 'C', 'N', 'HT1', 'HT2', 'HT3', 'H']
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


def get_ligand_name(filename):

    # finding ligand name...
    prefix, ext = os.path.splitext(filename)

    r_ligand = reader.open(filename)
    struct = r_ligand.next()
    if ext == '.pdb':
        idxname = 2
    elif ext == '.mol2':
        idxname = 7
    else:
        raise ValueError('File format for file_l should be .pdb or .mol2!')
    names = [row[idxname] for row in struct['ATOM']]
    name = list(set(names))

    if len(name) != 1:
        raise IOError('More than one ligand found in file %s'%file_l)
    else:
        name = name[0]
    return name

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

def run_antechamber(infile, outfile, at='gaff', c='gas', logfile='antechamber.log'):
    """ use H++ idea of running antechamber multiple times with bcc's 
charge method to estimate the appropriate net charge!!"""

    suffix, ext = os.path.splitext(infile)
    ext = ext[1:]

    max_net_charge = 30
    net_charge = [0]

    if c.lower() != 'none':
        for nc in range(max_net_charge):
            net_charge.extend([nc+1,-(nc+1)])

        for nc in net_charge:
            iserror = False
            command = 'antechamber -i %(infile)s -fi %(ext)s -o %(outfile)s -fo mol2 -at %(at)s -c %(c)s -nc %(nc)s -du y -pf y >> %(logfile)s'%locals()
            subprocess.check_call('echo "# command used: %(command)s" > %(logfile)s'%locals(), shell=True) # print command in logfile
            subprocess.check_call(command, shell=True)
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
       command = 'antechamber -i %(infile)s -fi %(ext)s -o %(outfile)s -fo mol2 -at %(at)s -du y -pf y > %(logfile)s'%locals()
       subprocess.check_output(command, shell=True)

def prepare_leap_config_file(script_name, file_r, files_l, file_rl, solvate=False, PBRadii=None, forcefield='leaprc.ff14SB', nna=0, ncl=0, distance=10.0, closeness=1.0, remove=None, model='TIP3P'):
 
    solvation_line = ""
    pbradii_lines = ""
    add_ions_lines = ""
    remove_water_lines = ""
    ligand_lines = ""
    ions_libraries_lines = ""

    tip3p_models = ['TIP3P', 'TIP3PF', 'POL3', 'QSPCFW']
    tip4p_models = ['TIP4P', 'TIP4PEW']
    spc_models = ['QSPCFW', 'SPC', 'SPCFW']

    if solvate:
        boxtype = model.upper() + 'BOX'
        solvation_line = "\nsolvateBox complex %s %.2f %.2f"%(boxtype,distance,closeness)

        if nna > 0:
            add_ions_lines += "\naddions complex Na+ %i"%nna
        if ncl > 0:
            add_ions_lines += "\naddions complex Cl- %i"%ncl

        if model.upper() in tip3p_models:
            suffix_ions_libraries = 'tip3p'
        elif model.upper() in tip4p_models:
            suffix_ions_libraries = 'tip4pew'
        elif model.upper() in spc_models:
            suffix_ions_libraries = 'spce'
        else:
            raise ValueError('Solvent model %s unknown, should be one of !'% \
(model,', '.joined(tip3p_models + tip4p_models + spc_models)))
    else:
        suffix_ions_libraries = 'tip3p'

    if PBRadii:
        pbradii_lines = "\nset default PBRadii %s"%PBRadii

    if remove:
        for idx in remove:
            remove_water_lines += "\nremove complex complex.%i"%idx

    # loadoff atomic_ions.lib
    if files_l:
        if isinstance(files_l, basestring):
            files_l = [files_l]
        for idx, file_l in enumerate(files_l):
            file_l_prefix, ext = os.path.splitext(file_l)
            file_l_prefix = os.path.basename(file_l_prefix)
            name = get_ligand_name(file_l)
            ligand_lines += "\n%(name)s = loadmol2 %(file_l)s\nloadamberparams %(file_l_prefix)s.frcmod"%locals()
        with open(script_name, 'w') as leapf:
            script ="""source %(forcefield)s
source leaprc.gaff
loadamberparams frcmod.ionsjc_%(suffix_ions_libraries)s
loadamberparams frcmod.ionslm_1264_%(suffix_ions_libraries)s%(ligand_lines)s%(pbradii_lines)s
complex = loadPdb %(file_rl)s%(solvation_line)s%(add_ions_lines)s%(remove_water_lines)s
saveAmberParm complex start.prmtop start.inpcrd
savePdb complex start.pdb
quit\n"""%locals()
            leapf.write(script)
    else:
        with open(script_name, 'w') as leapf:
            script ="""source %(forcefield)s
loadoff atomic_ions.lib
loadamberparams frcmod.ionsjc_%(suffix_ions_libraries)s
loadamberparams frcmod.ionslm_1264_%(suffix_ions_libraries)s%(pbradii_lines)s
complex = loadPdb %(file_r)s%(solvation_line)s%(add_ions_lines)s
saveAmberParm complex start.prmtop start.inpcrd
savePdb complex start.pdb
quit\n"""%locals()
            leapf.write(script)

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

def prepare_ligand(file_r, files_l, file_rl, charge_method='gas'):

    if isinstance(files_l, basestring):
        files_l = [files_l]

    mol2files_l = []

    shutil.copyfile(file_r, file_rl)
    for file_l in files_l:
        file_l_prefix, ext = os.path.splitext(file_l)
        file_l_prefix = os.path.basename(file_l_prefix)

        mol2file = file_l_prefix + '.mol2'
        run_antechamber(file_l, 'tmp.mol2', at='gaff', c=charge_method)

        shutil.move('tmp.mol2', mol2file)
        subprocess.check_output('parmchk -i %s -f mol2 -o %s.frcmod'%(mol2file, file_l_prefix), shell=True, executable='/bin/bash')
        subprocess.check_output('antechamber -i %s -fi mol2 -o %s.pdb -fo pdb > /dev/null'%(mol2file, file_l_prefix), shell=True, executable='/bin/bash')

        mol2files_l.append(mol2file)
        with open(file_rl, 'a') as ffrl:
            with open(file_l_prefix + '.pdb', 'r') as ffl:
                for line in ffl:
                    if line.startswith(('ATOM', 'HETATM')):
                        ffrl.write(line)
            ffrl.write('TER\n')

    return mol2files_l
