import fileinput
import shutil
import subprocess

def rearrange_atom_names(file_l, order):

    tmpfile = 'tmp.pdb'

    with open(file_l, 'r') as oldf:
        with open(tmpfile, 'w') as newf:
            line = oldf.next()
            atom_lines = []
            while not line.startswith(('ATOM', 'HETATM')):
                line = oldf.next()
            while line.startswith(('ATOM', 'HETATM')):
                atom_lines.append(line)
                try:
                    line = oldf.next()
                except StopIteration:
                    break
            new_idxs = []
            for line in atom_lines:
                atom_name = str(line[12:17]).strip()
                new_idxs.append(order.index(atom_name))

            new_atom_lines = []
            for idx in range(len(order)):
                new_atom_lines.append(None)

            for idx, line in enumerate(atom_lines):
                new_atom_lines[new_idxs[idx]] = line

            new_atom_lines = filter(lambda a: a, new_atom_lines)

            idx = 0
            for line in new_atom_lines:
                idx += 1
                newline = line[:6] + ' '*(5-len(str(idx))) + str(idx) + line[11:]
                newf.write(newline)

    shutil.move(tmpfile, file_l)

def give_unique_atom_names(file_l):

    tmpfile = 'tmp.pdb'
    with open(file_l, 'r') as oldf:
        newf = open(tmpfile, 'w')
        known_atom_types = []
        atom_numbers = []
        for line in oldf:
            if line.startswith(('ATOM', 'HETATM')):
                atom_type = line[13:14].strip()
                if atom_type not in known_atom_types:
                    known_atom_types.append(atom_type)
                    atom_number = '1'
                    atom_numbers.append(1)
                else:
                    idx = known_atom_types.index(atom_type)
                    atom_number = str(atom_numbers[idx]+1)
                    atom_numbers[idx] += 1
                newline = line[:14]+atom_number+(3-len(atom_number))*' '+line[17:]
                # write line
                newf.write(newline)
            elif line.startswith(('ENDMDL','END')):
                newf.write(line)
                known_atom_types = []
                atom_numbers = []
            else:
                newf.write(line)
        newf.close()

    shutil.move(tmpfile, file_l)

def remove_hydrogens(infile, outfile):

    atom_number = 0
    with open(infile, 'r') as oldf:
        newf = open(outfile, 'w')
        for line in oldf:
            if line.startswith(('ATOM','HETATM')):
                atom_type = line[13:14]
                if not atom_type == 'H':
                    if atom_number == 99999:
                        atom_number = 1
                    else:
                        atom_number += 1
                    newf.write(line[:6] + ' '*(5-len(str(atom_number))) + str(atom_number) + line[11:])
            else:
                newf.write(line)

        newf.close()

def format_lig_file(file_l):

    tmpfile = 'tmp.pdb'
    is_atom_on_last_line = False
    with open(file_l, 'r') as oldf:
        newf = open(tmpfile, 'w')
        for line in oldf:
            if line.startswith(('ATOM','HETATM','CONECT')):
                newline = line.replace('\n','')
                for resname in ['UNK','UNL']:
                    newline = newline.replace(resname,'LIG')
                newline = newline.replace('Cl','CL') 
                newline = newline.replace('Br','BR')
                newline = newline.replace('HETATM','ATOM  ')
                print >> newf, newline
                if not is_atom_on_last_line:
                    is_atom_on_last_line = True
            else:
                if is_atom_on_last_line:
                    print >> newf, 'ENDMDL'
                is_atom_on_last_line = False
        newf.close()

    shutil.move(tmpfile, file_l)

def create_reclig_file(file_r, file_l, file_rl):
    """Create a PDB file containing both receptor and ligand"""

    nends = [0, 0]
    # (A) check the number of ligand and receptors in output files
    for idx, filename in enumerate([file_r, file_l]):
        is_atom_on_last_line = False
        with open(filename, 'r') as ff:
           for line in ff:
               if line.startswith(('ATOM','HETATM','TER')):
                   if not is_atom_on_last_line:
                       nends[idx] += 1
                       is_atom_on_last_line = True
               else:
                       is_atom_on_last_line = False

    # (B) check if there is a unique receptor for multiple ligands (docking with fix receptor)
    if nends[0] == 1 and nends[1] > 1:
        unique_receptor = True
    elif nends[0] == nends[1]:
        unique_receptor = False
    else:
        raise ValueError('the numbers of ligands and receptors are not the same')

    # (C) build .pdb file with ligand + receptor
    with open(file_rl, 'w') as reclf:
        if unique_receptor:
            with open(file_l, 'r') as ligf:
                for idx in range(nends[1]):
                    reclf.write('MODEL     %s\n'%(idx+1))
                    # copy receptor structure
                    with open(file_r, 'r') as recf:
                        is_ter_on_last_line = False
                        for line in recf:
                            if line.startswith(('ATOM', 'HETATM')):
                                reclf.write(line)
                                is_ter_on_last_line = False
                            elif line.startswith('TER'):
                                reclf.write(line)
                                is_ter_on_last_line = True
                        if not is_ter_on_last_line:
                            reclf.write('TER\n')
                    # copy ligand structure
                    is_atom_on_last_line = False
                    while True:
                        line = ligf.next()
                        if line.startswith(('ATOM', 'HETATM')):
                            reclf.write(line)
                            is_atom_on_last_line = True
                        else:
                            if is_atom_on_last_line:
                                reclf.write('ENDMDL\n')
                                break
        else:
            with open(file_r, 'r') as recf:
                with open(file_l, 'r') as ligf:
                    for idx in range(nends[0]):
                        reclf.write('MODEL     %s\n'%(idx+1))
                        is_atom_on_last_line = False
                        is_ter_on_last_line = False
                        # copy receptor structure
                        while True:
                            line = recf.next()
                            if line.startswith(('ATOM', 'HETATM')):
                                reclf.write(line)
                                is_atom_on_last_line = True
                                is_ter_on_last_line = False
                            elif line.startswith('TER'):
                                reclf.write(line)
                                is_atom_on_last_line = False
                                is_ter_on_last_line = True
                            else:
                                if is_atom_on_last_line:
                                    reclf.write('TER\n')
                                    break
                                elif is_ter_on_last_line:
                                    break
                        is_atom_on_last_line = False
                        while True:
                            line = ligf.next()
                            if line.startswith(('ATOM', 'HETATM')):
                                reclf.write(line)
                                is_atom_on_last_line = True
                            else:
                                if is_atom_on_last_line:
                                    reclf.write('ENDMDL\n')
                                    break

# the following lines can be used in case we want unique atom numbers and residue numbers
# when the file .pdb containing the complex target+ligand is created. In this case, the
# atom numbers and residue numbers of the ligand won't start from 1 but from where the atom numbers
# and residue numbers of the receptor ended (not useful! See Khaled's files target-ligand.pdb).

# copy ligand structure
#line = ''
#while not line.startswith('END'):
#    line = ligf.next()
#    if line.startswith(('ATOM', 'HETATM')):
#        atomidx = str(int(line[6:11]) + natoms[0])
#        residx = str(int(line[22:26]) + nresidues[0])
#        reclf.write(line[:6] + ' '*(5-len(atomidx)) + atomidx + line[11:22]  + ' '*(4-len(residx)) + residx + line[26:])
