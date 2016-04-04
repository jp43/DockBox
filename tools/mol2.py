import sys
import shutil

def get_atoms_names(filename):

    atoms_names = []
    with open(filename, 'r') as pdbfile:
        is_structure = False
        for line in pdbfile:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
            elif is_structure:
                line_s = line.split()
                atoms_names.append(line_s[1])
    return atoms_names

def give_unique_atom_names(file_l):

    tmpfile = 'tmp.mol2'
    with open(file_l, 'r') as oldf:
        newf = open(tmpfile, 'w')
        known_atom_types = []
        atom_numbers = []
        for line in oldf:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
                newf.write(line)
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
                newf.write(line)
            elif is_structure:
                line_s = line.split()
                atom = line_s[1]
                atom_type = ''.join([ch for ch in atom if not ch.isdigit()])
                if atom_type not in known_atom_types:
                    known_atom_types.append(atom_type)
                    atom_number = '1'
                    atom_numbers.append(1)
                else:
                    idx = known_atom_types.index(atom_type)
                    atom_number = str(atom_numbers[idx]+1)
                    atom_numbers[idx] += 1
                new_atom_type = line_s[1]+atom_number
                newline = ' '*(7-len(line_s[0])) + line_s[0] + 2 * ' ' + new_atom_type + \
' '*(10-len(new_atom_type)) + line[19:]
                newf.write(newline)
            else:
                newf.write(line)

    shutil.move(tmpfile, file_l)

def pdb2mol2(pdbfile, mol2file, sample=None):

    if not sample:
        raise NotImplementedError("No sample provided!")

    with open(mol2file, 'w') as mf:
        with open(pdbfile, 'r') as pdbf:
            with open(sample, 'r') as sf:
                # go to first atom line in sample file 
                line_s = sf.next()
                mf.write(line_s) # write non-atom lines 
                while not line_s.startswith('@<TRIPOS>ATOM'):
                    line_s = sf.next()
                    mf.write(line_s)
                line_s = sf.next()
                # go to first atom line in PDB file
                line_pdb = pdbf.next()
                while not line_pdb.startswith(('ATOM','HETATM')):
                    line_pdb = pdbf.next()
                # print 
                while line_pdb.startswith(('ATOM','HETATM')):
                    if line_s.split()[1] != line_pdb.split()[2]:
                        raise ValueError('mol2 sample provided does not match with pdbfile')
                    coords = line_pdb[30:54].split()
                    # write newline
                    newline = line_s[:16] + ' '*(10-len(coords[0])) + str(coords[0]) + \
' '*(10-len(coords[1])) + str(coords[1]) + ' '*(10-len(coords[2])) + str(coords[2]) + line_s[46:]
                    mf.write(newline)
                    try:
                        line_pdb = pdbf.next()
                        line_s = sf.next()
                    except StopIteration:
                        break
                for line_s in sf:
                    mf.write(line_s)
