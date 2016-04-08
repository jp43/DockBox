import sys
import shutil

def get_atoms_names(filename):

    atoms_names = []
    with open(filename, 'r') as mol2f:
        is_structure = False
        for line in mol2f:
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

def get_coordinates(filename):

    coords = []
    with open(filename, 'r') as mol2f:
        is_structure = False
        for line in filename:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
            elif is_structure:
                line_s = line.split()
                coords.append(map(float,line_s[2:5]))
    return coords

def update_mol2_from_pdb(input_pdbfile, output_mol2file, sample_mol2file=None):
    """update_mol2_from_pdb(input_pdbfile, output_mol2file, sample_mol2file=None)

        Update the sample .mol2 file provided with the coordinates of the atoms 
from the pdbfile whose atom names match with the ones in the sample

!!!!!!WARNING: the function assumes unique atom names in both the pdb and mol2 files"""

    if not sample_mol2file:
        raise NotImplementedError("No .mol2 sample provided!")

    # load atom lines in PDB:
    atom_lines_pdb = []
    with open(input_pdbfile, 'r') as pdbf:
        for line in pdbf:
            if line.startswith(('ATOM','HETATM')):
                atom_lines_pdb.append(line)

    # load atom lines in mol2:
    atom_lines_mol2 = []
    with open(sample_mol2file, 'r') as mol2f:
        is_structure = False
        for line in mol2f:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
            elif is_structure:
                atom_lines_mol2.append(line)

    new_atom_lines_mol2 = []
    for line_mol2 in atom_lines_mol2:
        newline = line_mol2
        atom_name_mol2 = line_mol2.split()[1]        
        for line_pdb in atom_lines_pdb:
            atom_name_pdb = line_pdb[12:16].strip()
            if atom_name_pdb == atom_name_mol2:
                coords = [coord + '0' for coord in line_pdb[30:54].split()]
                newline = line_mol2[:16] + ' '*(10-len(coords[0])) + str(coords[0]) + \
' '*(10-len(coords[1])) + str(coords[1]) + ' '*(10-len(coords[2])) + str(coords[2]) + line_mol2[46:] 
        new_atom_lines_mol2.append(newline)

    idx = 0
    with open(output_mol2file, 'w') as mf:
        with open(sample_mol2file, 'r') as mol2f:
            is_structure = False
            for line in mol2f:
                if line.startswith('@<TRIPOS>ATOM'):
                    is_structure = True
                    mf.write(line)
                elif line.startswith('@<TRIPOS>'):
                    is_structure = False
                    mf.write(line)
                elif is_structure:
                    mf.write(new_atom_lines_mol2[idx])
                    idx += 1
                else:
                    mf.write(line)
