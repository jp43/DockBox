import os
import sys
import shutil

# recognized sections
known_section = ['MOLECULE', 'ATOM', 'BOND', 'SUBSTRUCTURE']

class Reader(object):

    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, 'r')

        # move to first molecule
        for line in self.file:
            if line.startswith('@<TRIPOS>MOLECULE'):
                break 

    @property
    def ligname(self):
        ff = open(self.filename, 'r')
        first_atom = False
        ligname = None
        for line in ff:
            if line.startswith('@<TRIPOS>ATOM'):
                first_atom = True
            elif first_atom:
                line_s = line.split()
                if not ligname:
                    ligname = line_s[-2]
                elif line_s[-2] != ligname:
                    raise ValueError('Ligand name not consistent between structures')
                first_atom = False
        ff.close()
        return ligname                

    @property
    def nmolecules(self):
        ff = open(self.filename, 'r')
        nmolecules = 0
        for line in ff:
            if line.startswith('@<TRIPOS>MOLECULE'):
                nmolecules += 1
        ff.close()
        return nmolecules

    def read(self):
        struct = self.next()
        while struct is not None:
            yield struct
            struct = self.next()

    def readlines(self):
        structs = []
        struct = self.next()
        while struct is not None:
            structs.append(struct)
            struct = self.next()
        return structs

    def next(self):
        struct = None
        for idx, line in enumerate(self.file):
            # initialize stucture
            if idx == 0:
                struct = {}
                section = 'MOLECULE'
                struct[section] = []
            # new molecule detected
            if line.startswith('@<TRIPOS>MOLECULE'):
                break
            elif line.startswith('@<TRIPOS>'):
                section = line[9:].strip()
                struct[section] = []
            elif section and line.strip():
                if section == 'ATOM':
                    atom = line.split()
                    struct[section].append(atom)
                else:
                    struct[section].append(line)
            elif section and not line.strip():
                struct[section].append(line)
                section = None
        return struct

    def close(self):
        self.file.close()

    def __iter__(self):
        return self.read()

    readline = next

class Writer(object):

    def write(self, filename, contents, mode='w', ligname=None, unique=False, mask=None):

        filename = self.get_filename(filename, contents)
        if not isinstance(contents, list):
            contents = [contents]

        for idx, struct in enumerate(contents):
            if unique:
                struct = give_unique_atom_names(struct, mask=mask)
            with open(filename[idx], mode) as ff:
                for section in known_section:
                    if section in struct:
                        ff.write('@<TRIPOS>'+section+'\n')
                        for line in struct[section]:
                            if section == 'ATOM':
                                if ligname:
                                    line[-2] = ligname
                                newline = '%7s %-5s    %9s %9s %9s %-5s    %2s %-5s     %8s\n'%tuple(line)
                            else:
                                newline = line
                            ff.write(newline)

    def get_filename(self, filename, contents):

        new_filename = []
        if isinstance(contents, list):
            suffix, ext = os.path.splitext(filename)
            for idx, struct in enumerate(contents):
                new_filename.append(suffix+str(idx+1)+ext)
        else:
            new_filename.append(filename)
        return new_filename


def give_unique_atom_names(struct, mask=None):

    known_atom_names = []
    atom_numbers = []

    for line in struct['ATOM']:
        if not mask or line[-4] in mask:
            atom = line[1]
            atom_name = ''.join([ch for ch in atom if not ch.isdigit()])
            if atom_name not in known_atom_names:
                known_atom_names.append(atom_name)
                atom_number = '1'
                atom_numbers.append(1)
            else:
                idx = known_atom_names.index(atom_name)
                atom_number = str(atom_numbers[idx]+1)
                atom_numbers[idx] += 1

    # generate new atom names with 3 characters
    new_atom_names = []
    natoms = len(known_atom_names)
    for idx in range(natoms):
        name = known_atom_names[idx] + str(atom_numbers[idx])
        nchars = len(name)
        if nchars > 5:
            raise ValueError("atom with more than 5 characters detected!")
        if nchars > 3 and nchars <= 5:
            nletters = len(known_atom_names[idx])
            nfigs = len(str(atom_numbers[idx]))
            if nfigs >= 3:
                raise ValueError("more than 99 atoms of the same type detected!")
            elif nletters == 2:
                new_atom_names.append(known_atom_names[idx][:1])
            elif nletters == 3:
                new_atom_names.append(known_atom_names[idx][:3-nfigs])
        elif nchars <= 3:
            new_atom_names.append(known_atom_names[idx])

    new_struct = struct

    new_known_atom_names = []
    new_atom_numbers = []
    for jdx, line in enumerate(struct['ATOM']):
        if not mask or line[-4] in mask:
            atom = line[1]
            atom_name = ''.join([ch for ch in atom if not ch.isdigit()])
            idx = known_atom_names.index(atom_name)
            atom_name = new_atom_names[idx]
            if atom_name not in new_known_atom_names:
                new_known_atom_names.append(atom_name)
                atom_number = '1'
                new_atom_numbers.append(1)
            else:
                idx = new_known_atom_names.index(atom_name)
                atom_number = str(new_atom_numbers[idx]+1)
                new_atom_numbers[idx] += 1
            new_struct['ATOM'][jdx][1] = atom_name+atom_number

    return new_struct

def get_ligand_name(filename):

    with open(filename, 'r') as mol2f:
        is_structure = False
        for line in mol2f:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
            elif is_structure:
                line = line.split()[7]
                break
    return line[0:3]

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

def change_ligand_name(file_l, newligname):

    tmpfile = 'tmp.mol2'
    with open(file_l, 'r') as oldf:
        newf = open(tmpfile, 'w')

        known_atom_types = []
        atom_numbers = []
        is_first_atom = True

        for line in oldf:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
                newf.write(line)
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
                newf.write(line)
            elif is_structure:
                line_s = line.rsplit(None, 2)
                ligname = line_s[-2]
                newline = line.replace(ligname,newligname)
                newf.write(newline)
            else:
                newf.write(line)

    shutil.move(tmpfile, file_l)

def get_coordinates(filename):

    coords = []
    with open(filename, 'r') as mol2f:
        is_structure = False
        for line in mol2f:
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
from the pdbfile whose atom names match with the ones in the sample (the function 
assumes unique atom names in both the pdb and mol2 files)"""

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
        line_mol2_start = line_mol2.rsplit(None, 7)[0]
        line_mol2_end = line_mol2.split(None, 5)[5]
        atom_name_mol2 = line_mol2.split()[1].lower()
        for line_pdb in atom_lines_pdb:
            atom_name_pdb = line_pdb[12:16].strip().lower()
            if atom_name_pdb == atom_name_mol2:
                coords = [coord + '0' for coord in line_pdb[30:54].split()]
                newline = line_mol2_start + ' '*(10-len(atom_name_mol2)) + ' '*(10-len(coords[0])) + str(coords[0]) + \
' '*(10-len(coords[1])) + str(coords[1]) + ' '*(10-len(coords[2])) + str(coords[2]) + ' ' + line_mol2_end
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
