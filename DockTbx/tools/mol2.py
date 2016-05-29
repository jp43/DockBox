import os
import sys
import shutil
import subprocess

from MolKit import Read
from PyBabel.atomTypes import AtomHybridization

from DockTbx.tools.babel import ArrangeHydrogens

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

    def write(self, filename, structs, mode='w', multi=False):

        if not isinstance(structs, list):
            structs = [structs]

        if not isinstance(filename, str):
            raise ValueError('filename should be a string!')

        if multi:
            suffix, ext = os.path.splitext(filename)
            filename = []
            for idx, struct in enumerate(structs):
                filename.append(suffix+str(idx+1)+ext)
        else:
            if len(structs) == 1:
                filename = [filename]
            else:
                raise ValueError(".mol2 writer not implemented to write a single file \
                    with multiple structures!")

        for idx, struct in enumerate(structs):
            with open(filename[idx], mode) as ff:
                for section in known_section:
                    if section in struct:
                        ff.write('@<TRIPOS>'+section+'\n')
                        for line in struct[section]:
                            if section == 'ATOM':
                                newline = '%7s %-5s    %9s %9s %9s %-5s    %2s %-5s     %8s\n'%tuple(line)
                            else:
                                newline = line
                            ff.write(newline)

def update_mol2file(inputfile, outputfile, ADupdate=None, multi=False, ligname=None, unique=False, mask=None):
    f = Reader(inputfile)
    structs = f.readlines()
    f.close()

    updated_structs = []
    for struct in structs:
        if ADupdate:
            struct = update_AD_output_from_original_struct(struct, ADupdate)
        if ligname:
            struct = update_ligand_name(struct, ligname)
        if unique:
            struct = give_unique_atom_names(struct, mask=mask)
        updated_structs.append(struct)

    Writer().write(outputfile, updated_structs, multi=multi)

def pdb2mol2(inputfile, outputfile, sample):
    # get atom lines in PDB:
    atom_lines_pdb = []
    with open(inputfile, 'r') as pdbf:
        for line in pdbf:
            if line.startswith(('ATOM','HETATM')):
                atom_lines_pdb.append(line)

    f = Reader(sample)
    structs = f.readlines()
    f.close()

    struct = structs[0]
    new_struct = struct
    for idx, line in enumerate(struct['ATOM']):
        atom_name_mol2 = line[1].lower()
        is_atom = False
        for line_pdb in atom_lines_pdb:
            atom_name_pdb = line_pdb[12:16].strip().lower()
            if atom_name_pdb == atom_name_mol2:
                if is_atom:
                    raise ValueError("Mol2 atom name already found in PDB file, your files should have unique atom names!")
                is_atom = True
                coords = [coord + '0' for coord in line_pdb[30:54].split()]
                for jdx in range(3):
                    new_struct['ATOM'][idx][jdx+2] = coords[jdx]
        if not is_atom:
            raise IOError("Mol2 atom name not found in PDB file, check your input files!")

    Writer().write(outputfile, new_struct)

def update_ligand_name(struct, ligname):
    new_struct = struct
    for idx, line in enumerate(struct['ATOM']):
        new_struct['ATOM'][idx][-2] = ligname

    ligname_p = new_struct['ATOM'][0][-2]
    if 'SUBSTRUCTURE' in struct:
        for idx, line in enumerate(struct['SUBSTRUCTURE']):
            new_struct['SUBSTRUCTURE'][idx] = line.replace(ligname_p, ligname)

    return new_struct

def is_unique_name(struct):
    known_atom_names = []
    for line in struct['ATOM']:
        atom_name = line[1]
        if atom_name not in known_atom_names:
            known_atom_names.append(atom_name)
        else:
            return False
    return True

def update_AD_output_from_original_struct(struct1, filename):
    # read original structure
    f = Reader(filename)
    struct2 = f.next()
    new_struct = struct2
    f.close()

    for idx, struct in enumerate([struct1, struct2]):
        if not is_unique_name(struct):
            raise ValueError("Mol2 structure (%i) should have unique atom names"%(idx+1))

    for idx, line2 in enumerate(struct2['ATOM']):
        for line1 in struct1['ATOM']:
            if line1[1] == line2[1]:
                for jdx in range(2,5):
                    new_struct['ATOM'][idx][jdx] = line1[jdx]

    return new_struct

def arrange_hydrogens(inputfile, outputfile):

    mol = Read(inputfile)
    base, ext = os.path.splitext(inputfile)

    # remove hydrogens from structure
    inputfile_noH = base + '_noH' + ext
    subprocess.check_output('babel -imol2 %s -omol2 %s -d &>/dev/null'%(inputfile,inputfile_noH), shell=True, executable='/bin/bash')
    molnoH = Read(inputfile_noH)

    allAtoms = mol.allAtoms
    allAtomsNoH = molnoH.allAtoms

    babel = AtomHybridization()

    babel.assignHybridization(allAtoms)
    babel.assignHybridization(allAtomsNoH)

    # get mol2 ids of all atoms and all hydrogens
    ff = Reader(inputfile)
    struct = ff.next()
    hat_mol2_ids = []
    at_mol2_ids = []
    for line in struct['ATOM']:
        atom_name = line[5]
        #print line[2]
        if atom_name[0] in ['H', 'h']:
            hat_mol2_ids.append(line[0])
        at_mol2_ids.append(line[0])

    hat_mol2_ids = map(int, hat_mol2_ids)
    at_mol2_ids = map(int, at_mol2_ids)

    # find out the heavy atom each hydrogen is bound with
    at_with_hat_mol2_ids = []
    for id in hat_mol2_ids:
        for line in struct['BOND']:
            line_s = line.split()
            origin_atom_id = int(line_s[1])
            target_atom_id = int(line_s[2])
            #print origin_atom_id
            if id == origin_atom_id:
                at_with_hat_mol2_ids.append(target_atom_id)
            elif id == target_atom_id:
                at_with_hat_mol2_ids.append(origin_atom_id)

    if len(at_with_hat_mol2_ids) != len(hat_mol2_ids):
        raise ValueError("Each hydrogen should have only one bound! Check you .mol2 file")

    #print hat_mol2_ids
    #print at_with_hat_mol2_ids

    addh = ArrangeHydrogens()
    # at_with_hat_idxs are the indices of atoms related to each hydrogen of hat
    hat, at_with_hat_idxs = addh.addHydrogens(allAtoms, allAtomsNoH)

    hat_coords = []
    hat_done_mol2_ids = []
    for idx, at_with_hat_idx in enumerate(at_with_hat_idxs):
        at_with_hat_mol2_id = at_mol2_ids[at_with_hat_idx]
        hat_mol2_ids_cur = [hat_mol2_ids[jdx] for jdx, id in enumerate(at_with_hat_mol2_ids) \
            if id == at_with_hat_mol2_id]
        kdx = 0
        while hat_mol2_ids_cur[kdx] in hat_done_mol2_ids:
            kdx += 1
        hat_done_mol2_ids.append(hat_mol2_ids_cur[kdx])
        hat_coords.append(hat[idx][0])

    for line in struct['ATOM']:
        id = int(line[0])
        if id in hat_done_mol2_ids:
            idx = hat_done_mol2_ids.index(id)
            line[2:5] = map("{:.4f}".format, hat_coords[idx])

    Writer().write(outputfile, struct)
    os.remove(inputfile_noH)

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

