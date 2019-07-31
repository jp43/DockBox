import os
import sys
import math
import numpy as np
import nwalign as nw

from mdkit.utility import mol2
from dockbox import pyqcprot

# prefix to identify ligand, target and isomer directories
ligdir_prefix = 'lig'
tardir_prefix = 'target'
isodir_prefix = 'isomer'

residues_3_to_1 = {'ALA': 'A',
'ARG': 'R',
'ASN': 'N',
'ASP': 'D',
'CYS': 'C',
'GLU': 'E',
'GLY': 'G',
'HIS': 'H',
'ILE': 'I',
'LEU': 'L',
'LYS': 'K',
'MET': 'M',
'PHE': 'F',
'PRO': 'P',
'GLN': 'Q',
'SER': 'S',
'SEC': 'U',
'THR': 'T',
'TRP': 'W',
'TYR': 'Y',
'VAL': 'V'}

equivalent_residues = {'CYM': 'CYS',
'LYN': 'LYS',
'ASH': 'ASP',
'CYX': 'CYS',
'GLH': 'GLU',
'HID': 'HIS',
'HIE': 'HIS',
'HIP': 'HIS'}

def get_total_residue_number(filename):
    indices = []
    nresidues = 0
    with open(filename, 'r') as pdbf:
        for line in pdbf:
            if line.startswith('ATOM'):
                resnum = line[22:26].strip()
                if resnum not in indices:
                    indices.append(resnum)
                    nresidues += 1
    return nresidues

def get_sequence_from_PDB(filename):
    indices = []
    sequence = ''

    with open(filename, 'r') as pdbf:
        for line in pdbf:
            if line.startswith('ATOM'):
                resnum = line[22:26].strip()
                resname = line[17:20].strip()
                if resname in equivalent_residues:
                    resname = equivalent_residues[resname]

                if resnum not in indices and resname in residues_3_to_1:
                    sequence += residues_3_to_1[resname]
                    indices.append(resnum)
    return sequence, indices

def get_residues_coordinates(filename, indices):
    indices_new = []
    coords = []

    with open(filename, 'r') as pdbf:
        for line in pdbf:
            if line.startswith('ATOM'):
                resnum = line[22:26].strip()
                resname = line[17:20].strip()
                atomname = line[12:16].strip()

                if resnum not in indices_new and resnum in indices:
                    coords.append([])
                    indices_new.append(resnum)
                #if resnum in indices and atomname[0] != 'H':
                if resnum in indices and atomname in ['CA', 'C', 'N', 'O']:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[-1].append([atomname, x, y, z])

    return coords, indices_new

def compute_rmsd(file1, file2, rotmat=np.eye(3), trans1=np.zeros(3), trans2=np.zeros(3)):
    """Compute RMSD between 2 poses"""

    for file in [file1, file2]:
        if isinstance(file, float) and math.isnan(file):
            return float('nan')

    # load coordinates of first pose (non-hydrogen atoms)
    coords1 = mol2.get_coordinates(file1, keep_h=False)
    coords1 = np.array(coords1)
    natoms = coords1.shape[0]

    coords1_rot = np.empty_like(coords1)
    for idx in range(natoms):
        coords1t = coords1[idx,:] + trans1
        coords1t = coords1t[:,np.newaxis]
        coords1_rot[idx,:] = np.dot(rotmat, coords1t).flatten() - trans2

    # load coordinates of second pose (non-hydrogen atoms)
    coords2 = mol2.get_coordinates(file2, keep_h=False)
    coords2 = np.array(coords2)

    rmsd = np.sqrt(np.sum((coords1_rot-coords2)**2)/natoms)
    return rmsd

def get_rmsd_rotation_and_translations(file1, file2):

    nres1 = get_total_residue_number(file1)
    nres2 = get_total_residue_number(file2)

    seq1, ind1 = get_sequence_from_PDB(file1)
    seq2, ind2 = get_sequence_from_PDB(file2)

    alignment = nw.global_align(seq1, seq2)

    nalign = len(alignment[0])
    nresidues_min = min(len(seq1), len(seq2))

    ind1new = []
    ind2new = []
    idx1, idx2 = 0, 0

    for idx in range(nalign):
        if (idx < nresidues_min) and seq1[idx] == seq2[idx] and seq1[idx] != '-':
            ind1new.append(ind1[idx1])
            ind2new.append(ind2[idx2])
        if (idx < len(seq1)) and seq1[idx] != '-':
            idx1 += 1
        if (idx < len(seq2)) and seq2[idx] != '-':
            idx2 += 1

    ind1 = ind1new
    ind2 = ind2new

    #TODO: add a threshold for the number of residues considered
    frac1 = len(ind1)*100.0/nres1
    frac2 = len(ind2)*100.0/nres2

    # get coordinates of specific residues 
    coords1, ind1 = get_residues_coordinates(file1, ind1)
    coords2, ind2 = get_residues_coordinates(file2, ind2)

    new_coords1 = []
    new_coords2 = []

    # check if there is consistency in atom names
    nresidues1 = len(coords1)
    for idx in range(nresidues1):
        coords1_res = coords1[idx]
        coords2_res = coords2[idx]

        atomnames1 = [item[0] for item in coords1_res]
        atomnames2 = [item[0] for item in coords2_res]
        if set(atomnames1) != set(atomnames2):
            sys.exit("Inconsistency found in residue %s in file %s and residue %s in file %s! Missing atom suspected..."%(ind1[idx],file1,ind2[idx],file2))

        # create new coordinates
        for an1, x1, y1, z1 in coords1_res:
            for an2, x2, y2, z2 in coords2_res:
                if an1 == an2:
                    new_coords1.append([x1, y1, z1])
                    new_coords2.append([x2, y2, z2])
                    break

    new_coords1 = np.array(new_coords1).T
    new_coords2 = np.array(new_coords2).T

    rotation = np.zeros(9)
    trans1 = -new_coords1[:,0]
    trans2 = -new_coords2[:,0]

    rmsd = pyqcprot.CalcRMSDRotationalMatrix(new_coords1, new_coords2, rotation, None)

    rotation = rotation.reshape((3, 3))
    trans1 += new_coords1[:,0]
    trans2 += new_coords2[:,0]

    return rotation, trans1, trans2

def get_rmsd_rotation_and_translations_all_targets(files_r):
    rmsd_rot_trans = {}

    for key1 in files_r:
        rmsd_rot_trans[key1] = {}

        for key2 in files_r:
            if key1 == key2:
                rotation = np.eye(3)
                trans1 = np.zeros(3)
                trans2 = np.zeros(3)
            else:
                file1 = files_r[key1]
                file2 = files_r[key2]
                rotation, trans1, trans2 = get_rmsd_rotation_and_translations(file1, file2)
            rmsd_rot_trans[key1][key2] = [rotation, trans1, trans2]

    return rmsd_rot_trans

def check_architecture(directory):
    """Check architecture %s*/%s*/%s* of specified directories"""%(ligdir_prefix,tardir_prefix,isodir_prefix)

    if os.path.isdir(directory):
        dir_split = directory.split('/')
        if dir_split[-1].startswith(isodir_prefix):
            isisomerID = True
            if len(dir_split) > 1 and dir_split[-2].startswith(tardir_prefix):
                istargetID = True
                if len(dir_split) > 2 and dir_split[-3].startswith(ligdir_prefix):
                    isligID = True
                else:
                    isligID = False
            elif len(dir_split) > 1 and dir_split[-2].startswith(ligdir_prefix):
                istargetID = False
                isligID = True
            else:
                istargetID = False
                isligID = False
        elif dir_split[-1].startswith(tardir_prefix):
            isisomerID = False
            istargetID = True
            if len(dir_split) > 1 and dir_split[-2].startswith(ligdir_prefix):
                isligID = True
            else:
                isligID = False
        elif dir_split[-1].startswith(ligdir_prefix):
            isisomerID = False
            istargetID = False
            isligID = True
        else:
            isisomerID = False
            istargetID = False
            isligID = False

    return isligID, istargetID, isisomerID

def get_IDs(directory, isligID, istargetID, isisomerID):
    """Get IDs of ligand target and isomer (if applicable) according to the current architecture."""

    if isisomerID:
        isomerID = directory.split('/')[-1]
        if istargetID:
            targetID = directory.split('/')[-2]
            if isligID:
                ligID = directory.split('/')[-3]
            else:
                ligID = None
        elif isligID:
            targetID = None
            ligID = directory.split('/')[-2]
        else:
            targetID = None
            ligID = None
    elif istargetID:
        isomerID = None
        targetID = directory.split('/')[-1]
        if isligID:
            ligID = directory.split('/')[-2]
        else:
            ligID = None
    elif isligID:
        isomerID = None
        targetID = None
        ligID = directory.split('/')[-1]
    else:
        isomerID = None
        targetID = None
        ligID = None

    return ligID, targetID, isomerID

def check_directories(dirs):
    if dirs != ['.']:
        iscwd = False
        for jdx, dir in enumerate(dirs):
            isligID, istargetID, isisomerID = check_architecture(dir)
            if jdx == 0:
                isligID_ref = isligID
                istargetID_ref = istargetID
                isisomerID_ref = isisomerID
            elif isligID != isligID_ref or istargetID != istargetID_ref or isisomerID != isisomerID_ref:
                raise ValueError("%s*/%s*/%s* architecture architecture inconsistent between folders!"%(ligdir_prefix,tardir_prefix,isodir_prefix))
    else:
        iscwd = True
        isligID = False
        istargetID = False
        isisomerID = False

    return iscwd, isligID, istargetID, isisomerID

