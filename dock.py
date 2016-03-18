import os
import sys
import subprocess
import shutil
import glob

import util.pdbtools as pdbt

required_programs = ['chimera', 'dms', 'sphgen_cpp', 'sphere_selector', 'showbox', 'grid', 'dock6']

default_settings = {'probe_radius': '1.4', 'minimum_sphere_radius': '1.4', 'maximum_sphere_radius': '4.0', \
'sphere_radius': '10.0', 'box_size': '10.0', 'grid_spacing': '0.3', 'attractive_exponent': '6', 'repulsive_exponent': '12', \
'max_orientations': '10000', 'num_scored_conformers': '10'}

known_settings = {'herg': {'center': '"3.966 8.683 11.093"'}, \
'herg-cut': {'center': '"3.966 8.683 11.093"'}, \
'herg-inactivated': {'center': '"0.000 0.000 -5.000"'}}

required_settings_names = ['center']

def write_docking_script(filename, input_file_r, input_file_l, config):

    write_shift_coordinates_script(config)
  
    # convert sdffile to PDBfile
    subprocess.check_call('babel %s %s 2>/dev/null'%(input_file_l,'lig.pdb'), shell=True)

    # write autodock script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash
set -e

# remove hydrogens from target
echo "delete element.H
write format pdb #0 target_noH.pdb" > removeH.cmd
chimera --nogui --silent %(input_file_r)s removeH.cmd
rm -rf removeH.cmd

# convert pdb to mol2
echo "write format mol2 #0 target.mol2" > pdb2mol.cmd
chimera --nogui --silent %(input_file_r)s pdb2mol.cmd
rm -rf pdb2mol.cmd

# generating receptor surface
dms target_noH.pdb -n -w %(probe_radius)s -v -o target_noH.ms

# generating spheres
echo "target_noH.ms
R
X
0.0
%(maximum_sphere_radius)s
%(minimum_sphere_radius)s
target_noH_site.sph" > INSPH
sphgen_cpp

# shift ligand coordintates
python shift_lig_coords.py lig.pdb lig_s.pdb %(center)s

# convert pdb to mol2
echo "write format mol2 #0 lig.mol2" > pdb2mol.cmd
chimera --nogui --silent lig_s.pdb pdb2mol.cmd
rm -rf pdb2mol.cmd
# selecting spheres
sphere_selector target_noH_site.sph lig.mol2 %(sphere_radius)s

# create box - the second argument in the file showbox.in
# is the radius used with sphere_selector
echo "Y
%(box_size)s
selected_spheres.sph
1
target_noH_box.pdb" > showbox.in
showbox < showbox.in

dock6path=`which dock6`
vdwfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/vdw_AMBER_parm99.defn'"`
flexfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/flex.defn'"`
flexdfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/flex_drive.tbl'"`

# create grid
echo "compute_grids yes
grid_spacing %(grid_spacing)s
output_molecule no
contact_score yes
energy_score yes
energy_cutoff_distance 9999
atom_model a
attractive_exponent %(attractive_exponent)s
repulsive_exponent %(repulsive_exponent)s
distance_dielectric yes
dielectric_factor 4
bump_filter yes
bump_overlap 0.75
receptor_file target.mol2
box_file target_noH_box.pdb
vdw_definition_file $vdwfile
score_grid_prefix grid
contact_cutoff_distance 4.5" > grid.in
grid -i grid.in

# flexible docking
echo "ligand_atom_file lig.mol2
limit_max_ligands no
skip_molecule no
read_mol_solvation no
calculate_rmsd no
use_database_filter no
orient_ligand yes
automated_matching yes
receptor_site_file selected_spheres.sph
max_orientations %(max_orientations)s
critical_points no
chemical_matching no
use_ligand_spheres no
use_internal_energy yes
internal_energy_rep_exp 12
flexible_ligand yes
user_specified_anchor no
limit_max_anchors no
min_anchor_size 5
pruning_use_clustering yes
pruning_max_orients 1000
pruning_clustering_cutoff 100
pruning_conformer_score_cutoff 100
use_clash_overlap yes
clash_overlap 0.5
write_growth_tree no
bump_filter yes
bump_grid_prefix grid
max_bumps_anchor 12
max_bumps_growth 12
score_molecules yes
contact_score_primary no
contact_score_secondary no
grid_score_primary yes
grid_score_secondary no
grid_score_rep_rad_scale 1
grid_score_vdw_scale 1
grid_score_es_scale 1
grid_score_grid_prefix grid
multigrid_score_secondary no
dock3.5_score_secondary no
continuous_score_secondary no
descriptor_score_secondary no
gbsa_zou_score_secondary no
gbsa_hawkins_score_secondary no
SASA_descriptor_score_secondary no
pbsa_score_secondary no
amber_score_secondary no
minimize_ligand yes
minimize_anchor yes
minimize_flexible_growth yes
use_advanced_simplex_parameters no
simplex_max_cycles 1
simplex_score_converge 0.1
simplex_cycle_converge 1.0
simplex_trans_step 1.0
simplex_rot_step 0.1
simplex_tors_step 10.0
simplex_anchor_max_iterations 1000
simplex_grow_max_iterations 1000
simplex_grow_tors_premin_iterations 0
simplex_random_seed 0
simplex_restraint_min no
atom_model all
vdw_defn_file $vdwfile
flex_defn_file $flexfile
flex_drive_file $flexdfile
ligand_outfile_prefix lig_out
write_orientations no
num_scored_conformers %(num_scored_conformers)s
write_conformations yes
cluster_conformations yes
cluster_rmsd_threshold 2.0
rank_ligands no" > dock6.in

dock6 -i dock6.in"""%dict(dict(locals()).items()+config.options['dock'].items())
        file.write(script)

def extract_docking_results(file_r, file_l, file_s, config):

    shutil.copyfile(config.input_file_r, file_r)

    # convert mol2 to pdb
    if config.extract == 'lowest':
        with open('pdb2mol.cmd', 'w') as ff:
            script = "write format pdb #0.1 %s" %file_l
            ff.write(script)
    elif config.extract == 'all':
        with open('pdb2mol.cmd', 'w') as ff:
            script = "write format pdb #0 %s" %file_l
            ff.write(script)

    subprocess.call('chimera --nogui --silent lig_out_scored.mol2 pdb2mol.cmd', shell=True)
    os.remove('pdb2mol.cmd')

    # save scores
    with open('lig_out_scored.mol2', 'r') as ffin:
        with open(file_s, 'w') as ffout:
            for line in ffin:
                if line.startswith('##########    Grid Score:'):
                    print >> ffout, line.split()[3]
                    if config.extract == 'lowest':
                        break

def cleanup(config):

    # remove map files
    for ff in glob.glob('grid*'):
        os.remove(ff)

    os.remove('selected_spheres.sph') 
    os.remove('target_noH.ms')

def write_shift_coordinates_script(config):

    with open('shift_lig_coords.py', 'w') as file:
        script ="""import sys
import numpy as np

# read PDB file

pdb = sys.argv[1]
newpdb = sys.argv[2]
center = map(float,(sys.argv[3]).split())

coords = []
with open(pdb, 'r') as pdbf:
    for line in pdbf:
        if line.startswith(('ATOM  ','HETATM')):
            coords.append(map(float,[line[30:38],line[38:46],line[46:54]]))

coords = np.array(coords)

def center_of_geometry(coords):
    cog = np.array([0.0, 0.0, 0.0])
    num_atom = coords.shape[0]
    for idx in xrange(num_atom):
        cog = cog + coords[idx]
    cog = cog / num_atom
    return cog

cog = center_of_geometry(coords)
coords = coords - (cog - center)

idx = 0
with open(newpdb, 'w') as pdbn:
    with open(pdb, 'r') as pdbf:
        for line in pdbf:
            if line.startswith(('ATOM  ','HETATM')):
                #coords_str = map(str, coords[idx]
                newline = line[:30]
                for jdx in range(3):
                    coord = ("%.3f"%coords[idx,jdx]).strip()
                    newline += ' '*(8-len(coord)) + coord
                if len(line) >= 54:
                    newline += line[54:-1]
                idx += 1
                print >> pdbn, newline
            else:
                pdbn.write(line)"""
        file.write(script)
