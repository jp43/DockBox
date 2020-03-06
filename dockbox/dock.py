import os
import sys
import method

import shutil
import subprocess
from glob import glob

from mdkit.amber import ambertools

from mdkit.utility import reader
from mdkit.utility import mol2
from mdkit.utility import utils

required_programs = ['chimera', 'dms', 'sphgen_cpp', 'sphere_selector', 'showbox', 'grid', 'dock6']

default_settings = {'probe_radius': '1.4', 'minimum_sphere_radius': '1.4', 'maximum_sphere_radius': '4.0', \
'grid_spacing': '0.3', 'extra_margin': '2.0', 'attractive_exponent': '6', 'repulsive_exponent': '12', \
'max_orientations': '10000', 'num_scored_conformers': '5000', 'nposes': '20', 'charge_method': 'gas', 'rmsd': '2.0', 'grid_dir': None}

class Dock(method.DockingMethod):

    def __init__(self, instance, site, options):

        super(Dock, self).__init__(instance, site, options)
        self.options['center'] = '\"' + ' '.join(map(str.strip, site[1].split(','))) + '\"' # set box center
        self.options['site'] = site[0]

        # set box size
        self.options['boxsize'] = map(float, map(str.strip, site[2].split(',')))
        self.options['sphgen_radius'] = str(max(self.options['boxsize'])/2)

        if self.options['site'] is None:
            self.options['dockdir'] = 'dock'
        else:
            self.options['dockdir'] = 'dock.' + self.options['site']

    def write_rescoring_script(self, filename, file_r, files_l):
        """Rescore using DOCK6 grid scoring function"""

        locals().update(self.options)
        self.write_script_ligand_prep()

        # cat mol2 files into a single mol2
        file_all_poses = 'poses.mol2'

        if self.options['charge_method']:
            amber_version = utils.check_amber_version()
            ambertools.run_antechamber(files_l[0], 'pose-1.mol2', at='sybyl', c=self.options['charge_method'], version=amber_version)
        else:
            shutil.copyfile(files_l[0], 'pose-1.mol2')

        for idx, file_l in enumerate(files_l):
            if idx > 0:
                if self.options['charge_method']:
                    # if not first one, do not regenerate the charges, copy charges generated the first time
                    coords_l = mol2.get_coordinates(file_l)
                    struct = mol2.Reader('pose-1.mol2').next()
                    struct = mol2.replace_coordinates(struct, coords_l)
                    mol2.Writer().write('pose-%i.mol2'%(idx+1), struct)
                else:
                    shutil.copyfile(file_l, 'pose-%i.mol2'%(idx+1))
            subprocess.check_output("cat pose-%i.mol2 >> %s"%(idx+1, file_all_poses), shell=True)
            if idx > 0:
                os.remove('pose-%i.mol2'%(idx+1))

        script ="""#!/bin/bash
set -e

# shift ligand coordinates
python prepare_ligand_dock.py pose-1.mol2 pose-1-centered.mol2 %(center)s\n"""%locals()

        if self.options['grid_dir'] is None:
            script += """\n# remove hydrogens from target
echo "delete element.H
write format pdb #0 target_noH.pdb" > removeH.cmd
chimera --nogui %(file_r)s removeH.cmd
rm -rf removeH.cmd

# prepare receptor (add missing h, add partial charges,...)
echo "import chimera
from DockPrep import prep

models = chimera.openModels.list(modelTypes=[chimera.Molecule])
prep(models)

from WriteMol2 import writeMol2
writeMol2(models, 'target.mol2')" > dockprep.py
chimera --nogui %(file_r)s dockprep.py

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

# selecting spheres within a user-defined radius (sphgen_radius)
sphere_selector target_noH_site.sph pose-1-centered.mol2 %(sphgen_radius)s

# create box - the second argument in the file showbox.in
# is the extra margin to also be enclosed to the box (angstroms)
echo "Y
%(extra_margin)s
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

echo "ligand_atom_file %(file_all_poses)s
limit_max_ligands no
skip_molecule no
read_mol_solvation no
calculate_rmsd no
use_database_filter no
orient_ligand no
use_internal_energy yes
internal_energy_rep_exp 12
flexible_ligand no
bump_filter no
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
amber_score_secondary no
minimize_ligand no
atom_model all
vdw_defn_file $vdwfile
flex_defn_file $flexfile
flex_drive_file $flexdfile
ligand_outfile_prefix poses_out
write_orientations no
num_scored_conformers 1
rank_ligands no" > dock6.in

dock6 -i dock6.in > dock.out\n"""%locals()
            file.write(script)
        else:
            # get directory where grid files are located
            grid_prefix = self.options['grid_dir'] + '/' + self.options['dockdir'] + '/grid'

            # check if grid file exists
            if os.path.isfile(grid_prefix+'.in'):
                # copy grid files to avoid opening the same file from multiple locations
                for gridfile in glob(grid_prefix+'*'):
                    basename = os.path.basename(gridfile)
                    shutil.copyfile(gridfile, basename)
            else:
                raise ValueError('No grid file detected in specified location %s'%self.options['grid_dir'])

        script += """\ndock6path=`which dock6`
vdwfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/vdw_AMBER_parm99.defn'"`
flexfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/flex.defn'"`
flexdfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/flex_drive.tbl'"`

echo "ligand_atom_file %(file_all_poses)s
limit_max_ligands no
skip_molecule no
read_mol_solvation no
calculate_rmsd no
use_database_filter no
orient_ligand no
use_internal_energy yes
internal_energy_rep_exp 12
flexible_ligand no
bump_filter no
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
amber_score_secondary no
minimize_ligand no
atom_model all
vdw_defn_file $vdwfile
flex_defn_file $flexfile
flex_drive_file $flexdfile
ligand_outfile_prefix poses_out
write_orientations no
num_scored_conformers 1
rank_ligands no" > dock6.in

dock6 -i dock6.in > dock.out\n"""%locals()

        # write DOCK6 rescoring script
        with open(filename, 'w') as file:
                file.write(script)

    def write_docking_script(self, filename, file_r, file_l):
        """Dock using DOCK6 flexible docking with grid scoring as primary score"""

        locals().update(self.options)
        self.write_script_ligand_prep()

        if self.options['charge_method']:
            amber_version = utils.check_amber_version()
            ambertools.run_antechamber(file_l, 'ligand-ref.mol2', at='sybyl', c=self.options['charge_method'], version=amber_version)
        else:
            shutil.copyfile(file_l, 'ligand-ref.mol2')

        script ="""#!/bin/bash
set -e

# shift ligand coordinates
python prepare_ligand_dock.py ligand-ref.mol2 ligand-ref-centered.mol2 %(center)s\n"""%locals()

        if self.options['grid_dir'] is None:
            script += """\n# remove hydrogens from target
echo "delete element.H
write format pdb #0 target_noH.pdb" > removeH.cmd
chimera --nogui %(file_r)s removeH.cmd
rm -rf removeH.cmd

# prepare receptor (add missing h, add partial charges,...)
echo "import chimera
from DockPrep import prep

models = chimera.openModels.list(modelTypes=[chimera.Molecule])
prep(models)

from WriteMol2 import writeMol2
writeMol2(models, 'target.mol2')" > dockprep.py
chimera --nogui %(file_r)s dockprep.py

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

# selecting spheres within a user-defined radius (sphgen_radius)
sphere_selector target_noH_site.sph ligand-ref-centered.mol2 %(sphgen_radius)s

# create box - the second argument in the file showbox.in
# is the extra margin to also be enclosed to the box (angstroms)
echo "Y
%(extra_margin)s
selected_spheres.sph
1
target_noH_box.pdb" > showbox.in
showbox < showbox.in

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

# create box - the second argument in the file showbox.in
# is the extra margin to also be enclosed to the box (angstroms)
echo "Y
%(extra_margin)s
selected_spheres.sph
1
target_noH_box.pdb" > showbox.in
showbox < showbox.in\n"""%locals()
        else:
            # get directory where grid files are located
            grid_prefix = self.options['grid_dir'] + '/' + self.options['dockdir'] + '/grid'

            # check if grid file exists
            if os.path.isfile(grid_prefix+'.in'):
                # copy grid files to avoid opening the same file from multiple locations
                for gridfile in glob(grid_prefix+'*'):
                    basename = os.path.basename(gridfile)
                    shutil.copyfile(gridfile, basename)
            else:
                raise ValueError('No grid file detected in specified location %s'%self.options['grid_dir'])

            sphfile = self.options['grid_dir'] + '/' + self.options['dockdir'] + '/selected_spheres.sph'
            # check if sphere file exists
            if os.path.isfile(sphfile):
                shutil.copyfile(sphfile, 'selected_spheres.sph')
            else:
                raise ValueError('No selected_spheres.sph file detected in specified location %s'%self.options['grid_dir'])

        script += """\ndock6path=`which dock6`
vdwfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/vdw_AMBER_parm99.defn'"`
flexfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/flex.defn'"`
flexdfile=`python -c "print '/'.join('$dock6path'.split('/')[:-2]) + '/parameters/flex_drive.tbl'"`

# flexible docking using grid score as primary score and no secondary score
echo "ligand_atom_file ligand-ref-centered.mol2
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
ligand_outfile_prefix poses_out
write_orientations no
num_scored_conformers %(num_scored_conformers)s
write_conformations no
cluster_conformations yes
cluster_rmsd_threshold %(rmsd)s
rank_ligands no" > dock6.in

dock6 -i dock6.in\n"""%locals()

        # write DOCK6 script
        with open(filename, 'w') as file:
                file.write(script)

    def extract_docking_results(self, file_s, input_file_r, input_file_l):
    
        # save scores
        if os.path.isfile('poses_out_scored.mol2'):
            with open('poses_out_scored.mol2', 'r') as ffin:
                with open(file_s, 'w') as ffout:
                    idx = 0
                    for line in ffin:
                        if line.startswith('##########    Grid Score:'):
                            ffout.write(line.split()[3]+'\n')
                            idx += 1
                        if idx == int(self.options['nposes']):
                            break 

            # create multiple mol2 files
            ligname = reader.open('poses_out_scored.mol2').ligname
            mol2.update_mol2file('poses_out_scored.mol2', 'pose-.mol2', ligname=ligname, multi=True, last=int(self.options['nposes']))
        else:
            open(file_s, 'w').close()

    def extract_rescoring_results(self, filename, nligands=None):

        with open(filename, 'a') as ff:
            with open('dock.out', 'r') as outf:
                for line in outf:
                    if line.strip().startswith('Grid Score:'):
                        line_s = line.split()
                        if len(line_s) > 2:
                            ff.write(line.split()[2]+'\n')
                        else:
                            ff.write('NaN\n')
                    elif line.strip().startswith('ERROR:  Conformation could not be scored.'):
                        ff.write('NaN\n')

    def write_script_ligand_prep(self):

        with open('prepare_ligand_dock.py', 'w') as file:
            script ="""import os
import sys
import numpy as np
import shutil

from mdkit.utility import utils
from mdkit.utility import mol2

# read mol2 file
mol2file = sys.argv[1]
new_mol2file = sys.argv[2]
center = map(float,(sys.argv[3]).split())

coords = np.array(mol2.get_coordinates(mol2file))
cog = utils.center_of_geometry(coords)
coords = coords - (cog - center)

idx = 0
with open(new_mol2file, 'w') as nmol2f:
    with open(mol2file, 'r') as mol2f:
        is_structure = False
        for line in mol2f:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
                nmol2f.write(line)
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
                nmol2f.write(line)
            elif is_structure:
                new_coords = [format(coord, '.4f') for coord in coords[idx]]
                newline = line[:16] + ' '*(10-len(new_coords[0])) + str(new_coords[0]) + \
' '*(10-len(new_coords[1])) + str(new_coords[1]) + ' '*(10-len(new_coords[2])) + str(new_coords[2]) + line[46:]
                nmol2f.write(newline)
                idx += 1
            else:
                nmol2f.write(line)"""%locals()
            file.write(script)
