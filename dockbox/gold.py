import os
import sys
import shutil
import subprocess

import method
import license

from glob import glob
from mdkit.utility import reader
from mdkit.utility import mol2

required_programs = ['gold_auto']

default_settings = {'nposes': '20'}

class Gold(method.DockingMethod):

    def __init__(self, instance, site, options):

        super(Gold, self).__init__(instance, site, options)

        # set box origin
        self.options['origin'] = ' '.join(map(str.strip, site[1].split(',')))

        # set box radius
        self.options['boxsize'] = map(float, map(str.strip, site[2].split(',')))
        self.options['radius'] = str(max(self.options['boxsize'])/2)

    def write_docking_script(self, filename, file_r, file_l):

        locals().update(self.options)

        dock_cmd = license.wrap_command("gold_auto gold.conf", 'gold') # cmd for docking

        # write autodock script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash

echo "  GOLD CONFIGURATION FILE

  AUTOMATIC SETTINGS
autoscale = 0.5

  POPULATION
popsiz = auto
select_pressure = auto
n_islands = auto
maxops = auto
niche_siz = auto

  GENETIC OPERATORS
pt_crosswt = auto
allele_mutatewt = auto
migratewt = auto

  FLOOD FILL
radius = %(radius)s
origin = %(origin)s
do_cavity = 1
floodfill_atom_no = 0
cavity_file = 
floodfill_center = point

  DATA FILES
ligand_data_file %(file_l)s %(nposes)s
param_file = DEFAULT
set_ligand_atom_types = 1
set_protein_atom_types = 0
directory = .
tordist_file = DEFAULT
make_subdirs = 0
save_lone_pairs = 1
fit_points_file = f it_pts.mol2
read_fitpts = 0

  FLAGS
internal_ligand_h_bonds = 0
flip_free_corners = 0
match_ring_templates = 0
flip_amide_bonds = 0
flip_planar_n = 1 flip_ring_NRR flip_ring_NHR
flip_pyramidal_n = 0
rotate_carboxylic_oh = flip
use_tordist = 1
postprocess_bonds = 1
rotatable_bond_override_file = DEFAULT
solvate_all = 1

  TERMINATION
early_termination = 1
n_top_solutions = 3
rms_tolerance = 1.5

  CONSTRAINTS
force_constraints = 0

  COVALENT BONDING
covalent = 0

  SAVE OPTIONS
save_score_in_file = 1
save_protein_torsions = 1

  FITNESS FUNCTION SETTINGS
initial_virtual_pt_match_max = 3
relative_ligand_energy = 1
gold_fitfunc_path = goldscore
start_vdw_linear_cutoff = 6
score_param_file = DEFAULT

  PROTEIN DATA
protein_datafile = %(file_r)s" > gold.conf

%(dock_cmd)s"""% locals()
            file.write(script)
    
    def extract_docking_results(self, file_s, input_file_r, input_file_l):

        for idx, mol2file in enumerate(sorted(glob('gold_soln_*_m1_*.mol2'))):
            output_mol2file = 'lig-%s.mol2'%(idx+1)
            mol2.update_mol2file(mol2file, output_mol2file, remove='LP')

