import os
import sys
import glob
import shutil
import socket
import glide
import subprocess

required_programs = ['prepwizard', 'ifd', 'ligprep', 'pdbconvert']

default_settings = {'rmsd': '0.18', 'cv_cutoff': '100.00', 'hbond_cutoff' : '-0.05', 'inner_box': '10.0', 'outer_box': '30.0'}

known_settings = {'herg': {'binding_site': '3.966,8.683,11.093'}, \
'herg-cut': {'binding_site': '3.966,8.683,11.093'}, \
'herg-inactivated': {'binding_site': '0.000,0.000,-5.000'}}

required_settings_names = ['binding_site']

def write_docking_script(filename, config):

    input_file_l = config.input_file_l
    input_file_r = config.input_file_r

    hostname = socket.gethostname()

    # convert .pdb to sdf
    subprocess.call('babel -ipdb ' + input_file_l  + ' -osdf lig.sdf 2> /dev/null', shell=True)

    # write vina script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash

# (A) Prepare receptor
prepwizard -WAIT -fix %(input_file_r)s target.mae

# (B) prepare ligand
ligprep -WAIT -isd lig.sdf -omae ligprep.mae -nz -ns -nr -nt -i 0

# (C) perform docking
echo "# Global Variables
# These variables affect the entire job, and must all appear
# before the first STAGE declaration. Multiple INPUT_FILE
# entries are supported, as are files containing multiple
# receptor structures. 
# If beginning with an existing Pose Viewer file, simply specify
# it as the INPUT_FILE (making sure the name ends in "_pv.mae")
# and ensure that the first GLIDE_DOCKING stage is commented out.
# The ligand used in producing the Pose Viewer file must also be
# provided to the second GLIDE_DOCKING stage, using the LIGAND_FILE
# keyword.
INPUT_FILE target.mae

# Protein Preparation
# Run a simple constrained minimization of the receptor # structure(s).
STAGE PPREP
  RMSD %(rmsd)s

# In order to temporarily remove the side chains of residues
# (i.e., mutate to Ala) that are blocking the binding site, uncomment
# the following STAGE line, and then specify the sidechains to be
# removed using either one of the two Methods described below.
#STAGE TRIM_SIDECHAINS
#  RESIDUES A:10B

# Glide Docking
# Perform the initial Glide docking, producing a ligand-receptor
# complex for each pose requested/found. If multiple receptor 
# structures are used, the requested number of poses will be 
# generated for each structure.
STAGE GLIDE_DOCKING
  BINDING_SITE coords %(binding_site)s
  RECEPTOR_SCALE 0.80
  RECEPTOR_CCUT 0.25
  LIGAND_FILE ligprep.mae
  LIGANDS_TO_DOCK all
  LIGAND_CCUT 0.15
  LIGAND_SCALE 1.00
  MULTI_LIG_CONF no
  CV_CUTOFF %(cv_cutoff)s
  HBOND_CUTOFF %(hbond_cutoff)s
  INNER_BOX %(inner_box)s
  MINIMUM_POSES 1
  MAX_POSESPERLIG 10
  AMIDE_MODE penal
  OUTER_BOX %(outer_box)s
  PRECISION %(precision)s

# Determine Residue to Refine
# Compile a list of all residues within the specified distance of any pose of the ligand.
STAGE COMPILE_RESIDUE_LIST
  DISTANCE_CUTOFF 5.0

# Prime Refinement
# Optimize the side chains of the residue list compiled previously, then minimize them along with the ligand.
STAGE PRIME_REFINEMENT
  NUMBER_OF_PASSES 1
  USE_MEMBRANE no
  FIX_ATOM_NAMES yes

# Sort and Filter
# Only retain poses with Prime Energies within the # specified range from the lowest energy pose.
STAGE SORT_AND_FILTER
  POSE_FILTER r_psp_Prime_Energy
  POSE_KEEP 30.0

# Sort and Filter
# Only retain the top number of poses specified.
STAGE SORT_AND_FILTER
  POSE_FILTER r_psp_Prime_Energy
  POSE_KEEP 20#

# Glide Docking
# Redock the ligand back into the newly optimized receptor, using default Glide settings.
STAGE GLIDE_DOCKING
  BINDING_SITE coords %(binding_site)s
  RECEPTOR_SCALE 1.00
  RECEPTOR_CCUT 0.25
  LIGAND_FILE ligprep.mae
  LIGANDS_TO_DOCK self
  LIGAND_CCUT 0.15
  LIGAND_SCALE 0.80
  MULTI_LIG_CONF yes
  CV_CUTOFF 0.0
  HBOND_CUTOFF 0.0
  INNER_BOX %(inner_box)s
  MAX_POSESPERLIG 1
  AMIDE_MODE penal
  OUTER_BOX %(outer_box)s
  PRECISION %(precision)s

# Scoring
# Compile the IFD Score, consisting of the GlideScore for the Glide Redocking plus 5 percent of the Prime Energy from the Prime Refinement.
STAGE SCORING
  SCORE_NAME r_psp_IFDScore
  TERM 1.0,r_i_glide_gscore,0
  TERM 0.05,r_psp_Prime_Energy,1
  REPORT_FILE report.csv" > dock.inp

ifd -WAIT -SUBHOST %(hostname)s -NPRIMECPU 1 -NGLIDECPU 1 dock.inp""" % dict(dict(locals()).items() + config.options['ifd'].items())
        file.write(script)

def extract_docking_results(config):

    # TODO: finish it! 
    return

    if config.extract == 'lowest':
        nflag = '-n 1 '
    elif config.extract == 'all':
        nflag = '-n :4 '
    
    subprocess.call('pdbconvert -imae dock-out.maegz -opdb dock-out.pdb -reorder_by_sequence %s'%nflag, shell=True)

    with open('rec-out.pdb', 'w') as fft:
        with open('lig-out.pdb', 'w') as ffl:
            with open('dock-out.pdb', 'r') as ffd:
                while line:
                    if line.startswith('MODEL'):
                        print >> ffl, line
                        print >> fft, line

    # write file with target
    with open('dock_sorted-1.pdb', 'r') as ffin:
        with open('target_out.pdb', 'w') as ffout:
            line = ffin.next()
            while not line.startswith('MODEL'):
                line = ffin.next()
            ffout.write(line)
            for line in ffin:
                ffout.write(line)

    # convert mol2 to pdb
    if config.extract == 'lowest':
        # write file with ligand
        with open('lig_out.pdb', 'w') as ffout:
            with open('dock_sorted-2.pdb', 'r') as ffin:
                line = ffin.next()
                while not line.startswith('MODEL'):
                    line = ffin.next()
                ffout.write(line)
                for line in ffin:
                    ffout.write(line)
    elif config.extract == 'all':
        # write file with ligand
        with open('lig_out.pdb', 'w') as ffout:
            idxs = []
            for pdbfile in glob.glob('dock_sorted-*.pdb'):
                idxs.append(int((pdbfile.split('-')[1]).split('.')[0]))
            maxidx = max(idxs)
            for idx in range(1,maxidx):
                with open('dock_sorted-%s.pdb'%(idx+1), 'r') as ffin:
                    line = ffin.next()
                    while not line.startswith('MODEL'):
                        line = ffin.next()
                    ffout.write('MODEL     %s\n'%idx)
                    for line in ffin:
                        if line.startswith('ENDMDL') or not line.startswith('END'):
                            ffout.write(line)
            ffout.write('END')

    # extract scores
    with open('dock.rept', 'r') as ffin:
        with open('scored.out', 'w') as ffout:
            line = ffin.next()
            while not line.startswith('===='):
                line = ffin.next()
            while True:
                line = ffin.next()
                if line.strip():
                    print >> ffout, line.split()[2]
                    if config.extract == 'lowest':
                        break
                else:
                    break

def cleanup(config):

    return
    for ff in glob.glob('dock_sorted*.pdb'):
        os.remove(ff)

