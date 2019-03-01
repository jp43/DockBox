import os
import sys
import shutil
import glide
import subprocess
import socket
import method
import license

required_programs = ['prepwizard', 'ifd', 'ligprep', 'pdbconvert']
default_settings = {'rmsd': '0.18', 'cv_cutoff': '100.00', 'hbond_cutoff' : '-0.05', 'poses_per_lig_dock1': '4', 'poses_per_lig_dock2': '5', 'precision': 'SP'}

class Ifd(method.DockingMethod):

    def __init__(self, instance, site, options):

        super(Ifd, self).__init__(instance, site, options)

        # set box center
        center = site[1] # set box
        self.options['grid_center'] = ','.join(map(str.strip, center.split(',')))

        # set box size
        boxsize = site[2]
        boxsize = map(float, boxsize.split(','))
        self.options['innerbox'] = '%.3f'%max(boxsize)

        self.tmpdirline = ""
        if 'tmpdir' in self.options:
            self.tmpdirline = "export SCHRODINGER_TMPDIR=%s"%self.options['tmpdir']

    def write_docking_script(self, filename, file_r, file_l):

        locals().update(self.options)

        prepwizard_cmd = license.eval("prepwizard -fix %(file_r)s target.mae"%locals(), 'schrodinger')

        hostname = socket.gethostname()
        ifd_cmd = license.eval("ifd -WAIT -SUBHOST %(hostname)s -NPRIMECPU 1 -NGLIDECPU 1 dock.inp"%locals(), 'schrodinger')

        tmpdirline = self.tmpdirline

        # write ifd script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash
%(tmpdirline)s

# (A) Prepare receptor
%(prepwizard_cmd)s

# (B) convert ligand to maestro format
structconvert -imol2 %(file_l)s -omae lig.mae

# perform docking
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
  BINDING_SITE coords %(grid_center)s
  RECEPTOR_SCALE 0.80
  RECEPTOR_CCUT 0.25
  LIGAND_FILE lig.mae
  LIGANDS_TO_DOCK all
  LIGAND_CCUT 0.15
  LIGAND_SCALE 1.00
  MULTI_LIG_CONF no
  CV_CUTOFF %(cv_cutoff)s
  HBOND_CUTOFF %(hbond_cutoff)s
  INNER_BOX %(innerbox)s
  MINIMUM_POSES 1
  MAX_POSESPERLIG %(poses_pre_fit)s
  AMIDE_MODE penal
  OUTER_BOX auto
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
  BINDING_SITE coords %(grid_center)s
  RECEPTOR_SCALE 1.00
  RECEPTOR_CCUT 0.25
  LIGAND_FILE lig.mae
  LIGANDS_TO_DOCK self
  LIGAND_CCUT 0.15
  LIGAND_SCALE 0.80
  MULTI_LIG_CONF yes
  CV_CUTOFF 0.0
  HBOND_CUTOFF 0.0
  INNER_BOX %(innerbox)s
  MAX_POSESPERLIG %(poses_per_lig)s
  AMIDE_MODE penal
  OUTER_BOX auto
  PRECISION %(precision)s

# Scoring
# Compile the IFD Score, consisting of the GlideScore for the Glide Redocking plus 5 percent of the Prime Energy from the Prime Refinement.
STAGE SCORING
  SCORE_NAME r_psp_IFDScore
  TERM 1.0,r_i_glide_gscore,0
  TERM 0.05,r_psp_Prime_Energy,1
  REPORT_FILE report.csv" > dock.inp

%(ifd_cmd)s""" % locals()
            file.write(script)

def extract_docking_results(config):
    sys.exit()
