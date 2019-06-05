import os
import sys
import glob
import shutil
import subprocess
import numpy as np

import method
import license

from mdkit.utility import reader
from mdkit.utility import mol2

required_programs = ['prepwizard', 'glide', 'glide_sort', 'pdbconvert']

default_settings = {'poses_per_lig': '10', 'pose_rmsd': '0.5', 'precision': 'SP', 'use_prepwizard': 'True'}

known_settings = {'precision': ['SP', 'XP'], 'use_prepwizard': ['true', 'false', 'yes', 'no']}

class Glide(method.DockingMethod):

    def __init__(self, instance, site, options):

        super(Glide, self).__init__(instance, site, options)

        # set box center
        center = site[1] # set box
        self.options['grid_center'] = ', '.join(map(str.strip, center.split(',')))

        # set box size
        boxsize = site[2]
        boxsize = map(str.strip, boxsize.split(','))
        self.options['innerbox'] = ', '.join(["%i"%int(float(boxsize[idx])) for idx in range(3)])

        outerbox = []
        for idx, xyz in enumerate(['x', 'y', 'z']):
            self.options['act'+xyz+'range'] = str("%.1f"%float(boxsize[idx]))
            outerbox.append(self.options['act'+xyz+'range'])

        self.options['outerbox'] = ', '.join(outerbox)

        self.tmpdirline = ""
        if 'tmpdir' in self.options:
            self.tmpdirline = "export SCHRODINGER_TMPDIR=%s"%self.options['tmpdir']

        if self.options['use_prepwizard'].lower() in ['yes', 'true']:
            self.use_prepwizard = True
        elif self.options['use_prepwizard'].lower() in ['no', 'false']:
            self.use_prepwizard = False
        else:
            raise ValueError("Value for use_prepwizard non recognized")
            
    def write_docking_script(self, filename, file_r, file_l):
        """ Write docking script for glide """
        locals().update(self.options)

        if self.use_prepwizard:
            # prepare protein cmd (the protein structure is already assumed to be minimized/protonated with prepwizard)
            prepwizard_cmd = license.wrap_command("prepwizard -fix %(file_r)s target.mae"%locals(), 'schrodinger')    
        else:
            prepwizard_cmd = "structconvert -ipdb %(file_r)s -omae target.mae"%locals()

        # prepare grid and docking cmd
        glide_grid_cmd = license.wrap_command("glide grid.in", 'schrodinger')
        glide_dock_cmd = license.wrap_command("glide dock.in", 'schrodinger')

        tmpdirline = self.tmpdirline
    
        # write glide script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash
%(tmpdirline)s

# (A) Prepare receptor
%(prepwizard_cmd)s

# (B) Prepare grid
echo "USECOMPMAE YES
INNERBOX %(innerbox)s
ACTXRANGE %(actxrange)s
ACTYRANGE %(actyrange)s
ACTZRANGE %(actzrange)s
GRID_CENTER %(grid_center)s
OUTERBOX %(outerbox)s
ENTRYTITLE target
GRIDFILE grid.zip
RECEP_FILE target.mae" > grid.in
%(glide_grid_cmd)s

# (C) convert ligand to maestro format
structconvert -imol2 %(file_l)s -omae lig.mae

# (D) perform docking
echo "WRITEREPT YES
USECOMPMAE YES
DOCKING_METHOD confgen
POSES_PER_LIG %(poses_per_lig)s
POSE_RMSD %(pose_rmsd)s
GRIDFILE $PWD/grid.zip
LIGANDFILE $PWD/lig.mae
PRECISION %(precision)s" > dock.in
%(glide_dock_cmd)s"""% locals()
            file.write(script)
 
    def extract_docking_results(self, file_s, input_file_r, input_file_l):
        """Extract Glide docking results""" 

        if os.path.exists('dock_pv.maegz'):
            # (1) cmd to extract results
            subprocess.check_output('glide_sort -r sort.rept dock_pv.maegz -o dock_sorted.mae', shell=True, executable='/bin/bash')

            # (2) convert to .mol2
            subprocess.check_output('mol2convert -n 2: -imae dock_sorted.mae -omol2 dock_sorted.mol2', shell=True, executable='/bin/bash')

            if os.path.exists('dock_sorted.mol2'):
                ligname = reader.open(input_file_l).ligname
                mol2.update_mol2file('dock_sorted.mol2', 'lig-.mol2', ligname=ligname, multi=True)
                # extract scores
                with open('dock.rept', 'r') as ffin:
                    with open(file_s, 'w') as ffout:
                        line = ffin.next()
                        while not line.startswith('===='):
                            line = ffin.next()
                        while True:
                            line = ffin.next()
                            if line.strip():
                                print >> ffout, line[43:51].strip()
                            else:
                                break

    def get_tmpdir_line(self):
        if self.options['tmpdir']:
            line = "export SCHRODINGER_TMPDIR=%(tmpdir)s"%locals()
        else:
            line = ""

    def write_rescoring_script(self, filename, file_r, files_l):
        """Rescore using Glide SP scoring function"""
        locals().update(self.options)

        files_l_joined = ' '.join(files_l)

        if self.use_prepwizard:
            # prepare protein cmd (the protein structure is already assumed to be minimized/protonated with prepwizard)
            prepwizard_cmd = license.wrap_command("prepwizard -fix %(file_r)s target.mae"%locals(), 'schrodinger')
        else:
            prepwizard_cmd = "structconvert -ipdb %(file_r)s -omae target.mae"%locals()


        # prepare grid and scoring cmd
        glide_grid_cmd = license.wrap_command("glide grid.in", 'schrodinger') # grid prepare
        glide_dock_cmd = license.wrap_command("glide dock.in", 'schrodinger') # docking command
        tmpdirline = self.tmpdirline

        with open(filename, 'w') as file:
            script ="""#!/bin/bash
%(tmpdirline)s
cat %(files_l_joined)s > lig.mol2

# (A) Prepare receptor
%(prepwizard_cmd)s

# (B) Prepare grid
echo "USECOMPMAE YES
INNERBOX %(innerbox)s
ACTXRANGE %(actxrange)s
ACTYRANGE %(actyrange)s
ACTZRANGE %(actzrange)s
GRID_CENTER %(grid_center)s
OUTERBOX %(outerbox)s
ENTRYTITLE target
GRIDFILE grid.zip
RECEP_FILE target.mae" > grid.in
%(glide_grid_cmd)s


# (C) convert ligand to maestro format
structconvert -imol2 lig.mol2 -omae lig.mae

# (D) perform rescoring
echo "WRITEREPT YES
USECOMPMAE YES
DOCKING_METHOD inplace
GRIDFILE $PWD/grid.zip
LIGANDFILE $PWD/lig.mae
PRECISION SP" > dock.in

%(glide_dock_cmd)s"""% locals()
            file.write(script)

    def extract_rescoring_results(self, filename, nligands=None):
        idxs = []
        scores = []

        if os.path.exists('dock.scor'):
            with open('dock.scor', 'r') as ffin:
                line = ffin.next()
                while not line.startswith('===='):
                    line = ffin.next()
                while True:
                    line = ffin.next()
                    if line.strip():
                        idxs.append(int(line[36:42].strip()))
                        scores.append(line[43:51].strip())
                    else:
                        break

            scores = np.array(scores)
            scores = scores[np.argsort(idxs)]
        else:
            scores = [ 'NaN' for idx in range(nligands)]

        with open(filename, 'w') as ffout:
            for sc in scores:
                print >> ffout, sc

