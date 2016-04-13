import os
import sys
import glob
import shutil
import subprocess

import method
import licence.check as chkl

required_programs = ['prepwizard', 'glide', 'ligprep', 'glide_sort', 'pdbconvert']

default_settings = {'poses_per_lig': '10', 'pose_rmsd': '0.5', 'precision': 'SP'}

class Glide(method.DockingMethod):

    def __init__(self, name, site, options):

        super(Glide, self).__init__(name, site, options)

        # set box center
        center = site[1] # set box
        self.options['grid_center'] = ', '.join(map(str.strip, center.split(',')))

        # set box size
        boxsize = site[2]
        boxsize = map(str.strip, boxsize.split(','))
        self.options['innerbox'] = ', '.join(map(str,map(int,map(float,boxsize))))

        outerbox = []
        for idx, xyz in enumerate(['x', 'y', 'z']):
            self.options['act'+xyz+'range'] = str(float(boxsize[idx]) + 10.0)
            outerbox.append(self.options['act'+xyz+'range'])

        self.options['outerbox'] = ', '.join(outerbox)

        self.tmpdirline = ""
        if 'tmpdir' in self.options:
            self.tmpdirline = "export SCHRODINGER_TMPDIR=%(tmpdir)s"%self.options['tmpdir']

    def write_docking_script(self, filename, file_r, file_l):
        """ Write docking script for glide """
        locals().update(self.options)

        # prepare protein cmd (the protein structure is already assumed to be minimized/protonated with prepwizard)
        prepwizard_cmd = chkl.eval("prepwizard -WAIT -noprotassign -nohtreat -noimpref %(file_r)s target.mae"%locals(), 'schrodinger')    

        # prepare grid cmd
        glide_grid_cmd = chkl.eval("glide -WAIT grid.in", 'schrodinger') # grid prepare
        structconvert_cmd = chkl.eval("structconvert -imol2 %(file_l)s -omae lig.mae"%locals(), 'schrodinger')
    
        # docking cmd
        glide_dock_cmd = chkl.eval("glide -WAIT dock.in", 'schrodinger') # docking command

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
%(structconvert_cmd)s

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
 
    def extract_docking_results(self, file_r, file_l, file_s, input_file_r, extract):
    
        #glide_sort_cmd = chkl.eval("glide_sort -r sort.rept dock_pv.maegz -o dock_sorted.mae", 'schrodinger') # cmd to extract results
        subprocess.call(chkl.eval("pdbconvert -brief -imae dock_sorted.mae -opdb dock_sorted.pdb", 'schrodinger', redirect='/dev/null'), shell=True)
        shutil.copyfile(input_file_r, file_r)
    
        if extract == 'lowest':
            # write file with ligand
            with open(file_l, 'w') as ffout:
                with open('dock_sorted-2.pdb', 'r') as ffin:
                    line = ffin.next()
                    while not line.startswith('MODEL'):
                        line = ffin.next()
                    ffout.write(line)
                    for line in ffin:
                        ffout.write(line)
        elif extract == 'all':
            # write file with ligand
            with open(file_l, 'w') as ffout:
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
    
        # extract scores
        with open('dock.rept', 'r') as ffin:
            with open(file_s, 'w') as ffout:
                line = ffin.next()
                while not line.startswith('===='):
                    line = ffin.next()
                while True:
                    line = ffin.next()
                    if line.strip():
                        print >> ffout, line.split()[2]
                        if extract == 'lowest':
                            break
                    else:
                        break

    def get_tmpdir_line(self):
        if self.options['tmpdir']:
            line = "export SCHRODINGER_TMPDIR=%(tmpdir)s"%locals()
        else:
            line = ""

    def write_rescoring_script(self, filename, file_r, file_l):
        """Rescore using Glide SP scoring function"""

        # prepare protein cmd (the protein structure is already assumed to be minimized/protonated with prepwizard)
        prepwizard_cmd = chkl.eval("prepwizard -WAIT -noprotassign -nohtreat -noimpref %(file_r)s target.mae"%locals(), 'schrodinger')

        # prepare grid cmd
        glide_grid_cmd = chkl.eval("glide -WAIT grid.in", 'schrodinger') # grid prepare
        structconvert_cmd = chkl.eval("structconvert -imol2 %(file_l)s -omae lig.mae"%locals(), 'schrodinger')

        # scoring cmd
        glide_dock_cmd = chkl.eval("glide -WAIT dock.in", 'schrodinger') # docking command
        tmpdirline = self.tmpdirline

        with open(filename, 'w') as file:
            script ="""#!/bin/bash
%(tmpdirline)s

if [ ! -f target.mae ]; then
  # (A) Prepare receptor
  %(prepwizard_cmd)s
fi

if [ ! -f grid.zip ]; then
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
fi

# (C) convert ligand to maestro format
%(structconvert_cmd)s

# (D) perform rescoring
echo "WRITEREPT NO
USECOMPMAE YES
DOCKING_METHOD inplace
GRIDFILE $PWD/grid.zip
LIGANDFILE $PWD/lig.mae
PRECISION SP" > dock.in

%(glide_dock_cmd)s"""% locals()
        file.write(script)
    
    def extract_rescoring_results(self, filename):
    
        with open(filename, 'a') as ff:
            with open('dock.scor', 'r') as outf:
                for line in outf:
                    if line.startswith('   1'):
                        print >> ff, line.split()[2]
    
    def cleanup(self):
        pass
        #for ff in glob.glob('dock_sorted*.pdb'):
        #    os.remove(ff)
