import os
import subprocess
import glob
import shutil
import method
from tools import mol2

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'prepare_dpf4.py', 'prepare_gpf4.py', 'autogrid4', 'autodock4', 'babel']

default_settings = {'ga_run': '100', 'spacing': '0.238'}

class ADBased(method.DockingMethod):

    def write_rescoring_script(self, filename, file_r, file_l):
        self.write_docking_script(filename, file_r, file_l, rescoring=True)

    def update_output_mol2files(self, sample=None):
        # number of mol2 files generated
        n_files_l = len(glob.glob('lig-*.mol2'))

        for idx in range(n_files_l):
            mol2file = 'lig-%s.mol2'%(idx+1)
            mol2.update_mol2file(mol2file, mol2file, ADupdate=sample, unique=True, mask=['h','H'])
            mol2.arrange_hydrogens(mol2file, 'tmp.mol2')
            shutil.move('tmp.mol2', mol2file)

    def write_fix_pdbqt_script(self):

        with open('fix_pdbqt.py', 'w') as file:
            script ="""import os
import sys
import shutil

input_file = sys.argv[1]

filename, ext = os.path.splitext(input_file)
file_tmp = filename + '_tmp.pdbqt'

lines_to_be_removed = []
has_branch_started = False
with open(input_file, 'r') as ff:
    for line in ff:
        if has_branch_started:
            has_branch_started = False
            branch_num = start_branch_line.split()[-1]
            if line.split()[1] != branch_num:
                lines_to_be_removed.append(start_branch_line)
                lines_to_be_removed.append('END' + start_branch_line)
        if line.startswith('BRANCH'):
            start_branch_line = line
            has_branch_started = True

if lines_to_be_removed:
    with open(input_file, 'r') as ff:
        with open(file_tmp, 'w') as of:
            for line in ff:
                if line.startswith(('BRANCH', 'ENDBRANCH')) and line in lines_to_be_removed:
                    pass
                else:
                    of.write(line)
    shutil.move(file_tmp, input_file)"""
            file.write(script)

class Autodock(ADBased):

    def __init__(self, instance, site, options):

        super(Autodock, self).__init__(instance, site, options)

        # set box center
        self.options['gridcenter'] = '\"' + ' '.join(map(str.strip, site[1].split(','))) + '\"'
 
        # set box size
        boxsize = map(float, map(str.strip, site[2].split(',')))
        spacing = float(options['spacing'])
        npts = []
        for size in boxsize:
            sz = int(size*1.0/spacing) + 1
            npts.append(str(sz)) # round to the integer above
        self.options['npts'] =  ','.join(npts)

        autogrid_options_names = ['spacing', 'npts', 'gridcenter']
        autodock_options_names = ['ga_run', 'ga_pop_size', 'ga_num_evals', 'ga_num_generations', 'outlev', 'seed']

        self.autogrid_options = {}
        for name in autogrid_options_names:
            if name in options:
                self.autogrid_options[name] = options[name]

        self.autodock_options = {}
        for name in autodock_options_names:
            if name in options:
                self.autodock_options[name] = options[name]

    def write_docking_script(self, filename, file_r, file_l, rescoring=False):

        # create flags with specified options for autogrid and autodock
        autogrid_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in self.autogrid_options.iteritems()])
        autodock_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in self.autodock_options.iteritems()])

        self.write_fix_pdbqt_script()

        if not rescoring:
            if 'ga_num_evals' not in self.options:
                ga_num_evals_lines="""prepare_dpf4.py -l lig.pdbqt -r target.pdbqt -o dock.dpf -p move=lig.pdbqt
ga_num_evals_flag=`python -c \"with open('dock.dpf') as ff:
    for line in ff:
        if line.startswith('torsdof'):
            torsion = int(line.split()[1])
            break
ga_num_evals = min(25000000, 987500 * torsion + 125000)
print \'-p ga_num_evals=%i\'%ga_num_evals\"`"""
            else:
                ga_num_evals_lines=""
 
            # write autodock script
            with open(filename, 'w') as file:
                script ="""#!/bin/bash
set -e
# generate .pdbqt files
prepare_ligand4.py -l %(file_l)s -o lig.pdbqt
python fix_pdbqt.py lig.pdbqt
prepare_receptor4.py -r %(file_r)s -o target.pdbqt

# run autogrid
prepare_gpf4.py -l lig.pdbqt -r target.pdbqt -o grid.gpf %(autogrid_options_flag)s
autogrid4 -p grid.gpf -l grid.glg

# prepare .dpf file
%(ga_num_evals_lines)s
prepare_dpf4.py -l lig.pdbqt -r target.pbdqt -o dock.dpf -p move=lig.pdbqt %(autodock_options_flag)s $ga_num_evals_flag

# run autodock
autodock4 -p dock.dpf -l dock.dlg"""% locals()
                file.write(script)
 
        else:
            # write autodock script for rescoring
            with open(filename, 'w') as file:
                script ="""#!/bin/bash
set -e
# generate .pdbqt files
prepare_ligand4.py -l %(file_l)s -o lig.pdbqt
python fix_pdbqt.py lig.pdbqt
if [ ! -f target.pdbqt ]; then
  prepare_receptor4.py -r %(file_r)s -o target.pdbqt
fi

# run autogrid
if [ ! -f grid.glg ]; then
  prepare_gpf4.py -l lig.pdbqt -r target.pdbqt -o grid.gpf %(autogrid_options_flag)s
  autogrid4 -p grid.gpf -l grid.glg
fi

# prepare .dpf file
if [ ! -f dock.dpf ]; then
  prepare_dpf4.py -l lig.pdbqt -r target.pbdqt -o dock.dpf -p move=lig.pdbqt %(autodock_options_flag)s $ga_num_evals_flag
  # construct new dock.dpf with rescoring options only
  sed -e "1,/about/w tmp.dpf" dock.dpf > /dev/null
  mv tmp.dpf dock.dpf
  echo 'epdb                                 # small molecule to be evaluated' >> dock.dpf
fi

# run autodock
autodock4 -p dock.dpf -l dock.dlg"""% locals()
                file.write(script)

    def extract_docking_results(self, file_s, input_file_r, input_file_l):
        """extract poses and scores from .dlg file"""

        if os.path.exists('dock.dlg'):
            with open('dock.dlg','r') as dlgf:
                with open(file_s, 'w') as sf: 
                    line = '' # initialize line
                    for line in dlgf:
                        if 'Estimated Free Energy of Binding' in line:
                            score = float(line.split()[8])
                            print >> sf, score
                        if 'CLUSTERING HISTOGRAM' in line:
                            break

            try:
                subprocess.check_output('babel -ad -ipdbqt dock.dlg -omol2 lig-.mol2 -m &>/dev/null', shell=True, executable='/bin/bash')
                self.update_output_mol2files(sample=input_file_l)
            except subprocess.CalledProcessError:
                pass

    def extract_rescoring_results(self, filename):
        """extract scores from .dlg file"""
        with open(filename, 'a') as ff:
            if os.path.exists('dock.dlg'):
                with open('dock.dlg', 'r') as outf:
                    has_fe_line = False
                    for line in outf:
                        if line.startswith('epdb: USER    Estimated Free Energy of Binding'):
                            print >> ff, line.split()[8]
                            has_fe_line = True
                    if not has_fe_line:
                        print >> ff, 'NaN'
            else:
                print >> ff, 'NaN'

    def cleanup(self):
        # remove map files
        for ff in glob.glob('*map*'):
            os.remove(ff)

        for ff in glob.glob('*pdbqt'):
            os.remove(ff)
