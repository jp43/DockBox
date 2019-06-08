import os
import sys
import subprocess
import shutil
from glob import glob

from mdkit.utility import mol2
import method

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'prepare_dpf4.py', 'prepare_gpf4.py', 'autogrid4', 'autodock4', 'babel']

default_settings = {'ga_run': '100', 'spacing': '0.3'}

class ADBased(method.DockingMethod):

    def write_rescoring_script(self, filename, file_r, file_l):
        self.write_docking_script(filename, file_r, file_l, rescoring=True)

    def update_output_mol2files(self, sample=None):
        # number of mol2 files generated
        n_files_l = len(glob('lig-*.mol2'))

        mgltools_path = subprocess.check_output('which prepare_ligand4.py', shell=True, executable='/bin/bash')
        mgltools_path = '/'.join(mgltools_path.split('/')[:-3]) 

        for idx in range(n_files_l):
            mol2file = 'lig-%s.mol2'%(idx+1)
            mol2.update_mol2file(mol2file, mol2file, ADupdate=sample, unique=True, mask=['h','H'])
            mol2.arrange_hydrogens(mol2file, 'tmp.mol2', path=mgltools_path)
            shutil.move('tmp.mol2', mol2file)

    def write_check_ligand_pdbqt_script(self, filename):

        with open(filename, 'w') as ff:
            content ="""import os
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
            ff.write(content)

    def write_check_ions_script(self, filename):

        with open(filename, 'w') as file:
            script = """import sys
import shutil
from tempfile import mkstemp

from mdkit.amber.ambertools import load_atomic_ions

# first all residues are supposed to be recognized
are_unrecognized_residues = False

# check if and which atoms were not recognized
unrecognized_residues = []
with open(sys.argv[2], 'r') as logf:
    for line in logf:
        if line.startswith('Sorry, there are no Gasteiger parameters available for atom'):
            are_unrecognized_residues = True
            resname = line.split()[-1].split(':')[0]
            resname = ''.join([i for i in resname if not i.isdigit()])
            unrecognized_residues.append(resname)

if are_unrecognized_residues:

    ions_amber = load_atomic_ions()
    print "No charges specified for ion(s) " + ', '.join(unrecognized_residues)
    print "Attributing formal charges..."

    # update .pdbqt file for the receptor
    fh, abs_path = mkstemp()

    with open(abs_path, 'w') as tempf:
        with open(sys.argv[1], 'r') as ff:

            for line in ff:
                is_ion = False

                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    if resname in unrecognized_residues:
                        assert resname in ions_amber
                        charge = "%.3f"%ions_amber[resname]
                        is_ion = True

                if is_ion:
                    tempf.write(line[:70] + ' '*(6-len(charge)) + charge + line[76:])
                else:
                    tempf.write(line)

    shutil.move(abs_path, sys.argv[1])"""
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
        #TODO: add treatment of ions for autogrid: http://autodock.scripps.edu/faqs-help/how-to/adding-new-atom-parameters-to-autodock

        # create flags with specified options for autogrid and autodock
        autogrid_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in self.autogrid_options.iteritems()])
        autodock_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in self.autodock_options.iteritems()])

        self.write_check_ligand_pdbqt_script('check_ligand_pdbqt.py')
        self.write_check_ions_script('check_ions.py')

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
            with open(filename, 'w') as ff:
                script ="""#!/bin/bash
set -e

MGLPATH=`which prepare_ligand4.py`
MGLPATH=`python -c "print '/'.join('$MGLPATH'.split('/')[:-3])"`
export PYTHONPATH=$PYTHONPATH:$MGLPATH

# prepare ligand
prepare_ligand4.py -l %(file_l)s -o lig.pdbqt
python check_ligand_pdbqt.py lig.pdbqt

# prepare receptor
prepare_receptor4.py -U nphs_lps_waters -r %(file_r)s -o target.pdbqt &> prepare_receptor4.log
python check_ions.py target.pdbqt prepare_receptor4.log

# run autogrid
prepare_gpf4.py -l lig.pdbqt -r target.pdbqt -o grid.gpf %(autogrid_options_flag)s
autogrid4 -p grid.gpf -l grid.glg

# prepare .dpf file
%(ga_num_evals_lines)s
prepare_dpf4.py -l lig.pdbqt -r target.pbdqt -o dock.dpf -p move=lig.pdbqt %(autodock_options_flag)s $ga_num_evals_flag

# run autodock
autodock4 -p dock.dpf -l dock.dlg"""% locals()
                ff.write(script)
 
        else:
            # write autodock script for rescoring
            with open(filename, 'w') as ff:
                script ="""#!/bin/bash
set -e

MGLPATH=`which prepare_ligand4.py`
MGLPATH=`python -c "print '/'.join('$MGLPATH'.split('/')[:-3])"`
export PYTHONPATH=$PYTHONPATH:$MGLPATH

# prepare ligand
prepare_ligand4.py -l %(file_l)s -o lig.pdbqt
python check_ligand_pdbqt.py lig.pdbqt

# prepare receptor only once
if [ ! -f target.pdbqt ]; then
  prepare_receptor4.py -U nphs_lps_waters -r %(file_r)s -o target.pdbqt > prepare_receptor4.log
  python check_ions.py target.pdbqt prepare_receptor4.log
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
                ff.write(script)

    def extract_docking_results(self, file_s, input_file_r, input_file_l):
        """Extract output structures in .mol2 formats"""

        try:
            subprocess.check_output('babel -ad -ipdbqt dock.dlg -omol2 lig-.mol2 -m &>/dev/null', shell=True, executable='/bin/bash')
            self.update_output_mol2files(sample=input_file_l)
            poses_extracted = True
        except:
            mol2files = glob('lig-*.mol2')
            if mol2files: # remove poses if exist
                for mol2file in mol2files:
                    os.remove(mol2file)
            poses_extracted = False

        if os.path.exists('dock.dlg'):
            with open('dock.dlg','r') as dlgf:
                with open(file_s, 'w') as sf:
                    line = '' # initialize line
                    for line in dlgf:
                        if line.startswith('DOCKED: USER    Estimated Free Energy of Binding'):
                            if poses_extracted:
                                score = float(line.split()[8])
                                print >> sf, score
                        if 'CLUSTERING HISTOGRAM' in line:
                            break

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

        to_be_removed = ['run_'+self.program+'.sh', 'prepare_receptor4.log', 'check_lig_pdbqt.py'] + list(glob('*map*')) + list(glob('*pdbqt'))

        for filename in to_be_removed:
            if os.path.isfile(filename):
                os.remove(filename)

