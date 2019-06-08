import os
import sys
from glob import glob
import shutil
import subprocess
import autodock

from mdkit.utility import mol2

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'vina', 'babel']

default_settings = {'cpu': '1', 'num_modes': '9', 'energy_range': '3'}

class Vina(autodock.ADBased):

    def __init__(self, instance, site, options):

        super(Vina, self).__init__(instance, site, options)

        center = map(str.strip, site[1].split(','))
        boxsize = map(str.strip, site[2].split(','))

        for idx, xyz in enumerate(['x', 'y', 'z']):
            self.options['center_'+xyz] = center[idx]
            self.options['size_'+xyz] = boxsize[idx]

    def write_docking_script(self, filename, file_r, file_l, rescoring=False):
        """write docking script for Vina"""

        locals().update(self.options)

        self.write_check_ligand_pdbqt_script('check_ligand_pdbqt.py')
        self.write_check_ions_script('check_ions.py')

        # write vina config file
        with open('vina.config', 'w') as cf:
            # write mandatory options
            print >> cf, 'receptor = target.pdbqt'
            print >> cf, 'ligand = lig.pdbqt'
            # write other options
            for key, value in self.options.iteritems():
                print >> cf, key + ' = ' + value

        # write vina script
        if not rescoring:
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

# run vina
vina --config vina.config 1> vina.out 2> vina.err"""% locals()
                ff.write(script)
        else:
            with open(filename, 'w') as ff:
                script ="""#!/bin/bash
set -e

MGLPATH=`which prepare_ligand4.py`
MGLPATH=`python -c "print '/'.join('$MGLPATH'.split('/')[:-3])"`
export PYTHONPATH=$PYTHONPATH:$MGLPATH

# prepare ligand
prepare_ligand4.py -l %(file_l)s -o lig.pdbqt
python check_ligand_pdbqt.py lig.pdbqt

if [ ! -f target.pdbqt ]; then
  prepare_receptor4.py -U nphs_lps_waters -r %(file_r)s -o target.pdbqt > prepare_receptor4.log
  python check_ions.py target.pdbqt prepare_receptor4.log
fi

# run vina
vina --score_only --config vina.config > vina.out"""% locals()
                ff.write(script)

    def extract_docking_results(self, file_s, input_file_r, input_file_l):
        """Extract output structures in .mol2 formats"""

        try:
            subprocess.check_output('babel -ipdbqt lig_out.pdbqt -omol2 lig-.mol2 -m &>/dev/null',shell=True, executable='/bin/bash')
            self.update_output_mol2files(sample=input_file_l)
            poses_extracted = True
        except:
            mol2files = glob('lig-*.mol2')
            if mol2files: # remove poses if exist
                for mol2file in mol2files:
                    os.remove(mol2file)
            poses_extracted = False

        if os.path.exists('lig_out.pdbqt'):
            # extract structures from .pdbqt file 
            with open('lig_out.pdbqt','r') as pdbqtf:
                with open(file_s, 'w') as sf:
                    for line in pdbqtf:
                        if line.startswith('REMARK VINA RESULT:'):
                            if poses_extracted:
                                score = float(line[19:].split()[0])
                                print >> sf, score

    def write_rescoring_script(self, filename, file_r, file_l):
        self.write_docking_script(filename, file_r, file_l, rescoring=True)
    
    def extract_rescoring_results(self, filename):

        with open(filename, 'a') as ff:
            with open('vina.out', 'r') as outf:
                for line in outf:
                    if line.startswith('Affinity:'):
                        print >> ff, line.split()[1]
        filenames = ['lig.pdbqt', 'target.pdbqt']
        for ff in filenames:
            if os.path.isfile(ff):
                os.remove(ff)

    def cleanup(self):

        to_be_removed = ['run_'+self.program+'.sh', 'lig.pdbqt', 'lig_out.pdbqt', 'prepare_receptor4.log', 'target.pdbqt', 'check_lig_pdbqt.py', 'vina.err']

        # remove map files
        for filename in to_be_removed:
            if os.path.isfile(filename):
                os.remove(filename)
