import os
import sys
import shutil
import subprocess
from glob import glob

from mdkit.utility import mol2
import autodock

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'vina', 'obabel']

default_settings = {'cpu': '1', 'num_modes': '9', 'energy_range': '3', 'seed': None}

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
            cf.write('receptor = target.pdbqt\n')
            cf.write('ligand = ligand.pdbqt\n')
            # write other options
            for key, value in self.options.iteritems():
                if value is not None:
                    cf.write(key+' = '+value+'\n')

        # write vina script
        if not rescoring:
            with open(filename, 'w') as ff:
                script ="""#!/bin/bash
set -e

MGLPATH=`which prepare_ligand4.py`
MGLPATH=`python -c "print '/'.join('$MGLPATH'.split('/')[:-3])"`
export PYTHONPATH=$PYTHONPATH:$MGLPATH

# prepare ligand
prepare_ligand4.py -l %(file_l)s -o ligand.pdbqt
python check_ligand_pdbqt.py ligand.pdbqt

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
prepare_ligand4.py -l %(file_l)s -o ligand.pdbqt
python check_ligand_pdbqt.py ligand.pdbqt

if [ ! -f target.pdbqt ]; then
  prepare_receptor4.py -U nphs_lps_waters -r %(file_r)s -o target.pdbqt > prepare_receptor4.log
  python check_ions.py target.pdbqt prepare_receptor4.log
fi

# run vina
vina --score_only --config vina.config > vina.out"""% locals()
                ff.write(script)

    def extract_docking_results(self, file_s, input_file_r, input_file_l):
        """Extract output structures in .mol2 formats"""

        poses_extracted = False
        if os.path.exists('ligand_out.pdbqt'):
            try:
                subprocess.check_output('obabel -ipdbqt ligand_out.pdbqt -omol2 -Opose-.mol2 -m &>/dev/null', shell=True, executable='/bin/bash')
                self.update_output_mol2files(sample=input_file_l)
                poses_extracted = True
            except:
                for mol2file in glob('pose-*.mol2'):
                    os.remove(mol2file)
                poses_extracted = False

        if poses_extracted:
            with open('ligand_out.pdbqt','r') as dlgf:
                with open(file_s, 'w') as sf:
                    for line in dlgf:
                        if line.startswith('REMARK VINA RESULT:'):
                            score = line[19:].split()[0].strip()
                            sf.write(score+'\n')
        else:
            open(file_s, 'w').close()

    def write_rescoring_script(self, filename, file_r, file_l):
        self.write_docking_script(filename, file_r, file_l, rescoring=True)
    
    def extract_rescoring_results(self, filename):
        with open(filename, 'a') as ff:
            with open('vina.out', 'r') as outf:
                for line in outf:
                    if line.startswith('Affinity:'):
                        ff.write(line.split()[1].strip()+'\n')
        filenames = ['ligand.pdbqt', 'target.pdbqt']
        for ff in filenames:
            if os.path.isfile(ff):
                os.remove(ff)
