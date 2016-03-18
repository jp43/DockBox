import os
import sys
import glob
import autodock
import shutil

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'vina', 'babel']

default_settings = {'cpu': '1'}

known_settings = {'herg': {'center_x': '3.966', 'center_y': '8.683', 'center_z': '11.093', 'size_x': '30.0', 'size_y': '30.0', 'size_z': '30.0'}, \
'herg-cut': {'center_x': '3.966', 'center_y': '8.683', 'center_z': '11.093', 'size_x': '30.0', 'size_y': '30.0', 'size_z': '30.0'}, \
'herg-inactivated': {'center_x': '0.000', 'center_y': '0.000', 'center_z': '-5.000', 'size_x': '30.0', 'size_y': '30.0', 'size_z': '30.0'}}

required_settings_names = ['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z']

def write_docking_script(filename, input_file_r, input_file_l, config):

    autodock.prepare_pdbfile(input_file_l, 'lig.pdb', config)

    # write vina config file
    with open('vina.config', 'w') as config_file:
        print >> config_file, 'receptor = target.pdbqt'
        print >> config_file, 'ligand = lig.pdbqt'
        for key, value in config.options['vina'].iteritems():
            print >> config_file, key + ' = ' + value

    # write vina script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash
set -e

# generate .pdbqt files
prepare_ligand4.py -l lig.pdb -o lig.pdbqt
prepare_receptor4.py -r %(input_file_r)s -o target.pdbqt

# run vina
vina --config vina.config > vina.out"""% locals()
        file.write(script)

def extract_docking_results(file_r, file_l, file_s, config):

    shutil.copyfile(config.input_file_r, file_r)

    if config.extract in ['lowest', 'all']:
        with open('lig_out.pdbqt','r') as pdbqtfile:
            with open(file_l, 'w') as lf:
                with open(file_s, 'w') as sf:
                    for line in pdbqtfile:
                        if line.startswith(('ATOM', 'HETATM')):
                            print >> lf, 'ATOM  ' + line[6:66]
                        elif line.startswith('REMARK VINA RESULT:'):
                            score = float(line[19:].split()[0])
                            print >> sf, score
                        elif line.startswith('ENDMDL'):
                            print >> lf, line.replace('\n','')
                            if config.extract == 'lowest':
                                break
    autodock.cleanup_pose(file_l, config)

def cleanup(config):

    for ff in glob.glob('*pdbqt'):
        os.remove(ff)
