import os
import sys
import glob
import autodock
import shutil

import util.pdbtools as pdbt

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'vina', 'babel']

default_settings = {'cpu': '1'}

def set_site_options(site, options):

    # set box center
    center = site['center'] # set box
    center = map(str.strip, center.split(','))

    boxsize = site['boxsize']
    boxsize = map(str.strip, boxsize.split(','))
    xyz = ['x', 'y', 'z']

    for idx, axis in enumerate(xyz):
        options['center_'+axis] = center[idx]
        options['size_'+axis] = boxsize[idx]

def write_docking_script(filename, input_file_r, input_file_l, options, rescoring=False):

    if not rescoring:
        autodock.prepare_pdbfile(input_file_l, 'lig.pdb')

    # write vina config file
    with open('vina.config', 'w') as config_file:
        print >> config_file, 'receptor = target.pdbqt'
        print >> config_file, 'ligand = lig.pdbqt'
        for key, value in options.iteritems():
            print >> config_file, key + ' = ' + value

    if rescoring:
        rescoring_flag = '--score_only'
    else:
        rescoring_flag = ''

    # write vina script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash
set -e

# generate .pdbqt files
prepare_ligand4.py -l lig.pdb -o lig.pdbqt
prepare_receptor4.py -r %(input_file_r)s -o target.pdbqt

# run vina
vina %(rescoring_flag)s --config vina.config > vina.out"""% locals()
        file.write(script)

def extract_docking_results(file_r, file_l, file_s, input_file_r, extract):

    shutil.copyfile(input_file_r, file_r)

    if extract in ['lowest', 'all']:
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
                            if extract == 'lowest':
                                break
    autodock.cleanup_pose(file_l)

def write_rescoring_script(filename, file_r, file_l, options):

    write_docking_script(filename, file_r, file_l, options, rescoring=True)

def extract_rescoring_results(filename):

    with open(filename, 'a') as ff:
        with open('vina.out', 'r') as outf:
            for line in outf:
                if line.startswith('Affinity:'):
                    print >> ff, line.split()[1]

    os.remove('lig.pdbqt')
    os.remove('target.pdbqt')
    
def cleanup():

    for ff in glob.glob('*pdbqt'):
        os.remove(ff)
