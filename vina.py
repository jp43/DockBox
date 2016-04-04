import os
import sys
import glob
import autodock
import shutil

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'vina', 'babel']

default_settings = {'cpu': '1'}

def write_docking_script(filename, input_file_r, input_file_l, site, options, rescoring=False):

    center = site[1]
    center = map(str.strip, center.split(','))

    boxsize = site[2]
    boxsize = map(str.strip, boxsize.split(','))

    autodock.prepare_file_ligand(input_file_l, 'lig.mol2', rescoring)

    # write vina config file
    with open('vina.config', 'w') as config_file:
        print >> config_file, 'receptor = target.pdbqt'
        print >> config_file, 'ligand = lig.pdbqt'
        print >> config_file, 'center_x = ' + center[0]
        print >> config_file, 'center_y = ' + center[1]
        print >> config_file, 'center_z = ' + center[2]
        print >> config_file, 'size_x = ' + float(boxsize[0])
        print >> config_file, 'size_y = ' + float(boxsize[1])
        print >> config_file, 'size_z = ' + float(boxsize[2])
        for key, value in options.iteritems():
            print >> config_file, key + ' = ' + value

    # write vina script
    if not rescoring:
        with open(filename, 'w') as file:
            script ="""#!/bin/bash
set -e
# generate .pdbqt files
prepare_ligand4.py -l lig.mol2 -C -B 'amide_guadinidium' -o lig.pdbqt
prepare_receptor4.py -r %(input_file_r)s -U 'waters_nonstdres' -o target.pdbqt

# run vina
vina --config vina.config > vina.out"""% locals()
            file.write(script)
    else:
        with open(filename, 'w') as file:
            script ="""#!/bin/bash
set -e
# generate .pdbqt files
prepare_ligand4.py -l lig.pdb -o lig.pdbqt
prepare_receptor4.py -r %(input_file_r)s -o target.pdbqt

# run vina
vina --score_only --config vina.config > vina.out"""% locals()
            file.write(script)

def extract_docking_results(file_r, file_l, file_s, input_file_r, site, extract):

    shutil.copyfile(input_file_r, file_r)

    # get center and box size to remove those poses that are out of the box
    center = map(float, site[1].split(',')) 
    boxsize = map(float, site[2].split(','))

    if extract in ['lowest', 'all']:
        with open('lig_out.pdbqt','r') as pdbqtfile:
            with open(file_l, 'w') as lf:
                with open(file_s, 'w') as sf:
                    coords = []
                    for line in pdbqtfile:
                        if line.startswith(('ATOM', 'HETATM')):
                            coords.append('ATOM  ' + line[6:66])
                        elif line.startswith('REMARK VINA RESULT:'):
                            score = float(line[19:].split()[0])
                        elif line.startswith('ENDMDL'):
                            # check if the pose is out of the box
                            isout = False # first suppose that it is not out
                            for coord in coords:
                                for idx, xyz in enumerate([coord[30:38], coord[38:46], coord[46:54]):
                                    if abs(float(xyz)-center[idx]) > boxsize[idx]*1./2: # the pose is out of the box
                                        isout = True
                                        break
                            if not isout:
                                for coord in coords:
                                    print >> lf, coord
                                print >> lf, 'ENDMDL'
                                print >> sf, score
                                if extract == 'lowest':
                                    break
                            coords = [] # reinitialize the coordinates
    autodock.cleanup_pose(file_l)

def write_rescoring_script(filename, file_r, file_l, site, options):

    write_docking_script(filename, file_r, file_l, site, options, rescoring=True)

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
