import os
import sys
import tempfile
import subprocess
import glob
import shutil
import fileinput
import numpy as np

import util.pdbtools as pdbt

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'prepare_dpf4.py', 'prepare_gpf4.py', 'autogrid4', 'autodock4', 'babel']

default_settings = {'ga_run': '100', 'spacing': '0.238'}

autogrid_options_names = ['npts', 'spacing', 'gridcenter']
autodock_options_names = ['ga_run', 'ga_pop_size', 'ga_num_evals', 'ga_num_generations', 'outlev']

def set_site_options(config):

    # set box center
    center = config.site['center'] # set box
    config.options['autodock']['gridcenter'] = '\"' + ','.join(map(str.strip, center.split(','))) + '\"'

    # set box size
    boxsize = config.site['boxsize']
    boxsize = map(float, map(str.strip, boxsize.split(',')))
    spacing = float(config.options['autodock']['spacing'])
    npts = []
    for size in boxsize:
         npts.append(str(int(size/spacing)))
    config.options['autodock']['npts'] =  ','.join(npts)

def write_docking_script(filename, input_file_r, input_file_l, config):

    prepare_pdbfile(input_file_l, 'lig.pdb', config)

    autogrid_options = {}
    for name in autogrid_options_names:
        if name in config.options['autodock']:
            autogrid_options[name] = config.options['autodock'][name]
    
    # create flag with specified options for autogrid
    autogrid_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in autogrid_options.iteritems()])

    autodock_options = {}
    for name in autodock_options_names:
        if name in config.options['autodock']:
            autodock_options[name] = config.options['autodock'][name]

    # create flag with specified options for autodock
    autodock_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in autodock_options.iteritems()])

    if 'ga_num_evals' not in config.options['autodock']:
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
prepare_ligand4.py -l lig.pdb -o lig.pdbqt
prepare_receptor4.py -r %(input_file_r)s -o target.pdbqt

# run autogrid
prepare_gpf4.py -l lig.pdbqt -r target.pdbqt -o grid.gpf %(autogrid_options_flag)s
autogrid4 -p grid.gpf -l grid.glg

# run autodock
%(ga_num_evals_lines)s
prepare_dpf4.py -l lig.pdbqt -r target.pbdqt -o dock.dpf -p move=lig.pdbqt %(autodock_options_flag)s $ga_num_evals_flag
autodock4 -p dock.dpf -l dock.dlg"""% locals()
        file.write(script)

def extract_docking_results(file_r, file_l, file_s, config):

    for ff in [file_s, file_l]:
        if os.path.isfile(ff):
            os.remove(ff)

    shutil.copyfile(config.input_file_r, file_r)

    hist = []
    with open('dock.dlg', 'r') as dlgfile:
        # (A) get the index of the most populated cluster
        line = dlgfile.next()
        while "CLUSTERING HISTOGRAM" not in line:
            line = dlgfile.next()
        nlines_to_skip = 8
        for idx in range(nlines_to_skip):
            dlgfile.next()
        while True:
            line = dlgfile.next()
            if line[0] == '_':
                break
            hist.append(int(line.split('|')[4].strip()))

        if config.extract == 'lowest':
            cluster_idxs = [np.argmax(np.array(hist))+1]
        elif config.extract == 'all': 
            cluster_idxs = (np.argsort(np.array(hist))+1).tolist()
            cluster_idxs =  cluster_idxs[::-1]

    for cluster_idx in cluster_idxs:
        with open('dock.dlg', 'r') as dlgfile:
            # (B) get the lowest binding free energy
            while "Cluster Rank = %i"%cluster_idx not in line:
                oldline = line # do backup 
                line = dlgfile.next()
            run_idx = int(oldline.split('=')[1])
            nlines_to_skip = 4
            for idx in range(nlines_to_skip):
                dlgfile.next()
            line = dlgfile.next()
            score = float(line[47:53])

            # save the binding free energy
            with open(file_s, 'a') as sf:
                print >> sf, score

            # (C) save the correct pose
            while "ATOM" not in line:
                line = dlgfile.next()

            with open(file_l, 'a') as lf:
                while "TER" not in line:
                    print >> lf, 'ATOM  ' + line[6:66]
                    line = dlgfile.next()
                print >> lf, "ENDMDL"

    if config.extract in ['lowest', 'all']:
        cleanup_pose(file_l, config)

def prepare_pdbfile(file_l, pdbfile, config):

    # convert sdffile to PDBfile
    subprocess.check_call('babel -isdf %s -opdb %s 2>/dev/null'%(file_l,pdbfile), shell=True)

    # give unique atom names
    pdbt.give_unique_atom_names(pdbfile)

def cleanup_pose(file_l, config):

    atom_names = []
    with open('lig.pdb', 'r') as pdbfile:
        for line in pdbfile:
            if line.startswith(('ATOM', 'HETATM')):
                atom_names.append(line[12:17].strip())

    # rearrange atoms in case if the original order was modified
    pdbt.rearrange_atom_names(file_l, atom_names)

    # remove/add hydrogens
    pdbt.format_lig_file('lig.out.pdb')
    subprocess.check_call('babel -ipdb lig.out.pdb -opdb ligtmp.out.pdb -d -h 2>/dev/null', shell=True)
    shutil.move('ligtmp.out.pdb','lig.out.pdb')

    # give unique atom names
    pdbt.give_unique_atom_names('lig.out.pdb')

def cleanup(config):
    # remove map files
    for ff in glob.glob('*map*'):
        os.remove(ff)

    for ff in glob.glob('*pdbqt'):
        os.remove(ff)
