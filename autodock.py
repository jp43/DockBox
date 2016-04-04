import os
import sys
import tempfile
import subprocess
import glob
import shutil
import fileinput
import numpy as np

import tools.PDB as pdbt
import tools.mol2 as mol2t

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'prepare_dpf4.py', 'prepare_gpf4.py', 'autogrid4', 'autodock4', 'babel']

default_settings = {'ga_run': '100', 'spacing': '0.238'}

autogrid_options_names = ['spacing']
autodock_options_names = ['ga_run', 'ga_pop_size', 'ga_num_evals', 'ga_num_generations', 'outlev', 'seed']

def write_docking_script(filename, input_file_r, input_file_l, site, options, rescoring=False):

    # set box center
    center = site[1] # set box
    gridcenter = '\"' + ','.join(map(str.strip, center.split(','))) + '\"'

    # set box size
    boxsize = site[2]
    boxsize = map(float, map(str.strip, boxsize.split(',')))
    spacing = float(options['spacing'])
    npts = []

    for size in boxsize:
        #if rescoring:
        #    sz = int((size*1.0 + 5)/spacing) # add 5A to avoid problems at boundaries
        sz = int(size*1.0/spacing) + 1
        npts.append(str(sz)) # round to the integer above
    npts =  ','.join(npts)

    prepare_file_ligand(input_file_l, 'lig.mol2', rescoring)

    autogrid_options = {}
    for name in autogrid_options_names:
        if name in options:
            autogrid_options[name] = options[name]
    
    # create flag with specified options for autogrid
    autogrid_options_flag = "-p npts=%(npts)s -p gridcenter=%(gridcenter)s "% locals()
    autogrid_options_flag += ' '.join(['-p ' + key + '=' + value for key, value in autogrid_options.iteritems()])

    autodock_options = {}
    for name in autodock_options_names:
        if name in options:
            autodock_options[name] = options[name]

    # create flag with specified options for autodock
    autodock_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in autodock_options.iteritems()])

    if not rescoring:
        if 'ga_num_evals' not in options:
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
prepare_ligand4.py -l lig.mol2 -U '' -C -B 'amide_guadinidium' -o lig.pdbqt
prepare_receptor4.py -r %(input_file_r)s -U 'waters_nonstdres' -o target.pdbqt

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
prepare_ligand4.py -l lig.mol2 -U '' -C -B 'amide_guadinidium' -o lig.pdbqt
if [ ! -f target.pdbqt ]; then
  prepare_receptor4.py -r %(input_file_r)s -U 'waters_nonstdres' -o target.pdbqt
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


def extract_docking_results(file_r, file_l, file_s, input_file_r, extract):

    for ff in [file_s, file_l]:
        if os.path.isfile(ff):
            os.remove(ff)

    shutil.copyfile(input_file_r, file_r)

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

        if extract == 'lowest':
            cluster_idxs = [np.argmax(np.array(hist))+1]
        elif extract == 'all': 
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

    if extract in ['lowest', 'all']:
        cleanup_pose(file_l)

def prepare_file_ligand(input_file_l, output_file_l, rescoring):

    if rescoring:
        mol2t.pdb2mol2(input_file_l, output_file_l, sample='../../autodock.site2/lig.mol2')
    else:
        # convert .sdf file to mol2
        subprocess.check_call('babel -isdf %s -omol2 %s 2>/dev/null'%(input_file_l, output_file_l), shell=True)

        # give unique atom names
        mol2t.give_unique_atom_names(output_file_l)

def write_rescoring_script(filename, file_r, file_l, site, options):

    write_docking_script(filename, file_r, file_l, site, options, rescoring=True)

def extract_rescoring_results(filename):

    with open(filename, 'a') as ff:
        with open('dock.dlg', 'r') as outf:
            for line in outf:
                if line.startswith('epdb: USER    Estimated Free Energy of Binding'):
                    print >> ff, line.split()[8]

    os.remove('lig.pdbqt')
    os.remove('target.pdbqt')

def cleanup_pose(file_l):

    # rearrange atoms in case if the original order was modified
    atoms_names = mol2t.get_atoms_names('lig.mol2')
    pdbt.rearrange_atom_names(file_l, atoms_names)

    # remove/add hydrogens
    #pdbt.format_lig_file('lig.out.pdb')
    #subprocess.check_call('babel -ipdb lig.out.pdb -opdb ligtmp.out.pdb -d -h 2>/dev/null', shell=True)
    #shutil.move('ligtmp.out.pdb','lig.out.pdb')

    ## give unique atom names
    #pdbt.give_unique_atom_names('lig.out.pdb')

def cleanup():
    # remove map files
    for ff in glob.glob('*map*'):
        os.remove(ff)

    for ff in glob.glob('*pdbqt'):
        os.remove(ff)
