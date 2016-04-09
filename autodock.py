import os
import sys
import tempfile
import subprocess
import glob
import shutil
import fileinput
import method
import numpy as np

import amber.minimz as mn
import tools.PDB as pdbt
import tools.mol2 as mol2t

required_programs = ['prepare_ligand4.py', 'prepare_receptor4.py', 'prepare_dpf4.py', 'prepare_gpf4.py', 'autogrid4', 'autodock4', 'babel']

default_settings = {'ga_run': '100', 'spacing': '0.238'}

class ADBased(method.DockingMethod):

    def prepare_ligand(self, input_file_l, mol2file, rescoring):
        """Prepare ligand structure for Autodock based docking methods, i.e. convert initial file to mol2 format"""

        #if rescoring:
        #    mol2t.pdb2mol2(input_file_l, output_file_l, sample='../../autodock.site3/lig.mol2')
        if not rescoring:
            # convert .sdf file to mol2
            subprocess.check_call('babel -isdf %s -omol2 %s 2>/dev/null'%(input_file_l, mol2file), shell=True)

            # give unique atom names
            mol2t.give_unique_atom_names(mol2file)

    def fix_poses(self, input_file_r, output_pdbfiles_l):
        """Convert output PDB files to mol2 and perform energy minimization on non-polar hydrogens"""

        input_mol2file = 'lig.mol2'
        # rearrange atoms in case the original order was modified
        atoms_names = mol2t.get_atoms_names(input_mol2file)

        n_output_files = len(output_pdbfiles_l)

        mol2files = []
        for idx, file_l in enumerate(output_pdbfiles_l):
            pdbt.rearrange_atom_names(file_l, atoms_names)
            mol2file = 'lig-%s.out.mol2'%idx
            mol2t.update_mol2_from_pdb(file_l, mol2file, sample_mol2file=input_mol2file) 
            #os.remove(file_l)
            mol2files.append(mol2file)

        # do energy minimization on ligand hydrogens
        mn.do_minimization(input_file_r, files_l=mol2files, restraints=":LIG & @H=", keep_hydrogens=True)

        # extract results from minimization and purge out
        new_poses = []
        for idx in range(n_output_files):
            mol2file = 'lig-%s.out.mol2'%idx
            #shutil.move('minimz/' + mol2file,  mol2file)
            new_poses.append(mol2file)
        #shutil.rmtree('minimz')

        return new_poses

    def write_docking_script(filename, file_r, file_l, rescoring=True):
        pass

    def write_rescoring_script(self, filename, file_r, file_l):
        self.write_docking_script(filename, file_r, file_l, rescoring=True)

class Autodock(ADBased):

    def __init__(self, name, site, options):

        super(Autodock, self).__init__(name, site, options)

        # set box center
        self.options['gridcenter'] = '\"' + ','.join(map(str.strip, site[1].split(','))) + '\"'
 
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

        self.prepare_ligand(file_l, 'lig.mol2', rescoring)

        # create flags with specified options for autogrid and autodock
        autogrid_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in self.autogrid_options.iteritems()])
        autodock_options_flag = ' '.join(['-p ' + key + '=' + value for key, value in self.autodock_options.iteritems()])

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
prepare_ligand4.py -l lig.mol2 -C -B 'amide_guadinidium' -o lig.pdbqt
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
prepare_ligand4.py -l lig.mol2 -C -B 'amide_guadinidium' -o lig.pdbqt
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

    def extract_docking_results(self, file_r, file_l, file_s, input_file_r, extract):
    
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

        output_pdbfiles_l = []
        for kdx, cluster_idx in enumerate(cluster_idxs):
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

                pdbfile = 'lig-%s.out.pdb'%kdx
                with open(pdbfile, 'w') as lf:
                    while "TER" not in line:
                        print >> lf, 'ATOM  ' + line[6:66]
                        line = dlgfile.next()
                output_pdbfiles_l.append(pdbfile)
 
        if extract in ['lowest', 'all']:
            self.fix_poses(input_file_r, output_pdbfiles_l)

    def extract_rescoring_results(self, filename):
    
        with open(filename, 'a') as ff:
            with open('dock.dlg', 'r') as outf:
                for line in outf:
                    if line.startswith('epdb: USER    Estimated Free Energy of Binding'):
                        print >> ff, line.split()[8]
    
        os.remove('lig.pdbqt')
        os.remove('target.pdbqt')
    
    def cleanup(self):
        # remove map files
        for ff in glob.glob('*map*'):
            os.remove(ff)
    
        for ff in glob.glob('*pdbqt'):
            os.remove(ff)
