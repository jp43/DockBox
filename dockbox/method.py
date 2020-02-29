import os
import sys
import stat
import shutil
import subprocess

from glob import glob

from mdkit.amber import minimization
from mdkit.utility import mol2

import configure

class DockingMethod(object):

    def __init__(self, instance, site, options):
        """Initialize docking instance"""

        self.instance = instance
        self.site = site
        self.options = options

        self.program = self.__class__.__name__.lower()

    def run_docking(self, file_r, file_l, minimize_options=None, cleanup=0, prepare_only=False, skip_docking=False):
        """Run docking one (file per ligand and receptor)"""

        curdir = os.getcwd()
        # find name for docking directory
        if 'name' in self.options:
            dockdir = self.options['name']
        else:
            dockdir = self.instance

        if self.site[0]:
            dockdir += '.' + self.site[0]

        if not skip_docking:
            # create directory for docking (remove directory if exists)
            shutil.rmtree(dockdir, ignore_errors=True)
            os.mkdir(dockdir)
        os.chdir(dockdir)

        if not skip_docking:
            print "Starting docking with %s..."%self.program.capitalize()
            print "The following options will be used:"
            options_info = ""
            for key, value in self.options.iteritems():
                options_info += str(key) + ': ' + str(value) + ', '
            print options_info[:-2]

            # (A) run docking
            script_name = "run_" + self.program + ".sh"
            self.write_docking_script(script_name, file_r, file_l)
            os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

            if prepare_only:
                return
            try:
                # try running docking procedure
                subprocess.check_output('./' + script_name + " &> " + self.program + ".log", shell=True, executable='/bin/bash')
            except subprocess.CalledProcessError as e:
                print e
                print "Error: check %s file for more details!"%(dockdir+'/'+self.program+'.log')
                os.chdir(curdir)
                return

        if prepare_only:
            return

        # (B) extract docking results
        self.extract_docking_results('score.out', file_r, file_l)

        # (C) cleanup poses (minimization, remove out-of-box poses)
        if minimize_options['minimization']:
            self.backup_files('origin')
            self.minimize_extracted_poses(file_r, 'score.out', cleanup=cleanup, **minimize_options)
        self.remove_out_of_range_poses('score.out')

        # (D) remove intermediate files if required
        if cleanup >= 1:
            self.cleanup()

        os.chdir(curdir)
        print "Docking with %s done."%self.program.capitalize()

    def run_rescoring(self, file_r, files_l):
        """Rescore multiple ligands on one receptor"""

        curdir = os.getcwd()
        # get name of rescoring from instance
        rescordir = self.instance
        if self.site[0]:
            rescordir += '.' + self.site[0]

        # overwrite previous directory if exists
        shutil.rmtree(rescordir, ignore_errors=True)
        os.mkdir(rescordir)

        # change directory
        os.chdir(rescordir)

        mol2files = files_l
        if self.program in configure.single_run_scoring_programs or (self.program == 'colvar' and self.options['type'] == 'sasa'):
            # if the program rescores in one run, provides a list of files
            mol2files = [mol2files]

        if mol2files:
            # iterate over all the poses
            for idx, file_l in enumerate(mol2files):
                # (A) write script
                script_name = "run_scoring_" + self.program + ".sh"
                self.write_rescoring_script(script_name, file_r, file_l)
                os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

                # (B) run scoring method
                try:
                    subprocess.check_output('./' + script_name + ' &> ' + self.program + '.log', shell=True, executable='/bin/bash')
                except subprocess.CalledProcessError as e:
                    print e.output
                    pass

                # (C) extract rescoring results
                if self.program in configure.single_run_scoring_programs:
                    nligands = len(file_l)
                    self.extract_rescoring_results('score.out', nligands=nligands)
                else:
                    self.extract_rescoring_results('score.out')
        else:
            # if no files provided, create an empty score.out file
            open('score.out', 'w').close()

        os.chdir(curdir)
        return scordir + '/score.out'

    def get_output_mol2files(self):
        """Get output mol2files sorted by pose ranking after docking"""

        filenames_idxs = []
        for filename in glob('lig-*.mol2'):
            suffix, ext = os.path.splitext(filename)
            filenames_idxs.append(int(suffix.split('-')[-1]))
        filenames_idxs = sorted(filenames_idxs)

        mol2files = []
        for idx in filenames_idxs:
            mol2files.append('lig-%s.mol2'%idx)
        return mol2files

    def backup_files(self, dir):
        """Do a backup of output mol2files""" 

        mol2files = self.get_output_mol2files()
        shutil.rmtree(dir, ignore_errors=True)
        os.mkdir(dir)
        for filename in mol2files:
            shutil.copyfile(filename, dir+'/'+filename) 

    def remove_scores_from_scorefile(self, file_s, indices, nligands=None):
        """Remove scores of bad poses (failed minimization, out of the box...) from score.out"""
        if os.path.exists(file_s):
            new_content = []
            with open(file_s, 'r') as sf:
                for idx, line in enumerate(sf):
                    if idx not in indices:
                        new_content.append(idx)
            if nligands:
                # consistency check
                assert nligands == idx+1, "number of ligand mol2files should be equal to number of lines in score.out"
            with open(file_s, 'w') as sf:
                for line in new_content:
                    sf.write(line)

    def minimize_extracted_poses(self, file_r, file_s, cleanup=0, **minimize_options):
        """Perform AMBER minimization on extracted poses"""

        mol2files = self.get_output_mol2files()
        if mol2files:
            # do energy minimization on ligand
            minimization.do_minimization_after_docking(file_r, mol2files, keep_hydrogens=True, charge_method=minimize_options['charge_method'],\
ncyc=minimize_options['ncyc'], maxcyc=minimize_options['maxcyc'], cut=minimize_options['cut'], amber_version=minimize_options['amber_version'])

            failed_idxs = []
            # extract results from minimization and purge out
            for idx, filename_before_min in enumerate(mol2files):
                suffix, ext = os.path.split(filename_before_min)
                filename = 'em/' + suffix + '-out' + ext
                if os.path.isfile(filename): # the minimization succeeded
                    shutil.copyfile(filename, filename_before_min)

                else: # the minimization failed
                    os.remove(filename_before_min)
                    failed_idxs.append(idx)

            # remove scores of failed poses
            self.remove_scores_from_scorefile(file_s, failed_idxs, nligands=len(mol2files))

            if failed_idxs:
                # display warning message
                failed_mol2files = [mol2files[idx] for idx in failed_idxs]
                print "Warning: minimization of poses %s failed, poses were removed!"%(', '.join(failed_mol2files))

        if cleanup >= 1:
            # if cleanup is more than 1, remove EM directory
            shutil.rmtree('em', ignore_errors=True)

    def remove_out_of_range_poses(self, file_s):
        """Get rid of poses which were predicted outside the box"""

        mol2files = self.get_output_mol2files()
        if mol2files:
            sitename, center, boxsize = self.site
            # get values of docking box center and boxsize
            center = map(float, center.split(','))
            boxsize = map(float, boxsize.split(','))
 
            out_of_range_idxs = []
            for jdx, filename in enumerate(mol2files):
                is_out = False
                for coord in mol2.get_coordinates(filename):
                    for idx, value in enumerate(coord):
                        # check if the pose is out of the box
                        if abs(value - center[idx]) > boxsize[idx]*1./2:
                            is_out = True
                            break
                    if is_out:
                        os.remove(filename)
                        out_of_range_idxs.append(jdx)
                        break
            # remove scores of failed poses
            self.remove_scores_from_scorefile(file_s, out_of_range_idxs, nligands=len(mol2files))

            if out_of_range_idxs:
                # display warning message
                out_of_range_mol2files = [mol2files[idx] for idx in out_of_range_idxs]
                print "Warning: poses %s were found out of the box, poses were removed!"%(', '.join(out_of_range_mol2files))

    def cleanup(self):
        """Remove all intermediate files"""
        for filename in glob('*'):
            if not filename.startswith('lig-') and filename != 'score.out':
                os.remove(filename)

    def write_rescoring_script(self, script_name, file_r, file_l):
        pass

    def extract_rescoring_results(self, filename):
        pass

    def write_docking_script(self, script_name, file_r, file_l):
        pass

    def extract_docking_results(self, file_r, file_l, file_s, input_file_r):
        pass

class ScoringMethod(DockingMethod):

    def run_docking(self, file_r, file_l, minimize=False, cleanup=0, extract_only=False):
        pass

    def remove_out_of_range_poses(self, file_s):
        pass

    def minimize_extracted_poses(self, file_r):
        pass
