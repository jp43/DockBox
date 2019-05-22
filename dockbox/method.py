import os
import sys
import stat
from glob import glob
import shutil
import subprocess
import setconf

from mdkit.amber import minimization
from mdkit.utility import mol2
from mdkit.amber import clustering

class DockingMethod(object):

    def __init__(self, instance, site, options):

        self.instance = instance
        self.site = site
        self.options = options

        self.program = self.__class__.__name__.lower()

    def run_docking(self, file_r, file_l, minimize_options=None, cleanup=0, cutoff_clustering=0.0, prepare_only=False, skip_docking=False):
        """Run docking on one receptor (file_r) and one ligand (file_l)"""

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
                print "Check %s file for more details!"%(dockdir+'/'+self.program+'.log')
                pass

        if prepare_only:
            return

        # (B) extract docking results
        self.extract_docking_results('score.out', file_r, file_l)

        # (C) cleanup poses (minimization, remove out-of-box poses)
        if minimize_options['minimization']:
            self.backup_files('origin')
            self.minimize_extracted_poses(file_r, 'score.out', cleanup=cleanup, **minimize_options)
        self.remove_out_of_range_poses('score.out')

        if cutoff_clustering != 0.0:
            self.remove_duplicates('score.out', cutoff=cutoff_clustering)

        # (D) remove intermediate files if required
        if cleanup >= 1:
            self.cleanup()

        os.chdir(curdir)
        print "Docking with %s done."%self.program.capitalize()

    def run_rescoring(self, file_r, files_l):
        """Rescore multiple ligands on one receptor"""

        curdir = os.getcwd()
        # find name for scoring directory
        if 'name' in self.options:
            scordir = self.options['name']
        else:
            scordir = self.instance

        if self.site[0]:
            scordir += '.' + self.site[0]
        shutil.rmtree(scordir, ignore_errors=True)
        os.mkdir(scordir)

        # change directory
        os.chdir(scordir)

        if self.program in setconf.single_run_scoring_programs or (self.program == 'colvar' and self.options['type'] == 'sasa'):
            # if the program rescores in one run, provides a list of files
            files_l = [files_l]

        if files_l:
            # iterate over all the poses
            for idx, file_l in enumerate(files_l):
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

                # (C) extract docking results
                if self.program in setconf.single_run_scoring_programs:
                    nligands = len(files_l[0])
                    self.extract_rescoring_results('score.out', nligands=nligands)
                else:
                    self.extract_rescoring_results('score.out')
        else:
            # if no files provided, create an empty score.out file
            open('score.out', 'w').close()

        os.chdir(curdir)
        return scordir + '/score.out'


    def get_output_files_l(self):

        poses_idxs = []
        for filename in glob('lig-*.mol2'):
            poses_idxs.append(int((filename.split('.')[-2]).split('-')[-1]))
        poses_idxs = sorted(poses_idxs)

        files_l = []
        for pose_idx in poses_idxs:
            files_l.append('lig-%s.mol2'%pose_idx)

        return files_l

    def backup_files(self, dir):

        files_l = self.get_output_files_l()
        shutil.rmtree(dir, ignore_errors=True)
        os.mkdir(dir)
        for file_l in files_l:
            shutil.copyfile(file_l, dir+'/'+file_l) 

    def minimize_extracted_poses(self, file_r, file_s, cleanup=0, **minimize_options):
        """Perform AMBER minimization on extracted poses"""

        files_l = self.get_output_files_l()
        nfiles_l = len(files_l)

        # get minimization options
        charge_method = minimize_options['charge_method']
        ncyc = minimize_options['ncyc']
        maxcyc = minimize_options['maxcyc']
        cut = minimize_options['cut']
        amber_version = minimize_options['amber_version']

        if files_l:
            # do energy minimization on ligand
            minimization.do_minimization_after_docking(file_r, files_l, keep_hydrogens=True, charge_method=charge_method,\
ncyc=ncyc, maxcyc=maxcyc, cut=cut, amber_version=amber_version)

        failed_idxs = []
        # extract results from minimization and purge out
        for idx in range(nfiles_l):
            mol2file = 'lig-%s-out.mol2'%(idx+1)
            if os.path.isfile('em/'+mol2file): # the minimization succeeded
                shutil.copyfile('em/'+mol2file, 'lig-%s.mol2'%(idx+1))
            else: # the minimization failed
                os.remove('lig-%s.mol2'%(idx+1))
                failed_idxs.append(idx)

        if files_l:
            if os.path.exists(file_s):
                with open(file_s, 'r') as sf:
                    with open('score.tmp.out', 'w') as sft:
                        for idx, line in enumerate(sf):
                            if idx not in failed_idxs:
                                sft.write(line)
                shutil.move('score.tmp.out', file_s)

        if cleanup >= 1:
            shutil.rmtree('em', ignore_errors=True)

    def remove_out_of_range_poses(self, file_s):
        """Get rid of poses which were predicted outside the box"""

        files_l = self.get_output_files_l()

        center = map(float, self.site[1].split(','))
        boxsize = map(float, self.site[2].split(','))

        out_of_range_idxs = []
        for jdx, file_l in enumerate(files_l):
            isout = False
            coords = mol2.get_coordinates(file_l)
            for kdx, coord in enumerate(coords):
                for idx, xyz in enumerate(coord):
                    # check if the pose is out of the box
                    if abs(float(xyz)-center[idx]) > boxsize[idx]*1./2:
                        isout = True
                        break
            if isout:
                #print file_l, "out"
                os.remove(file_l)
                out_of_range_idxs.append(jdx)

        if files_l:
            if os.path.exists(file_s):
                with open(file_s, 'r') as sf:
                    with open('score.tmp.out', 'w') as sft:
                        for idx, line in enumerate(sf):
                            if idx not in out_of_range_idxs:
                                sft.write(line)
                shutil.move('score.tmp.out', file_s)

    def remove_duplicates(self, file_s, cutoff=0.0):

        files_l = self.get_output_files_l()
        files_r = [file_r for idx in range(len(files_l))]

        nfiles_l = len(files_l)
        if nfiles_l > 1:
            # cluster poses
            clustering.cluster_poses(files_r, files_l, cutoff=cutoff, cleanup=True)

            with open('clustering/info.dat', 'r') as ff:
                for line in ff:
                    if line.startswith('#Representative frames:'):
                        rep_structures = map(int, line.split()[2:])

            for idx, file_l in enumerate(files_l):
                if idx+1 not in rep_structures:
                    os.remove(file_l)

            if os.path.exists(file_s):
                with open(file_s, 'r') as sf:
                    with open('score.tmp.out', 'w') as sft:
                        for idx, line in enumerate(sf):
                            if idx+1 in rep_structures:
                                sft.write(line)

                shutil.move('score.tmp.out', file_s)

    def write_rescoring_script(self, script_name, file_r, file_l):
        pass

    def extract_rescoring_results(self, filename):
        pass

    def write_docking_script(self, script_name, file_r, file_l):
        pass

    def extract_docking_results(self, file_r, file_l, file_s, input_file_r):
        pass

    def cleanup(self):
        pass

class ScoringMethod(DockingMethod):

    def run_docking(self, file_r, file_l, minimize=False, cleanup=0, extract_only=False):
        pass

    def remove_out_of_range_poses(self, file_s):
        pass

    def minimize_extracted_poses(self, file_r):
        pass
