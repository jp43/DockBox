import sys
import os
import shutil
import time
import subprocess

from DockTbx import multi

class Rescoring(object):

    def __init__(self, config, args):

        self.is_rescoring = self.is_yesno_option(config, 'DOCKING', 'rescoring')

        if self.is_rescoring:
            self.config = multi.MultiProgramScoring(config)
            self.instances = self.config.instances
            self.site = self.config.site

        self.rescore_only = False
        if args.rescore_only:
            self.rescore_only = True

    def is_yesno_option(self, config, section, option, default=False):

        if config.has_option(section, option):
            yesno = config.get(section, option).lower()
            if yesno == 'yes':
                return True
            elif yesno == 'no':
                return False
            else:
                raise ValueError("option %s should be yes or no!"%option)
        else:
            return default

    def run(self, file_r):
        """Run rescoring on docking poses"""
    
        tcpu1 = time.time()
        print "Starting rescoring..."
    
        # look for results folder
        if not os.path.isdir('poses'):
            raise IOError('no folder poses found!')
        else:
            with open('poses/info.dat') as inff:
                nposes = inff.next()
                nposes = map(int, nposes.strip().split())

        curdir = os.getcwd()
        workdir = 'rescoring'

        shutil.rmtree(workdir, ignore_errors=True)
        os.mkdir(workdir)
        os.chdir(workdir)
    
        for kdx in range(len(self.site)):
            site = self.site['site'+str(kdx+1)]

            # iterate over rescoring instances
            for instance, program, options in self.instances:
                # get complex filenames
                files_l = [os.path.abspath('../poses/lig-%s.mol2'%idx) for idx in range(nposes[kdx], nposes[kdx+1])]

                # get docking class
                DockingClass = getattr(sys.modules[program], program.capitalize())

                DockingInstance = DockingClass(instance, site, options)
                outputfile = DockingInstance.run_rescoring(file_r, files_l)

                # cat output in file
                subprocess.check_output('cat %s >> %s'%(outputfile,instance+'.score'), shell=True, executable='/bin/bash')

        os.chdir(curdir)
        tcpu2 = time.time()
        print "Rescoring done. Total time needed: %i s" %(tcpu2-tcpu1)
