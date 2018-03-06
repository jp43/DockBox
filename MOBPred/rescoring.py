import os
import sys
import shutil
import time
import subprocess

import setup

class Rescoring(object):

    def __init__(self, config, args):

        self.is_rescoring = self.is_yesno_option(config, 'DOCKING', 'rescoring')

        if self.is_rescoring:
            self.config = setup.ScoringSetup(config)
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

    def run(self, file_r, posedir):
        """Run rescoring on docking poses"""
    
        tcpu1 = time.time()
        print "Starting rescoring..."
    
        # look for results folder
        if not os.path.isdir(posedir):
            raise IOError('no folder %s found!'%posedir)
        else:
            with open(posedir+'/info.dat') as inff:
                nposes = inff.next()
                nposes = nposes[1:] # the first character is a # sign
                nposes = map(int, nposes.split(','))

        curdir = os.getcwd()
        workdir = 'rescoring'
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        os.chdir(workdir)
        # iterate over rescoring instances
        for instance, program, options in self.instances:

            # possibility of renaming the folder and output file 
            if 'name' in options:
                name = options['name']
            else:
                name = instance

            # remove old scoring file
            if os.path.isfile(name+'.score'):
                os.remove(name+'.score')

            for kdx in range(len(self.site)):
                site = self.site['site'+str(kdx+1)]

                # get complex filenames
                files_l = [os.path.abspath('../'+posedir+'/lig-%s.mol2'%idx) for idx in range(nposes[kdx], nposes[kdx+1])]

                # get docking class
                DockingClass = getattr(sys.modules[program], program.capitalize())

                DockingInstance = DockingClass(instance, site, options)
                outputfile = DockingInstance.run_rescoring(file_r, files_l)

                # cat output in file (cat instead of copying because of the binding sites)
                subprocess.check_output('cat %s >> %s'%(outputfile,name+'.score'), shell=True, executable='/bin/bash')

        os.chdir(curdir)
        tcpu2 = time.time()
        print "Rescoring done. Total time needed: %i s" %(tcpu2-tcpu1)
