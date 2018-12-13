#!/usr/bin/python
from __future__ import with_statement
import os
import sys
import shutil
import argparse
import ConfigParser
import time
import pandas as pd
from glob import glob
import subprocess
import setup
from mdtools.utility import mol2

class DockingConfig(object):

    def __init__(self, args, task='docking'):

        # check if config file exist
        if not os.path.exists(args.config_file):
            raise ValueError("Config file %s not found!"%(args.config_file))

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file)

        # prepare ligand file
        file_l = os.path.abspath(args.input_file_l)
        new_file_l = os.path.basename(file_l)
        pref, ext = os.path.splitext(new_file_l)
        new_file_l = pref + '_u' + ext # new ligand file with unique names for every atom

        # create a ligand file with unique atom names
        mol2.update_mol2file(file_l, new_file_l, unique=True, ligname='LIG')
        self.input_file_l = os.path.abspath(new_file_l)

        # check if ligand file exists
        if not os.path.exists(self.input_file_l):
            raise IOError("File %s not found!"%(self.input_file_l))

        self.input_file_r = os.path.abspath(args.input_file_r)

        # check if receptor file exists
        if not os.path.exists(self.input_file_r):
            raise IOError("File %s not found!"%(self.input_file_r))

        if task == 'docking':
            self.docking = setup.DockingSetup(config)
            self.rescoring = setup.RescoringSetup(config)
        elif task == 'scoring':
            self.scoring = setup.ScoringSetup(config)
        else:
            raise ValueError("Task should be one of docking or scoring")

class Scoring(object):

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="""runscore : score with multiple software --------
Requires one file for the ligand (1 struct.) and one file for the receptor (1 struct.)""")

        parser.add_argument('-l',
            type=str,
            dest='input_file_l',
            required=True,
            help = 'Ligand coordinate file(s): .mol2')

        parser.add_argument('-r',
            type=str,
            dest='input_file_r',
            required=True,
            help = 'Receptor coordinate file(s): .pdb')

        parser.add_argument('-f',
            dest='config_file',
            required=True,
            help='config file containing docking parameters')

        return parser

    def run_scoring(self):
        """Run scoring on original poses provided"""

        parser = self.create_arg_parser()
        args = parser.parse_args()

        print "Setting up parameters..."
        config = DockingConfig(args, task='scoring')

        tcpu1 = time.time()
        file_r = config.input_file_r
        config_s = config.scoring

        print "Starting scoring..."
        for kdx in range(len(config_s.site)):
            site = config_s.site['site'+str(kdx+1)]

            # iterate over rescoring instances
            for instance, program, options in config_s.instances:

                # get docking class
                ScoringClass = getattr(sys.modules[program], program.capitalize())
 
                ScoringInstance = ScoringClass(instance, site, options)
                outputfile = ScoringInstance.run_rescoring(config.input_file_r, [config.input_file_l], cleanup=config_s.cleanup)
 
        tcpu2 = time.time()
        print "Scoring done. Total time needed: %i s" %(tcpu2-tcpu1)

    def run_rescoring(self, config, args):
        """Run rescoring on docking poses"""

        tcpu1 = time.time()

        file_r = config.input_file_r
        config_r = config.rescoring
        posedir = args.posedir

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
        print "Starting rescoring..."
        # iterate over rescoring instances
        for instance, program, options in config_r.instances:

            # possibility of renaming the folder and output file 
            if 'name' in options:
                name = options['name']
            else:
                name = instance

            # remove old scoring file
            if os.path.isfile(name+'.score'):
                os.remove(name+'.score')

            for kdx in range(len(config_r.site)):
                site = config_r.site['site'+str(kdx+1)]

                # get complex filenames
                files_l = [os.path.abspath('../'+posedir+'/lig-%s.mol2'%idx) for idx in range(nposes[kdx], nposes[kdx+1])]
                # get docking class
                ScoringClass = getattr(sys.modules[program], program.capitalize())

                ScoringInstance = ScoringClass(instance, site, options)
                outputfile = ScoringInstance.run_rescoring(file_r, files_l, cleanup=config_r.cleanup)

                # cat output in file (cat instead of copying because of the binding sites)
                subprocess.check_output('cat %s >> %s'%(outputfile,name+'.score'), shell=True, executable='/bin/bash')

        os.chdir(curdir)
        tcpu2 = time.time()
        print "Rescoring done. Total time needed: %i s" %(tcpu2-tcpu1)

class Docking(object):

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="""rundock : dock with multiple software --------
Requires one file for the ligand (1 struct.) and one file for the receptor (1 struct.)""")

        parser.add_argument('-l',
            type=str,
            dest='input_file_l',
            required=True,
            help = 'Ligand coordinate file(s): .mol2')

        parser.add_argument('-r',
            type=str,
            dest='input_file_r',
            required=True,
            help = 'Receptor coordinate file(s): .pdb')

        parser.add_argument('-f',
            dest='config_file',
            required=True,
            help='config file containing docking parameters')

        parser.add_argument('-d',
            dest='posedir',
            default='poses',
            help='Directory containing poses to rescore (should be used with rescore_only option)')

        parser.add_argument('-prepare_only',
            dest='prepare_only',
            action='store_true',
            help='Only prepare scripts for docking (do not run docking)')

        parser.add_argument('-rescore_only',
            dest='rescore_only',
            action='store_true',
            default=False,
            help='Run rescoring only')

        parser.add_argument('-skip_docking',
            dest='skip_docking',
            action='store_true',
            default=False,
            help='Skip docking (used for debugging minimization step prior to rescoring)')

        return parser

    def finalize(self, config):
        """create directory containing all the poses found!"""

        config_d = config.docking

        resultdir = 'poses'
        shutil.rmtree(resultdir, ignore_errors=True)
        os.mkdir(resultdir)

        nposes = [1] # number of poses involved for each binding site
        sh = 1 # shift of model

        info = {}
        features = ['program', 'nposes', 'firstidx', 'site']
        for ft in features:
            info[ft] = []

        for kdx in range(len(config_d.site)):
            bs = config_d.site['site'+str(kdx+1)] # current binding site
            for name, program, options in config_d.instances:
                # find name for docking directory
                instdir = '%s'%name
                if bs[0]:
                    instdir += '.' + bs[0]                
                poses_idxs = []
                for filename in glob(instdir+'/lig-*.mol2'):
                    poses_idxs.append(int((filename.split('.')[-2]).split('-')[-1]))
                poses_idxs = sorted(poses_idxs)
                nposes_idxs = len(poses_idxs)

                for idx, pose_idx in enumerate(poses_idxs):
                    shutil.copyfile(instdir+'/lig-%s.mol2'%pose_idx, resultdir+'/lig-%s.mol2'%(idx+sh))

                # update info
                info['program'].append(name)
                info['nposes'].append(nposes_idxs)
                info['firstidx'].append(sh)
                info['site'].append(bs[0])

                # update shift
                sh += nposes_idxs
            nposes.append(sh)

        # write info
        info = pd.DataFrame(info)
        info[features].to_csv(resultdir+'/info.dat', index=False)

        # insert line at the beginning of the info file
        with open(resultdir+'/info.dat', 'r+') as ff:
            content = ff.read()
            ff.seek(0, 0)
            line = '#' + ','.join(map(str,nposes))+'\n'
            ff.write(line.rstrip('\r\n') + '\n' + content)

        # copy receptor in folder
        shutil.copyfile(config.input_file_r, resultdir+'/rec.pdb')

    def run_docking(self, config, args):
        """Running docking simulations using each program specified..."""

        tcpu1 = time.time()

        config_d = config.docking
        # iterate over all the binding sites
        for kdx in range(len(config.docking.site)):
            for instance, program, options in config.docking.instances: # iterate over all the instances

                # get docking class
                DockingClass = getattr(sys.modules[program], program.capitalize())

                # create docking instance and run docking
                DockingInstance = DockingClass(instance, config.docking.site['site'+str(kdx+1)], options)
                DockingInstance.run_docking(config.input_file_r, config.input_file_l, minimize_options=config_d.minimize, \
cleanup=config_d.cleanup, cutoff_clustering=config_d.cutoff_clustering, prepare_only=args.prepare_only, skip_docking=args.skip_docking)

        if args.prepare_only:
            return

        tcpu2 = time.time()
        print "Docking procedure done. Total time needed: %i s" %(tcpu2-tcpu1)

    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args()    

        print "Setting up parameters..."
        config = DockingConfig(args)

        # run docking
        if not args.rescore_only:
            self.run_docking(config, args)

        if args.prepare_only:
            return

        if not args.rescore_only:
            # create folder with poses
            self.finalize(config)

        # run rescoring
        if config.rescoring.is_rescoring:
            Scoring().run_rescoring(config, args)
