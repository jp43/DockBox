#!/usr/bin/python
from __future__ import with_statement

import os
import sys
import shutil
import argparse
import ConfigParser
import time

from glob import glob
import pandas as pd
import subprocess

from mdkit.utility import mol2
from mdkit.amber.ambertools import load_PROTON_INFO
from mdkit.amber.ambertools import load_atomic_ions
import setconf

class DockingConfig(object):

    def __init__(self, args, task='docking'):

        # check if config file exist
        if os.path.exists(args.config_file):
            config = ConfigParser.SafeConfigParser()
            config.read(args.config_file)
        else:
            raise ValueError("Config file %s not found!"%(args.config_file))

        # check if ligand file exists
        if not os.path.isfile(args.input_file_l):
            raise IOError("File %s not found!"%(args.input_file_l))

        file_l_abs = os.path.abspath(args.input_file_l)
        pref, ext = os.path.splitext(file_l_abs)
        if ext != '.mol2':
            raise IOError("Ligand file provided with -l option should be in .mol2 format! %s format detected!"%ext)

        nligands = int(subprocess.check_output('fgrep -c "@<TRIPOS>ATOM" %s'%file_l_abs, shell=True))
        if nligands == 0:
            raise IOError("No ligand detected in %s, check your file again!"%args.input_file_l)
        elif nligands > 1:
            raise IOError("More than one ligand detected in %s. Only one structure per ligand file is allowed!"%args.input_file_l)

        # new ligand file with unique names for every atom
        new_file_l = pref + '_dbx' + ext

        # create a ligand file with unique atom names
        mol2.update_mol2file(file_l_abs, new_file_l, unique=True, ligname='LIG')
        self.input_file_l = os.path.abspath(new_file_l)

        if task == 'docking':
            self.docking = setconf.DockingSetup(config)
            self.rescoring = setconf.RescoringSetup(config)
        elif task == 'scoring':
            self.scoring = setconf.ScoringSetup(config)
        else:
            raise ValueError("Task should be one of docking or scoring")

        self.check_pdbfile(args.input_file_r)

    def check_pdbfile(self, filename):

        # check if receptor file exists
        if not os.path.isfile(filename):
            raise IOError("File %s not found!"%(filename))

        proton_info = load_PROTON_INFO()
        ions_info = load_atomic_ions()

        # check if the .pdb file is valid
        with open(filename, 'r') as pdbf:
            is_end_line = False
            for line in pdbf:
                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    if resname in ions_info:
                        for instance, program, options in self.docking.instances:
                            if program not in setconf.programs_handling_ions:
                                sys.exit("Ion %s found in structure %s! DockBox is not configured to apply %s with ions!" %(resname, filename, program))
                    elif resname not in proton_info or line.startswith('HETATM'):
                        sys.exit('Unrecognized residue %s found in %s! The .pdb file should \
only contains one protein structure with standard residues (with possibly ions)!'%(resname, filename))
                    elif is_end_line:
                        sys.exit("More than one structure detected in .pdb file! Check your file again!")
                elif line.startswith('END'):
                    is_end_line = True

        self.input_file_r = os.path.abspath(filename)

class Scoring(object):

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="""runscore : score in-place with multiple software --------
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
                outputfile = ScoringInstance.run_rescoring(config.input_file_r, [config.input_file_l])
 
        tcpu2 = time.time()
        print "Scoring done. Total time needed: %i s" %(tcpu2-tcpu1)

    def run_rescoring(self, config, args):
        """Run rescoring on docking poses"""

        tcpu1 = time.time()

        file_r = config.input_file_r
        config_r = config.rescoring
        posedir = 'poses'

        # look for results folder
        if os.path.isdir(posedir):
            with open(posedir+'/info.dat') as inff:
                nposes = inff.next()
                nposes = nposes[1:] # the first character is a # sign
                nposes = map(int, nposes.split(','))
        else:
            raise IOError('no folder %s found!'%posedir)

        curdir = os.getcwd()
        workdir = 'rescoring'
        if not os.path.exists(workdir):
            print "Creating rescoring folder..."
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
                outputfile = ScoringInstance.run_rescoring(file_r, files_l)

                # cat output in file (cat instead of copying because of the binding sites)
                subprocess.check_output('cat %s >> %s'%(outputfile,name+'.score'), shell=True, executable='/bin/bash')

                if config.docking.cleanup >= 1:
                    shutil.rmtree(os.path.dirname(outputfile), ignore_errors=True)

        os.chdir(curdir)
        tcpu2 = time.time()
        print "Rescoring done. Total time needed: %i s" %(tcpu2-tcpu1)

class Docking(object):

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="""rundbx : dock and rescore with multiple programs --------
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

        parser.add_argument('-prepare_only',
            dest='prepare_only',
            action='store_true',
            help='Only prepare scripts for docking (does not run docking)')

        parser.add_argument('-rescore_only',
            dest='rescore_only',
            action='store_true',
            default=False,
            help='Run rescoring only')

        parser.add_argument('-skip_docking',
            dest='skip_docking',
            action='store_true',
            default=False,
            help=argparse.SUPPRESS)

        return parser

    def finalize(self, config):
        """create directory containing all the poses found!"""

        config_d = config.docking

        posedir = 'poses'
        shutil.rmtree(posedir, ignore_errors=True)
        os.mkdir(posedir)

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
                    shutil.copyfile(instdir+'/lig-%s.mol2'%pose_idx, posedir+'/lig-%s.mol2'%(idx+sh))

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
        info[features].to_csv(posedir+'/info.dat', index=False)

        # insert line at the beginning of the info file
        with open(posedir+'/info.dat', 'r+') as ff:
            content = ff.read()
            ff.seek(0, 0)
            line = '#' + ','.join(map(str,nposes))+'\n'
            ff.write(line.rstrip('\r\n') + '\n' + content)

        # copy receptor in folder
        shutil.copyfile(config.input_file_r, posedir+'/rec.pdb')

    def do_final_cleanup(self, config):

        shutil.rmtree('poses', ignore_errors=True)
        for dir in glob('*'):
            scorefile = dir+'/score.out'

            if os.path.isfile(scorefile):
                for filename in glob(dir+'/*'):
                    if scorefile != filename:
                        if os.path.isfile(filename):
                            os.remove(filename)
                        elif os.path.isdir(filename):
                            shutil.rmtree(filename, ignore_errors=True)
                        else:
                            raise ValueError('Folder or file not recognized for %s'%filename)
        os.remove(config.input_file_l)


    def run_docking(self, config, args):
        """Running docking simulations using each program specified..."""

        tcpu1 = time.time()

        config_d = config.docking
        # iterate over all the binding sites
        for kdx in range(len(config_d.site)):
            for instance, program, options in config_d.instances: # iterate over all the instances

                # get docking class
                DockingClass = getattr(sys.modules[program], program.capitalize())

                # create docking instance and run docking
                DockingInstance = DockingClass(instance, config_d.site['site'+str(kdx+1)], options)
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

        # final cleanup if needed
        if config.docking.cleanup == 2:
            self.do_final_cleanup(config)
