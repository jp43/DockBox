#!/global/software/python-2.7.2/bin/python
from __future__ import with_statement

import sys
import os
import subprocess
import shutil
import argparse
import ConfigParser
import stat
import time

import multi
import consensus
import numpy as np

class DockingConfig(object):

    def __init__(self, args):

        # check if config file exist
        if not os.path.exists(args.config_file):
            raise ValueError("Config file %s not found!"%(args.config_file))

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file)

        self.input_file_l = os.path.abspath(args.input_file_l[0])
        # check if ligand file exists
        if not os.path.exists(self.input_file_l):
            raise IOError("File %s not found!"%(self.input_file_l))

        self.input_file_r = os.path.abspath(args.input_file_r[0])
        # check if receptor file exists
        if not os.path.exists(self.input_file_r):
            raise IOError("File %s not found!"%(self.input_file_r))

        self.docking = multi.MultiProgramDocking(config)
        self.consensus = consensus.ConsensusDocking(config, args)

        self.cleanup = self.is_yesno_option(config, 'cleanup')
        self.extract_only = args.extract_only

    def is_yesno_option(self, config, option, default=False):

        if config.has_option('DOCKING', option):
            yesno = config.get('DOCKING', option).lower()
            if yesno == 'yes':
                return True
            elif yesno == 'no':
                return False
            else:
                raise ValueError("option %s should be yes or no!"%option)
        else:
            return default

class DockingWorker(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Run Docking")

        parser.add_argument('-l',
            type=str,
            dest='input_file_l',
            required=True,
            nargs=1,
            help = 'Ligand coordinate file(s): .pdb, .sdf, .smi')

        parser.add_argument('-r',
            type=str,
            dest='input_file_r',
            required=True,
            nargs=1,
            help = 'Receptor coordinate file(s): .pdb')

        parser.add_argument('-f',
            dest='config_file',
            required=True,
            help='config file containing docking parameters')

        parser.add_argument('--consensus_only',
            dest='consensus_only',
            action='store_true',
            default=False,
            help='Run consensus only')

        parser.add_argument('--extract_only',
            dest='extract_only',
            action='store_true',
            default=False,
            help='Extract structures only (used for debugging)')

        return parser


    def run_docking(self, config):
        """Running docking simulations using each program specified..."""

        tcpu1 = time.time()

        # iterate over all the binding sites
        for kdx in range(len(config.docking.site)):

            for instance, program, options in config.docking.instances: # iterate over all the instances

                # get docking class
                DockingClass = getattr(sys.modules[program], program.capitalize())

                # create docking instance and run docking
                DockingInstance = DockingClass(instance, config.docking.site['site'+str(kdx+1)], options)
                DockingInstance.run_docking(config.input_file_r, config.input_file_l, config.docking.extract, \
config.cleanup, extract_only=config.extract_only)

        tcpu2 = time.time()
        print "Docking procedure done. Total time needed: %i s" %(tcpu2-tcpu1)

    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args()    

        print "Setting up parameters..."
        config = DockingConfig(args)

        # run docking
        if not config.consensus.only:
            self.run_docking(config)

        config.consensus.find_consensus(config.docking.instances, config.input_file_r, config.docking.site)

