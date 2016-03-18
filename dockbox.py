#!/global/software/python-2.7.2/bin/python

from __future__ import with_statement

import sys
import os
import subprocess
import tempfile
import shutil
import argparse
import ConfigParser
import logging
import stat
import time
import prep
import consensus
import numpy as np
import prep

import util.pdbtools as pdbt

known_programs = ['autodock', 'vina', 'dock', 'glide', 'moe']
known_extract_options = ['all', 'lowest', 'none']

class DockingConfigError(Exception):
    pass

class DockingConfig(object):

    def __init__(self, args):

        # check if config file exist
        if not os.path.exists(args.config_file):
            raise DockingConfigError("Config file %s not found!"%(args.config_file))

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file)

        if config.has_option('DOCKING', 'program'):
            self.programs = []
            programs = config.get('DOCKING', 'program').lower()
            programs = map(str.strip, programs.split(','))
            for program in programs:
                if program not in known_programs:
                    raise DockingConfigError("program option should be one of " + ", ".join(known_programs))
                self.programs.append(program)
        else:
            self.programs = ['autodock']

        for program in self.programs:
            __import__(program)

        if config.has_option('GENERAL', 'system'):
            self.system = config.get('GENERAL', 'system').lower()
        else:
            self.system = None

        options = {}

        for program in self.programs:
            # check if all needed executables are available
            if hasattr(sys.modules[program], 'default_settings'):
                required_programs = getattr(sys.modules[program], 'required_programs')
                for exe in required_programs:
                    try:
                        subprocess.check_call('which %s > /dev/null'%exe, shell=True)
                    except subprocess.CalledProcessError:
                        raise DockingConfigError('Executable %s needed for docking with %s is not found in your PATH! \
Make sure the program has been installed!'%(exe,program))

            options[program] = {} 
            # load default parameters
            if hasattr(sys.modules[program], 'default_settings'):
                default_settings = getattr(sys.modules[program], 'default_settings')
                for key, value in default_settings.iteritems():
                    options[program][key] = value

            # load preset parameters
            if hasattr(sys.modules[program], 'known_settings'):
                known_settings = getattr(sys.modules[program], 'known_settings')
                if self.system in known_settings:
                    # print "System known with %s: "%program + self.system
                    for key, value in known_settings[self.system].iteritems():
                        options[program][key] = value

            # check config file (would possibly overwrite preset parameters)
            if config.has_section(program.upper()):
               config_d = dict(config.items(program.upper()))
               for key, value in config_d.iteritems():
                   options[program][key] = value

            # check if all required options have been set
            if hasattr(sys.modules[program], 'required_settings_names'):
                required_settings_names = getattr(sys.modules[program], 'required_settings_names')
                for option in required_settings_names:
                    if option not in options[program]:
                        raise DockingConfigError("option %s required to use %s. Check manual!" %(option,program.capitalize()))

        # save options
        self.options = options

        self.input_file_l = os.path.abspath(args.input_file_l[0])
        # check if ligand file exists
        if not os.path.exists(self.input_file_l):
            raise DockingConfigError("File %s not found!"%(self.input_file_l))

        self.input_file_r = os.path.abspath(args.input_file_r[0])
        # check if receptor file exists
        if not os.path.exists(self.input_file_r):
            raise DockingConfigError("File %s not found!"%(self.input_file_r))

        if config.has_section('LIGPREP'):
            ligprep = dict(config.items('LIGPREP'))
        else:
            ligprep = {}

        # set defaults options if they are not mentionned in the config file
        for key, value in prep.ligprep_default_options.items():
            if key not in ligprep:
                ligprep[key] = value
        self.ligprep = ligprep

        # get value of option cleanup
        if config.has_option('DOCKING', 'cleanup'):
            yesno = config.get('DOCKING', 'cleanup').lower()
            if yesno == 'yes':
                cleanup = True
            elif yesno == 'no':
                cleanup = False
            else:
                raise DockingConfigError("option consensus should be yes or no!")
        else:
            cleanup = False
        self.cleanup = cleanup

        # get value of option extract
        if config.has_option('DOCKING', 'extract'):
            extract = config.get('DOCKING', 'extract').lower()
            if not extract in known_extract_options:
                raise DockingConfigError("Extract option should be one of " + ", ".join(known_extract_options))
        else:
            extract = 'lowest'
        self.extract = extract

        # set options for consensus docking
        self.initialize_consensus_options(config)

    def initialize_consensus_options(self, config):

        options = {}
        if config.has_option('DOCKING', 'consensus'):
            yesno = config.get('DOCKING', 'consensus').lower()
            if yesno == 'yes':
                import consensus
                options['consensus'] = True
            elif yesno == 'no':
                options['consensus'] = False
            else:
                raise DockingConfigError("option consensus should be yes or no!")
        else:
            options['consensus'] = False

        if options['consensus']:
            # load default parameters
            if hasattr(sys.modules['consensus'], 'default_settings'):
                default_settings = getattr(sys.modules['consensus'], 'default_settings')
                for key, value in default_settings.iteritems():
                    options[key] = value
            # check config file (would possibly overwrite preset parameters)
            if config.has_section('CONSENSUS'):
                config_c = dict(config.items('CONSENSUS'))
                for key, value in config_c.iteritems():
                    options[key] = value
  
        self.consensus = options

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

        return parser

    def prepare_structures(self, config):

        curdir = os.getcwd()
        prepdir = 'prep'
        if os.path.isdir(prepdir):
            shutil.rmtree(prepdir)
        os.mkdir(prepdir)
        os.chdir(prepdir)

        print "Preparing ligand..."
        prep.prepare_ligand(config)

        # TODO: add possible receptor preparation using prepwizard
        # copy the file containing the receptor
        print "Preparing receptor..."
        shutil.copyfile(config.input_file_r, 'rec.pdb')
        config.input_file_r = os.path.abspath('rec.pdb')
        
        os.chdir(curdir)

    def run_docking(self, config):
        """Running docking simulations using each program specified..."""

        curdir = os.getcwd()

        for program in config.programs:
            for idx, file_l in enumerate(config.input_files_l):
                print "Starting docking with %s..."%program.capitalize()
                print "The following options will be used:"
                print config.options[program]

                if config.nisomers == 1:
                    prgdir = program
                else:
                    prgdir = program + '.' + str(idx+1)

                if os.path.isdir(prgdir):
                    shutil.rmtree(prgdir)
                os.mkdir(prgdir)
                os.chdir(prgdir)

                # (A) write docking script
                script_name = "run_" + program + ".sh"
                write_docking_script = getattr(sys.modules[program], 'write_docking_script')
                write_docking_script(script_name, config.input_file_r, file_l, config)
                os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

                # (B) execute docking procedure
                subprocess.check_call("./" + script_name + " &> " + program + ".log", shell=True, executable='/bin/bash')

                # (C) extract docking results
                extract_docking_results = getattr(sys.modules[program], 'extract_docking_results')
                extract_docking_results('rec.out.pdb', 'lig.out.pdb', 'score.out', config)

                # (D) remove intermediate files if required
                if config.cleanup:
                    cleanup = getattr(sys.modules[program], 'cleanup')
                    cleanup(config)
 
                os.chdir(curdir)
                print "Docking with %s done."%program.capitalize()

        if config.consensus['consensus']:
            consensus.find_consensus(config)

    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args()    

        tcpu1 = time.time()

        print "Setting up parameters..."
        config = DockingConfig(args)

        # prepare ligand (ligprep, 2D -> 3D, provide unique atom names)
        self.prepare_structures(config)

        # run docking
        self.run_docking(config)

        tcpu2 = time.time()
        print "Docking procedure done. Total time needed: %i s" %(tcpu2-tcpu1)

if __name__ == '__main__':
    DockingWorker().run()
