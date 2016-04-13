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
import ligprep
import glob

import numpy as np
import moe

class PrepDocking(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Prepare files for Docking or Virtual Screening")

        parser.add_argument('-l',
            type=str,
            dest='input_file_l',
            required=True,
            nargs='+',
            help = 'Ligand coordinate file(s): .sdf, .smi')

        parser.add_argument('-r',
            type=str,
            dest='input_file_r',
            required=True,
            nargs='+',
            help = 'Receptor coordinate file(s): .pdb')

        parser.add_argument('-f',
            type=str,
            dest='config_file',
            help='config file: .ini')

        parser.add_argument('--findsites',
            dest='findsites',
            action='store_true',
            default=False,
            help='Find possible binding sites (use MOE\'s site finder)')

        parser.add_argument('--minplb',
            type=int,
            dest='minplb',
            default=1.0,
            help='Minimum PLB value (MOE) to select the binding sites (requires findsites option). Default: 1.0 ')

        parser.add_argument('--nsitesmax',
            type=int,
            dest='nsitesmax',
            default=0,
            help='Maximum number of binding sites kept (when 0 all the binding sites are retained, requires findsites option). Default: 0')

        parser.add_argument('--inplace',
            dest='inplace',
            action='store_true',
            default=False,
            help='Do not create new directories when single input files are specified')

        parser.add_argument('--lpflags',
            type=str,
            default="-W e,-ph,7.0,-pht,2.0 -epik -r 1 -bff 14",
            dest='ligprep_flags',
            help='Ligprep (Schrodinger) flags for ligand preparation. Default: "-W e,-ph,7.0,-pht,2.0 -epik -r 1 -bff 14"')

        parser.add_argument('--nolp',
            dest='no_ligprep',
            action='store_true',
            default=False,
            help='No ligand preparation with ligprep')

        parser.add_argument('--pwzdflags',
            type=str,
            default="-fix -fillsidechains -pH \'neutral\'",
            dest='prepwizard_flags',
            help='Prepwizard (Schrodinger) flags for protein preparation. Default: "-fix -fillsidechains -pH \'neutral\'"')

        parser.add_argument('--nopwzd',
            dest='no_prepwizard',
            action='store_true',
            default=False,
            help='No protein preparation with prepwizard')

        parser.add_argument('--skip',
            dest='skip',
            action='store_true',
            default=False,
            help='Skip protein and ligand preparation (used for debugging)')

        return parser

    def initialize(self, args):
        # check if config file exist
        if args.config_file:
            if not os.path.exists(args.config_file):
                raise ValueError("Config file %s not found!"%(args.config_file))
        else:
            print "Warning: no config file provided! Please provide config file if you want it to be updated!"

        # get full names of ligand input files
        self.input_file_l = []
        if args.input_file_l:
            for file_l in args.input_file_l:
                # check if ligand file exists
                if not os.path.exists(file_l):
                    raise ValueError("File %s not found!"%(self.input_file_l))
                self.input_file_l.append(os.path.abspath(file_l))
        self.nfiles_l = len(self.input_file_l)
            
        # get full names of receptor input files
        self.input_file_r = []
        if args.input_file_r:
            for file_r in args.input_file_r:
                # check if ligand file exists
                if not os.path.exists(file_r):
                    raise ValueError("File %s not found!"%(self.input_file_r))
                self.input_file_r.append(os.path.abspath(file_r))
        self.nfiles_r = len(self.input_file_r)

    def prepare_vs(self, args):

        if args.skip:
            # get ligand file names
            self.files_prep_l = []
            for idx in range(self.nfiles_l):
                output_mol2files = []
                ligpdir = 'lig-prep' + str(idx+1)
                suffix, ext = os.path.splitext(glob.glob(ligpdir+'/*.mol2')[0])
                for idx in range(len(glob.glob(suffix[:-1] + '*.mol2'))):
                    output_mol2files.append(os.path.abspath(suffix[:-1] + '%s.mol2'%(idx+1))) 
                self.files_prep_l.append(output_mol2files)
            # get receptor file names
            self.files_prep_r = []
            for idx in range(self.nfiles_r):
                recpdir = 'rec-prep' + str(idx+1)
                self.files_prep_r.append(glob.glob(recpdir+'/*.pdb')[0])
        else:
            # remove existing directories
            for ligdir in glob.glob('lig-prep*'):
                shutil.rmtree(ligdir)
            for recdir in glob.glob('rec-prep*'):
                shutil.rmtree(recdir)
            # prepare structures
            self.prepare_structures(args)

        # (A) LIGAND level
        for idx in range(self.nfiles_l):
            if self.nfiles_l == 1 and not args.inplace:
                ligdir = '.'
            else:
                ligdir = 'lig' + str(idx+1)
                if os.path.isdir(ligdir):
                    shutil.rmtree(ligdir)
                os.mkdir(ligdir)

            # (B) RECEPTOR level
            for jdx in range(self.nfiles_r):
                recdir = ligdir
                if self.nfiles_r == 1 and not args.inplace:
                    recdir += '/.'
                else:
                    recdir += '/rec' + str(jdx+1)
                    if os.path.isdir(recdir):
                        shutil.rmtree(recdir)
                    os.mkdir(recdir)

                # (C) ISOMER level
                file_r = self.files_prep_r[jdx]
                recpdir = 'rec-prep%s'%(jdx+1)
                for kdx, file_l in enumerate(self.files_prep_l[idx]):
                    suffix = os.path.splitext(file_l)[0]
                    isodir = recdir + '/isomer' + str(kdx+1)
                    shutil.rmtree(isodir, ignore_errors=True)
                    os.mkdir(isodir)
                    if args.config_file:
                        if args.findsites:
                            table = np.loadtxt(recpdir+'/sitefinder.log')
                            if len(table.shape) == 1:
                                table = table[np.newaxis,:]
                            # TODO: include the size of the ligand to compute boxsize
                            new_config_file = isodir + '/config.ini'
                            self.update_site_config_file(new_config_file, args.config_file, table)
                        else:
                            # copy the config file in the directory
                            shutil.copyfile(args.config_file, isodir + '/config.ini')
                    # copy the files for ligand and receptor in the corresponding dir
                    shutil.copyfile(file_l, isodir + '/lig.mol2')
                    shutil.copyfile(file_r, isodir +'/rec.pdb')

    def generate_mol2_files(self, file_l, args):

        suffix, ext = os.path.splitext(file_l)
        if ext == '.sdf':
            input_format_flag = '-isdf'
        elif ext in ['.smi', '.txt']:
            input_format_flag = '-ismi'
        else:
            raise IOError("Format %s not recognized!"%(ext[1:]))

        # generate multiple .mol2 files using babel
        output_mol2file_model = suffix + '_.mol2'
        subprocess.check_call('babel %s %s -omol2 %s -m 2>/dev/null'%(input_format_flag,file_l,output_mol2file_model),shell=True)

        output_mol2files = []
        for idx in range(len(glob.glob('*.mol2'))):
            output_mol2files.append(os.path.abspath(suffix + '_%s.mol2'%(idx+1)))
        return output_mol2files

    def prepare_structures(self, args):

        curdir = os.getcwd()
        print "Preparing ligands..."
        self.files_prep_l = []
        for idx, file_l in enumerate(self.input_file_l):
            # create new directory
            ligpdir = 'lig-prep' + str(idx+1)
            shutil.rmtree(ligpdir, ignore_errors=True)
            os.mkdir(ligpdir)
            os.chdir(ligpdir) # change directory
            # (A) Run Schrodinger's ligprep
            if not args.no_ligprep:
                new_file_l = ligprep.prepare_ligand(file_l, args.ligprep_flags)
            else:
                new_file_l = file_l
            # (B) Generate mol2files using babel
            mol2files = self.generate_mol2_files(new_file_l, args)
            self.files_prep_l.append(mol2files)
            os.chdir(curdir)

        print "Preparing receptors..."
        self.files_prep_r = []
        for idx, file_r in enumerate(self.input_file_r):
            recpdir = 'rec-prep' + str(idx+1)
            shutil.rmtree(recpdir, ignore_errors=True)
            os.mkdir(recpdir)
            os.chdir(recpdir)
            # (A) Run Schrodinger's Prepwizard
            if not args.no_prepwizard:
                new_file_r = ligprep.prepare_receptor(file_r, args.prepwizard_flags)
            else:
                new_file_r = os.path.basename(file_r)
                shutil.copyfile(file_r, new_file_r)
            # (B) Run MOE's site finder
            if args.findsites:
                self.find_binding_sites(new_file_r, args)
            self.files_prep_r.append(new_file_r)
            os.chdir(curdir)

    def update_site_config_file(self, new_config_file, config_file, table):

        shutil.copyfile(config_file, new_config_file)

        # create tmp config file name from original config file
        tmp_config_file = list(os.path.splitext(new_config_file))
        tmp_config_file.insert(1,'_tmp')
        tmp_config_file = ''.join(tmp_config_file)

        # remove section 'SITE' and option site in DOCKING section of config file if exists
        with open(tmp_config_file, 'w') as tmpf:
            with open(new_config_file, 'r') as newf:
                isdock = False
                sitesection = False
                docksection = False
                for line in newf:
                    # check if still in section SITE*
                    if line.startswith(('[SITE','[LIGPREP]')):
                        sitesection = True
                    if sitesection and line.startswith('[') and not line.startswith(('[SITE','[LIGPREP]')): # new section has been reached
                        sitesection = False
                    # check if still in section DOCKING
                    if line.startswith('[DOCKING]'):
                        docksection = True
                        isdock = True
                    if docksection and line.startswith('[') and not line.startswith('[DOCKING]'): # new section has been reached
                        docksection = False
                     # check if option line in section DOCKING
                    if line.strip().startswith('site') and docksection:
                        siteline = True
                    else:
                        siteline = False
                    if not sitesection and not siteline:
                        tmpf.write(line)
        shutil.move(tmp_config_file, new_config_file)

        # add new sections 'SITE' and option site
        with open(tmp_config_file, 'w') as tmpf:
            with open(new_config_file, 'r') as newf:
                for line in newf:
                    tmpf.write(line)
                    if line.startswith('[DOCKING]'):
                        tmpf.write('site = ' + ', '.join(['site%s'%int(line[0]) for line in table])+'\n')
                for line in table:
                    section = 'SITE' + str(int(line[0]))
                    center_conf = ', '.join(map(str, line[2:5].tolist()))
                    boxsize_conf = ', '.join(map(str, [2*line[5] for idx in range(3)]))
                    newsite_section = """
[%(section)s]
center = %(center_conf)s
boxsize = %(boxsize_conf)s"""% locals()

                    tmpf.write(newsite_section+'\n')

        shutil.move(tmp_config_file, new_config_file)

    def find_binding_sites(self, pdbfile, args):

        # (A) write script
        script_name = 'find_binding_sites.sh'
        moe.write_sitefinder_script(script_name, pdbfile, args)
        os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

        # (B) execute script
        subprocess.check_call("./" + script_name + " &> sitefinder.log", shell=True, executable='/bin/bash')

    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args()

        tcpu1 = time.time()
        self.initialize(args)

        # prepare virtual screening
        self.prepare_vs(args)
        tcpu2 = time.time()

if __name__ == '__main__':
    PrepDocking().run()
