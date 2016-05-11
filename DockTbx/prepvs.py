#!/usr/bin/python
from __future__ import with_statement

import sys
import os
import subprocess
import shutil
import argparse
import ConfigParser
import stat
import time
import glob
import numpy as np

from DockTbx import moe
from DockTbx.prep import ligprep
from DockTbx.amber import minimz as mn

class PrepDocking(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Prepare files for Docking or Virtual Screening")

        parser.add_argument('-l',
            type=str,
            dest='input_file_l',
            nargs='+',
            help = 'Ligand coordinate file(s): .sdf, .smi')

        parser.add_argument('-r',
            type=str,
            dest='input_file_r',
            nargs='+',
            help = 'Receptor coordinate file(s): .pdb')

        parser.add_argument('-f',
            type=str,
            dest='config_file',
            help='config file: .ini')

        parser.add_argument('-findsites',
            dest='findsites',
            action='store_true',
            default=False,
            help='Find possible binding sites (use MOE\'s site finder)')

        parser.add_argument('-minplb',
            type=float,
            dest='minplb',
            default=1.0,
            help='Minimum PLB value (MOE) to select the binding sites (requires findsites option). Default: 1.0 ')

        parser.add_argument('-nsitesmax',
            type=int,
            dest='nsitesmax',
            default=0,
            help='Maximum number of binding sites kept (when 0 all the binding sites are retained, requires findsites option). Default: 0')

        parser.add_argument('-inplace',
            dest='inplace',
            action='store_true',
            default=False,
            help='Do not create new directories when single input files are specified')

        parser.add_argument('-ligpflags',
            type=str,
            default="-W e,-ph,7.0,-pht,2.0 -epik -r 1 -bff 14",
            dest='ligprep_flags',
            help='Ligprep (Schrodinger) flags for ligand preparation. Default: "-W e,-ph,7.0,-pht,2.0 -epik -r 1 -bff 14"')

        parser.add_argument('-noligp',
            dest='no_ligprep',
            action='store_true',
            default=False,
            help='No ligand preparation with ligprep')

        parser.add_argument('-pwzdflags',
            type=str,
            default="-fix -fillsidechains -pH \'neutral\'",
            dest='prepwizard_flags',
            help='Prepwizard (Schrodinger) flags for protein preparation. Default: "-fix -fillsidechains -pH \'neutral\'"')

        parser.add_argument('-nopwzd',
            dest='no_prepwizard',
            action='store_true',
            default=False,
            help='No protein preparation with prepwizard')

        parser.add_argument('-update',
            dest='update',
            action='store_true',
            default=False,
            help='Update directories and files only (used for debugging)')

        return parser

    def initialize(self, args):

        # check if config file exist
        if args.config_file:
            if not os.path.exists(args.config_file):
                raise ValueError("Config file %s not found!"%(args.config_file))
            else:
                self.config_file = args.config_file
        else:
            self.config_file = None
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

    def cleanup(self):
        """Cleanup directories used for preparation"""
        # remove existing directories

        for ligdir in glob.glob('lig-prep*'):
            shutil.rmtree(ligdir)

        for recdir in glob.glob('rec-prep*'):
            shutil.rmtree(recdir)

    def make_dirs(self, root, type, ndirs, index, args):

        dirname = {'ligand': 'lig', 'receptor': 'rec', 'isomer': 'iso'}
        dir = root
        if ndirs == 1 and args.inplace:
            dir += '/.'
        else:
            dir += '/' + dirname[type] + str(index+1)
            shutil.rmtree(dir, ignore_errors=True)
            os.mkdir(dir)
        return dir

    def get_prepared_filenames_from_old_run(self):

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

    def prepare_vs(self, args):
        """Prepare files and directories for Virtual Screening"""

        if args.update:
            self.get_prepared_filenames_from_old_run()
        else:
            self.cleanup()
            # prepare structures
            self.prepare_structures(args)

        # ligand and receptor files provided
        if self.nfiles_l > 0 and self.nfiles_r > 0:
            # LIGAND level
            for idx in range(self.nfiles_l):
                ligdir = self.make_dirs('.', 'ligand', self.nfiles_l, idx, args)
                # RECEPTOR level
                for jdx in range(self.nfiles_r):
                    recdir = self.make_dirs(ligdir, 'receptor', self.nfiles_r, jdx, args)
                    # get prepared ligand filename
                    file_r = self.files_prep_r[jdx]
                    recpdir = 'rec-prep%s'%(jdx+1)
                    # ISOMER level
                    for kdx, file_l in enumerate(self.files_prep_l[idx]):
                        workdir = self.make_dirs(recdir, 'isomer', self.files_prep_l[idx], kdx, args)
                        if self.config_file:
                            self.update_config_file(workdir, config_file, recpdir, args)
                        # copy the files for ligand and receptor in the corresponding working dir
                        shutil.copyfile(file_l, workdir + '/lig.mol2')
                        shutil.copyfile(file_r, workdir +'/rec.pdb')
                        with open(workdir + '/vs.info', 'w') as ff:
                            print >> ff, 'Location of original ligand file: ' + self.input_file_l[idx]
                            print >> ff, 'Location of original receptor file: ' + self.input_file_r[jdx]

        # no ligand files provided
        elif self.nfiles_l == 0 and self.nfiles_r > 0:
            # RECEPTOR level
            for jdx in range(self.nfiles_r):
                workdir = self.make_dirs('.', 'receptor', self.nfiles_r, jdx, args)
                file_r = self.files_prep_r[jdx]
                recpdir = 'rec-prep%s'%(jdx+1)
                if self.config_file:
                    self.update_config_file(workdir, config_file, recpdir, args)
                shutil.copyfile(file_r, workdir +'/rec.pdb')
                with open(workdir + '/vs.info', 'w') as ff:
                    print >> ff, 'Location of original receptor file: ' + self.input_file_r[jdx]

        # no receptor files provided
        elif self.nfiles_l > 0 and self.nfiles_r == 0:
            # LIGAND level
            for idx in range(self.nfiles_l):
                ligdir = self.make_dirs('.', 'ligand', self.nfiles_l, idx, args)
                # ISOMER level
                for kdx, file_l in enumerate(self.files_prep_l[idx]):
                    workdir = self.make_dirs(ligdir, 'isomer', self.files_prep_l[idx], kdx, args)
                    if self.config_file:
                        self.update_config_file(workdir, config_file, recpdir, args)
                    shutil.copyfile(file_l, workdir + '/lig.mol2')
                    with open(workdir + '/vs.info', 'w') as ff:
                        print >> ff, 'Location of original ligand file: ' + self.input_file_l[idx]

        # no receptor and ligand files provided
        else:
            pass

    def generate_mol2file(self, file_l, args):

        suffix, ext = os.path.splitext(file_l)
        if ext == '.sdf':
            input_format_flag = '-isdf'
        elif ext in ['.smi', '.txt']:
            input_format_flag = '-ismi'
        else:
            raise IOError("Format %s not recognized!"%(ext[1:]))

        # generate multiple .mol2 files using babel
        output_mol2file_model = suffix + '_.mol2'
        subprocess.check_output('babel %s %s -omol2 %s -m 2>/dev/null'%(input_format_flag, file_l, output_mol2file_model), shell=True, executable='/bin/bash')

        output_mol2files = []
        for idx in range(len(glob.glob('*.mol2'))):
            output_mol2files.append(os.path.abspath(suffix + '_%s.mol2'%(idx+1)))

        return output_mol2files

    def prepare_structures(self, args):

        curdir = os.getcwd()

        if self.nfiles_l > 0:
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
                new_file_l = os.path.basename(file_l)
                shutil.copyfile(file_l, new_file_l)

            # (B) Generate mol2file using babel
            mol2files = self.generate_mol2file(new_file_l, args)
            self.files_prep_l.append(mol2files)
            os.chdir(curdir)

        if self.nfiles_r > 0:
            print "Preparing receptors..."
        self.files_prep_r = []

        for idx, file_r in enumerate(self.input_file_r):
            recpdir = 'rec-prep' + str(idx+1)
            shutil.rmtree(recpdir, ignore_errors=True)
            os.mkdir(recpdir)
            os.chdir(recpdir)

            # (A) clean-up PDBfile
            removed_residues = []
            new_file_r = os.path.basename(file_r)
            atoms_info = mn.load_PROTON_INFO()
            with open(file_r, 'r') as pdbi:
                with open(new_file_r, 'w') as pdbo:
                    for line in pdbi:
                        if line.startswith('ATOM'):
                            resname = line[17:20]
                            if resname in atoms_info:
                                pdbo.write(line)
                            elif resname not in removed_residues:
                                removed_residues.append(resname)
                                print "Removed unknown residue %s"%line[17:20]
                        elif line.startswith(('TER','END')):
                            pdbo.write(line)

            # (B) Run Schrodinger's Prepwizard
            if not args.no_prepwizard:
                new_file_r = ligprep.prepare_receptor(new_file_r, args.prepwizard_flags)

            # (C) Run MOE's site finder
            if args.findsites:
                self.find_binding_sites(new_file_r, args)

            self.files_prep_r.append(os.path.abspath(new_file_r))
            os.chdir(curdir)


    def update_config_file(self, workdir, config_file, recpdir, args):

        new_config_file = workdir+'/config.ini'

        if args.findsites:
            self.update_config_file_site_options(new_config_file, config_file, recpdir)

        elif not os.path.abspath(workdir) == os.getcwd():
            shutil.copyfile(config_file, new_config_file)


    def update_config_file_site_options(self, new_config_file, config_file, recpdir):
        """Update the config file provided to """

        # update config file
        table = np.loadtxt(recpdir+'/sitefinder.log')
        if len(table.shape) == 1:
            table = table[np.newaxis,:]

        # create tmp config file name from original config file
        tmp_config_file = list(os.path.splitext(new_config_file))
        tmp_config_file.insert(1,'_tmp')
        tmp_config_file = ''.join(tmp_config_file)

        # remove section 'SITE' and option site in DOCKING section of config file if exists
        with open(tmp_config_file, 'w') as tmpf:
            with open(config_file, 'r') as newf:
                isdock = False
                sitesection = False
                docksection = False
                for line in newf:
                    # check if still in section SITE*
                    if line.startswith('[SITE'):
                        sitesection = True
                    if sitesection and line.startswith('[') and not line.startswith('[SITE'): # new section has been reached
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
                        tmpf.write('site = ' + ', '.join(['site%s'%int(line[0]) for line in sitetable])+'\n')
                for line in sitetable:
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
        """Write and execute MOE script to localize possible binding sites"""

        # (A) write script
        script_name = 'find_binding_sites.sh'
        moe.write_sitefinder_script(script_name, pdbfile, args)
        os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

        # (B) execute script
        subprocess.check_output("./" + script_name + " &> sitefinder.log", shell=True, executable='/bin/bash')

    def run(self):
        """Run Virtual Screening preparation"""
        parser = self.create_arg_parser()
        args = parser.parse_args()

        tcpu1 = time.time()
        self.initialize(args)

        # prepare virtual screening
        self.prepare_vs(args)
        tcpu2 = time.time()

if __name__ == '__main__':
    PrepDocking().run()
