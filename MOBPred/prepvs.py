#!/usr/bin/python
import os
import sys
import shutil
import argparse
import ConfigParser
from glob import glob
import pandas as pd

from prep import ligprep

class PrepDocking(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Prepare files for Docking or Virtual Screening")

        parser.add_argument('-l',
            type=str,
            dest='input_files_l',
            nargs='+',
            help = 'Ligand coordinate file(s): .sdf, .smi')

        parser.add_argument('-r',
            type=str,
            dest='input_files_r',
            nargs='+',
            help = 'Receptor coordinate file(s): .pdb')

        parser.add_argument('-f',
            type=str,
            dest='config_file',
            help='config file: .ini')

        parser.add_argument('-lpflags',
            type=str,
            default="-ph 7.0 -pht 2.0 -i 2 -s 8 -t 4",
            dest='lpflags',
            help='Ligprep (Schrodinger) flags for ligand preparation. Default: "-ph 7.0 -pht 2.0 -i 2 -s 8 -t 4"')

        parser.add_argument('-ligprep',
            dest='use_ligprep',
            action='store_true',
            default=False,
            help='Prepare compounds using ligprep')

        parser.add_argument('-site',
            dest='site',
            type=str,
            default=None,
            help='Update binding sites info in config file from file')

        parser.add_argument('-noprep',
            dest='noprep',
            action='store_true',
            default=False,
            help='No structure preparation, update directories and files only (used debbuging)')

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

    def cleanup(self):
        """Cleanup directories used for preparation"""
        # remove existing directories
        pass

    def prepare_vs(self, args):
        """Prepare files and directories for Virtual Screening"""

        for dir in glob('lig*'):
            if dir[3:].isdigit():
                shutil.rmtree(dir)

        info_l = []
        for key_l, value_l in self.files_l.iteritems():
            for key_r, value_r in self.files_r.iteritems():
                for idx, file_l in enumerate(value_l['isomers']):
                    workdir = key_l + '/' + key_r + '/isomer' + str(idx+1)
                    os.makedirs(workdir)
                    if self.config_file:
                        self.update_config_file(workdir, args.config_file, key_r, args)
                    # copy the files for ligand and receptor in the corresponding working dir
                    shutil.copyfile(file_l, workdir+'/ligand.mol2')
                    shutil.copyfile(value_r['filename'], workdir+'/protein.pdb')
            info_l.append("%s: %s %s"%(key_l,value_l['name'],value_l['filename']))
        info_l = sorted(info_l)

        info_r = []
        for key_r, value_r in self.files_r.iteritems(): 
            info_r.append("%s: %s"%(key_r,value_r['filename']))
        info_r = sorted(info_r)

        with open('INFO.vs', 'w') as infof:
            infof.write("LIGANDS:\n\n")
            infof.write('\n'.join(info_l))
            infof.write('\n\n')
            infof.write("TARGETS:\n\n")
            infof.write('\n'.join(info_r))
            infof.write('\n')

    def prepare_compounds(self, args):
        curdir = os.getcwd()

        if not args.noprep:
            shutil.rmtree('ligprep', ignore_errors=True)
            os.mkdir('ligprep')
            print "Preparing compounds with ligprep..."

        files_l = {}
        for idx, file_l in enumerate(args.input_files_l):
            label = 'lig' + (3-len(str(idx+1)))*'0' + str(idx+1)
            files_l[label] = {}
            if not os.path.exists(file_l):
                raise ValueError("File %s not found!"%(file_l))
            files_l[label]['filename'] = os.path.abspath(file_l)
            basename, ext = os.path.splitext(file_l)
            with open(file_l) as ff:
               if ext == '.sdf':
                   name_l = ff.next().strip()
               elif ext == '.smi':
                   name_l = ff.next().split()[-1]
            files_l[label]['name'] = name_l

            if not args.noprep:
                dir_l = 'ligprep/' + label
                os.mkdir(dir_l)
                os.chdir(dir_l)
                files_l_mol2 = ligprep.prepare_ligand(files_l[label]['filename'], args.lpflags)
                files_l[label]['isomers'] = files_l_mol2
                os.chdir(curdir)
            else:
                dir_l = 'ligprep/' + label
                files_l_mol2 = glob(dir_l+'/*_prep_*.mol2')
                suffix, ext = os.path.splitext(files_l_mol2[0])
                files_l_mol2_s = []
                for jdx in range(len(files_l_mol2)):
                    files_l_mol2_s.append(os.path.abspath(suffix[:-1]+'%s.mol2'%(jdx+1)))
                files_l[label]['isomers'] = files_l_mol2_s

        self.files_l = files_l

    def prepare_targets(self, args):

        files_r = {}
        if args.input_files_r:
            for idx, file_r in enumerate(args.input_files_r):
                label = 'target' + (3-len(str(idx+1)))*'0' + str(idx+1)
                files_r[label] = {}
                # check if target file exists
                if not os.path.exists(file_r):
                    raise ValueError("File %s not found!"%(file_r))
                files_r[label]['filename'] = os.path.abspath(file_r)
        self.files_r = files_r

    def update_config_file(self, workdir, config_file, label_r, args):

        new_config_file = workdir + '/config.ini'

        if args.site:
            self.update_config_file_site_options(new_config_file, config_file, label_r, args.site)
        else:
            # if site option is not provided, copy the original file
            shutil.copyfile(config_file, new_config_file)

    def update_config_file_site_options(self, new_config_file, config_file, label_r, csvfile):
        """Update binding site parameters in config file"""

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

        df = pd.read_csv(csvfile)
        recidx = label_r[6:]
        rows = df[df['recID']==int(recidx)]
        # add new sections 'SITE' and option site
        with open(tmp_config_file, 'w') as tmpf:
            with open(new_config_file, 'r') as newf:
                for line in newf:
                    tmpf.write(line)
                    if line.startswith('[DOCKING]'):
                        tmpf.write('site = ' + ', '.join(['site%s'%int(row[1]['siteID']) for row in rows.iterrows()])+'\n')
                for row in rows.iterrows():
                    section = 'SITE' + str(int(row[1]['siteID']))
                    center_conf = row[1]['center']
                    boxsize_conf = row[1]['size']

                    newsite_section = """
[%(section)s]
center = %(center_conf)s
boxsize = %(boxsize_conf)s"""% locals()

                    tmpf.write(newsite_section+'\n')
        shutil.move(tmp_config_file, new_config_file)

    def run(self):
        """Run Virtual Screening preparation"""
        parser = self.create_arg_parser()
        args = parser.parse_args()

        self.initialize(args)
        self.prepare_compounds(args)
        self.prepare_targets(args)

        # prepare virtual screening
        self.prepare_vs(args)

if __name__ == '__main__':
    PrepDocking().run()
