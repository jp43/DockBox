import os
import sys
import argparse
import shutil
import numpy as np
import subprocess
import glob
import time

from sklearn import preprocessing

import pandas as pd
from MOBPred.amber import clustr

class DockAnalysis(object):

    def initialize(self, args):

        self.dirs = []
        for dir in args.dirs:
            posedir = dir+'/poses'
            if not os.path.isdir(dir):
                raise ValueError('Directory %s not found'%dir)
            elif not os.path.isdir(posedir):
                raise ValueError('Folder %s not found'%posedir)
            else:
                self.dirs.append(os.path.abspath(dir))
        self.ndirs = len(self.dirs)
        self.get_files_location(args)

        self.scoring_functions = []
        if args.scoring_functions:
            self.scoring_functions = args.scoring_functions 

        self.colvar = []
        if args.colvar:
            self.colvar = args.colvar

        if args.rmsd_file:
            self.rmsd_file = args.rmsd_file

    def get_files_location(self, args):

        self.files_l = []
        self.files_r = []

        self.sites = []
        self.programs = []

        self.instances = []
        self.dirs_l = []

        # needed to get results of rescoring
        self.indices_per_dir = []
        self.indices_per_instance = []

        nposes_per_instance = {}

        for jdx, dir in enumerate(self.dirs):
            posedir = dir + '/poses'

            # get location of poses and receptor files
            with open(posedir+'/info.dat') as inff:
                inff.next()
                inff.next()

                index_per_dir = 0
                for line in inff:
                     program, nposes, firstidx, site = line.split()
                     firstidx = int(firstidx)
                     nposes = int(nposes)

                     instance = program
                     if site != 'None':
                         instance += '.' + site

                     if not (args.instances and instance not in args.instances):
                         poses_idxs = range(firstidx, firstidx+nposes)

                         for index, idx in enumerate(poses_idxs):
                             self.files_l.append(os.path.abspath(posedir+'/lig-%s.mol2'%idx))
                             self.files_r.append(os.path.abspath(posedir+'/rec.pdb'))
                             self.sites.append(site)

                             self.programs.append(program)
                             self.instances.append(instance)
                             self.indices_per_instance.append(index)

                             self.dirs_l.append(dir)
                             self.indices_per_dir.append(idx-1)

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="Run docking analysis")

        parser.add_argument('-w',
            type=str,
            dest='dirs',
            nargs='+',
            default=['.'],
            help='Working directories used for analysis')

        parser.add_argument('-rmsd',
            type=float,
            dest='rmsd',
            default=2.0,
            help='RMSD cutoff for clustering. Default: 2.0')

        parser.add_argument('-instances',
            type=str,
            dest='instances',
            nargs='+',
            help='Choose instances to be used for consensus docking')

        parser.add_argument('-np',
            type=int,
            dest='np',
            default=1,
            help='Keep only clusters predicted by at least np softwares)')

        parser.add_argument('-s',
            nargs='+',
            dest='scoring_functions',
            type=str)

        parser.add_argument('-colvar',
            nargs='+',
            dest='colvar',
            type=str)

        parser.add_argument('-extract_only',
            action='store_true',
            dest='extract_results_only',
            default=False,
            help='Extract results only!!!!')

        parser.add_argument('-cleanup',
            action='store_true',
            dest='cleanup',
            default=False,
            help='Cleanup intermediate files')

        parser.add_argument('-add_rmsd',
            dest='rmsd_file',
            type=str)

        return parser

    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()

        tcpu1 = time.time()
        print "Starting docking analysis..."

        # initialize parameters
        self.initialize(args)

        if not args.extract_results_only:
            # performs clustering
            clustr.do_clustering(self.files_r, self.files_l, cutoff=args.rmsd, cleanup=True)
        self.extract_results('clustr/info.dat', args)

        tcpu2 = time.time()
        print "Analysis procedure done. Total time needed: %i s" %(tcpu2-tcpu1)
   
    def extract_results(self, filename, args):

        features = ['pose_idx', 'cluster_idx', 'programs', 'population', 'instance', 'score', 'file', 'is_representative']
        features.extend(self.scoring_functions+self.colvar)

        results = {}
        for column in features:
            results[column] = []

        cluster_idx = 0
        with open(filename, 'r') as ff:
            for line in ff:
                if not line.startswith('#'):
                    indices = [i for i, x in enumerate(line.strip()) if x == 'X']
                    population = len(indices)

                    instances_cluster = np.array(self.programs)[indices].tolist()
                    programs = list(set(instances_cluster))
                    nprograms = len(programs)
                    cluster_idx += 1

                    for index in indices:
                        results['pose_idx'].append(int(index+1))
                        results['cluster_idx'].append(int(cluster_idx))
                        results['programs'].append(','.join(programs))
                        results['population'].append(population)

                        instance = self.instances[index]
                        results['instance'].append(instance)

                        with open(self.dirs_l[index]+'/'+instance+'/score.out', 'r') as sout:
                            for kdx, line in enumerate(sout):
                                if kdx == self.indices_per_instance[index]:
                                    results['score'].append(line.replace('\n',''))
                        
                        for sf in self.scoring_functions + self.colvar:
                            with open(self.dirs_l[index]+'/rescoring/'+sf+'.score', 'r') as sout:
                                for kdx, line in enumerate(sout):
                                    if kdx == self.indices_per_dir[index]:
                                        results[sf].append(float(line.replace('\n','')))
                        results['file'].append(self.files_l[index])

                elif line.startswith('#Representative frames:'):
                    rep_structures = map(int, line.split()[2:])
                    for index in results['pose_idx']:
                        results['is_representative'].append(index in rep_structures)
 
        dataset = pd.DataFrame(results)
        dataset = dataset.sort_values('pose_idx')

        if hasattr(self, 'rmsd_file'):
            rmsd = np.loadtxt(self.rmsd_file)
            dataset['rmsd'] = rmsd[1:,1]

        if self.scoring_functions:
            # get rid of nan values
            condition_on_score = ' and '.join([name + ' < 0 ' for name in self.scoring_functions if not np.isnan(dataset[name].values).all()])
            dataset = dataset.query(condition_on_score).copy()
            # compute rescaled scores
            for sf in self.scoring_functions:
                if np.isnan(dataset[sf].values).all():
                    dataset[sf + '_s'] = dataset[sf]
                else:
                    dataset[sf + '_s'] = preprocessing.scale(dataset[sf].apply(lambda x: x))
            dataset['score_multi_s'] = dataset[[sf + '_s' for sf in self.scoring_functions]].sum(axis=1)
            cols = features + [sf + '_s' for sf in self.scoring_functions] + ['score_multi_s']

        if hasattr(self, 'rmsd_file'):
            cols.append('rmsd')

        if args.np > 1:
            dataset = dataset[dataset['programs'].apply(lambda x: len(x.split(','))>=args.np)]

        clusters_best_score = {}
        if self.scoring_functions:
            f = {'score_multi_s': 'min'}
            #dataset_big_clusters = dataset[dataset['cluster_idx']<=10]
            df_groupby = dataset.groupby('cluster_idx')
            df_clusters = df_groupby.agg(f).reset_index()
            df_clusters_best_score = df_clusters.sort_values('score_multi_s').head(5)
            for idx, row in enumerate(df_clusters_best_score.iterrows()):
                clusters_best_score[int(row[1]['cluster_idx'])] = (idx+1, row[1]['score_multi_s'])
        cwd = os.getcwd()

        #shutil.rmtree('rep_poses_big_clusters', ignore_errors=True)
        #os.mkdir('rep_poses_big_clusters')

        shutil.rmtree('rep_poses_best_scores', ignore_errors=True)
        os.mkdir('rep_poses_best_scores')

        for name, group in dataset.groupby('cluster_idx'):
            #if int(name) <= 5:
            #    dirname = 'big_cluster_%s'%int(name)
            #    shutil.rmtree(dirname, ignore_errors=True)
            #    os.mkdir(dirname)
            #    for row in group.iterrows():
            #        filename = row[1]['file']
            #        filename_rel = os.path.relpath(filename, cwd)
            #        filename_rel_l = filename_rel.split('/')
            #        filename_rel_l.remove('poses')
            #        filename_rel_r = '/'.join(filename_rel_l[:-1] + ['rec.pdb'])
            #        new_filename_rel_l = '-'.join(filename_rel_l) 
            #        new_filename_l = dirname  + '/' + new_filename_rel_l
            #        shutil.copyfile(filename, new_filename_l)
            #        if row[1]['is_representative']:
            #            shutil.copyfile(filename, 'rep_poses_big_clusters/' + new_filename_rel_l)
            #            if os.path.isfile(filename_rel_r):
            #                new_filename_rel_l_suf, ext = os.path.splitext(new_filename_rel_l)
            #                shutil.copyfile(filename_rel_r, 'rep_poses_big_clusters/'+new_filename_rel_l_suf+'.pdb')

            if name in clusters_best_score:
                dirname = 'cluster_best_score_%s'%clusters_best_score[name][0]
                shutil.rmtree(dirname, ignore_errors=True)
                os.mkdir(dirname)
                min_score_str = "%.3f"%clusters_best_score[name][1]

                for row in group.iterrows():
                    filename = row[1]['file']
                    filename_rel = os.path.relpath(filename, cwd)
                    filename_rel_l = filename_rel.split('/')
                    filename_rel_l.remove('poses')
                    filename_rel_r = '/'.join(filename_rel_l[:-1] + ['rec.pdb'])

                    new_filename_rel_l = '-'.join(filename_rel_l)
                    new_filename_l = dirname  + '/' + new_filename_rel_l
                    shutil.copyfile(filename, new_filename_l)

                    if clusters_best_score[name][1] == row[1]['score_multi_s']:
                        shutil.copyfile(filename, dirname + '/rep_pose.mol2')
                        population = row[1]['population']
                        programs = ', '.join(row[1]['programs'].split(','))
                        info = """Directory name: %(dirname)s
Cluster index: %(name)s
Population: %(population)s
Programs: %(programs)s
Best score: %(min_score_str)s
Representative pose (rep_pose.mol2): %(filename)s\n"""%locals()
                        if os.path.isfile(filename_rel_r):
                            filename_r = os.path.abspath(filename_rel_r)
                            new_filename_rel_l_suf, ext = os.path.splitext(new_filename_rel_l)
                            shutil.copyfile(filename_rel_r, dirname + '/rec.pdb')
                            info += "Receptor file (rec.pdb): %(filename_r)s\n"%locals()
                        with open(dirname+'/info.dat', 'w') as ff:
                            ff.write(info)

        if self.scoring_functions:
            cols_to_format = self.scoring_functions + [sf + '_s' for sf in self.scoring_functions] + ['score_multi_s'] + self.colvar
            dataset[cols_to_format] = dataset[cols_to_format].applymap(lambda x: '{0:.3f}'.format(x))

        dataset = dataset[cols]
        dataset.to_csv('poses.csv', index=False)

if __name__ == '__main__':
    DockAnalysis().run() 
