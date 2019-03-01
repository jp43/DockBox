import os
import sys
import argparse
import shutil
import numpy as np
import subprocess
import glob
import time
import pandas as pd

from sklearn import preprocessing
from mdtools.amber import clustr

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
                     program, nposes, firstidx, site = line.strip().split(',')
                     firstidx = int(firstidx)
                     nposes = int(nposes)

                     instance = program
                     if site:
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
        self.extract_results('clustering/info.dat', args)

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
                        results['file'].append(os.path.relpath(self.files_l[index]))

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
            dataset['score_multi'] = dataset[[sf + '_s' for sf in self.scoring_functions]].sum(axis=1)
            features = features + [sf + '_s' for sf in self.scoring_functions] + ['score_multi']

        if hasattr(self, 'rmsd_file'):
            features.append('rmsd')

        if args.np > 1:
            dataset = dataset[dataset['programs'].apply(lambda x: len(x.split(','))>=args.np)]

        self.save_dataset('poses.csv', dataset, features)
        if self.scoring_functions:
            self.select_best_poses('clusters.csv', dataset, features)

        #dataset = dataset[features]
        #dataset.to_csv('poses.csv', index=False)

    def save_dataset(self, filename, dataset, features):

        saved_dataset = dataset.copy()

        # format dataset to make it more user friendly
        if self.scoring_functions:
            cols_to_format = self.scoring_functions + [sf + '_s' for sf in self.scoring_functions] + ['score_multi'] + self.colvar
            saved_dataset[cols_to_format] = saved_dataset[cols_to_format].applymap(lambda x: '{0:.3f}'.format(x))

        saved_dataset = saved_dataset[features]
        saved_dataset.to_csv(filename, index=False)

    def select_best_poses(self, filename, dataset, features):

        f = {'score_multi': 'min', 'population': 'max', 'programs': 'min'}
        for sf in self.scoring_functions:
            f[sf] = 'min'
        clusters = dataset.groupby('cluster_idx').agg(f).reset_index()

        score_multi_norm = np.exp(-clusters['score_multi']).sum()
        clusters['P_score'] = np.exp(-clusters['score_multi'])/score_multi_norm
        clusters['P_pro'] = clusters['programs'].apply(lambda x: len(x.split(',')))/len(self.scoring_functions)
        clusters['P'] = clusters['P_score']*clusters['P_pro']

        cols_to_format = ['P_score', 'P_pro', 'P']
        clusters[cols_to_format] = clusters[cols_to_format].applymap(lambda x: '{0:.3e}'.format(x))

        cols_to_format = self.scoring_functions + ['score_multi']
        clusters[cols_to_format] = clusters[cols_to_format].applymap(lambda x: '{0:.3f}'.format(x))

        clusters = clusters[['cluster_idx', 'population', 'programs'] + self.scoring_functions + ['P_score', 'P_pro', 'P']]
        clusters.to_csv(filename, index=False)

if __name__ == '__main__':
    DockAnalysis().run() 
