import os
import sys
import argparse
import shutil
import numpy as np
import subprocess
import glob
import time

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

    def get_files_location(self, args):

        self.files_l = []
        self.files_r = []
        self.sites = []
        self.programs_l = []
        self.instances_l = []
        self.relative_idxs_l = []
        self.ranks = []
        self.dirs_l = []

        self.scores = {}
        nposes_per_instance = {}

        for jdx, dir in enumerate(self.dirs):
            posedir = dir + '/poses'
            # get location of poses and receptor files
            with open(posedir+'/info.dat') as inff:
                inff.next()
                inff.next()
                for line in inff:
                     program, nposes, firstidx, site = line.split()
                     firstidx = int(firstidx)
                     nposes = int(nposes)
                     if site == 'None':
                         instance = program
                     else:
                         instance = program + '.' + site
                     if not (args.instances and instance not in args.instances):
                         poses_idxs = range(firstidx,firstidx+nposes)
                         for index, idx in enumerate(poses_idxs):
                             self.files_l.append(os.path.abspath(posedir+'/lig-%s.mol2'%idx))
                             self.files_r.append(os.path.abspath(dir+'/rec.pdb'))
                             self.sites.append(site)
                             if not instance in self.instances_l:
                                 nposes_per_instance[instance] = 1
                             else:
                                 nposes_per_instance[instance] += 1
                             self.ranks.append(nposes_per_instance[instance])
                             self.programs_l.append(program)
                             self.instances_l.append(instance)
                             self.relative_idxs_l.append(index)
                             self.dirs_l.append(dir)
        if args.ninstances:
            self.ninstances = args.ninstances
        else:
            self.ninstances = len(set(self.instances_l))
        self.nposes = len(self.files_l)

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="Run docking analysis")

        parser.add_argument('-d',
            type=str,
            dest='dirs',
            nargs='+',
            default=['.'],
            help='Directories used for analysis')

        parser.add_argument('-cutoff',
            type=float,
            dest='cutoff',
            default=2.0,
            help='RMSD cutoff for clustering. Default: 2.0')

        parser.add_argument('-instances',
            type=str,
            dest='instances',
            nargs='+',
            help='Choose instances to be used for consensus docking')

        parser.add_argument('-ni',
            type=int,
            dest='ninstances',
            help='Number of instances')

        parser.add_argument('-overlap',
            type=int,
            dest='min_overlap',
            default=1,
            help='Do not keep small clusters (predicted at least min_overlap softwares)')

        parser.add_argument('-extract-only',
            action='store_true',
            dest='extract_results_only',
            default=False,
            help='Extract results only!!!!')

        parser.add_argument('-cleanup',
            action='store_true',
            dest='cleanup',
            default=False,
            help='Cleanup intermediate files')

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
            clustr.do_clustering(self.files_r, self.files_l, cutoff=args.cutoff, cleanup=False)
        self.extract_results('clustr/info.dat', args)

        tcpu2 = time.time()
        print "Analysis procedure done. Total time needed: %i s" %(tcpu2-tcpu1)
   
    def extract_results(self, filename, args):

        features = ['pose_idx', 'cluster_idx', 'heterogeneity', 'programs', 'population', 'instance', 'score', 'file', 'is_representative']

        results = {}
        for column in features:
            results[column] = []

        cluster_idx = 0
        with open(filename, 'r') as ff:
            for line in ff:
                if not line.startswith('#'):
                    indices = [i for i, x in enumerate(line.strip()) if x == 'X']
                    population = len(indices)
                    instances_cluster = np.array(self.programs_l)[indices].tolist()
                    programs = list(set(instances_cluster))
                    nprograms = len(programs)
                    cluster_idx += 1
                    for index in indices:
                        results['pose_idx'].append(index+1)
                        results['cluster_idx'].append(cluster_idx)
                        results['heterogeneity'].append(nprograms*100./self.ninstances)
                        results['programs'].append(','.join(programs))
                        results['population'].append(population)
                        instance = self.instances_l[index]
                        results['instance'].append(instance)
                        with open(self.dirs_l[index]+'/'+instance+'/score.out', 'r') as sout:
                            for kdx, line in enumerate(sout):
                                if kdx == self.relative_idxs_l[index]:
                                    results['score'].append(line.replace('\n',''))
                        results['file'].append(self.files_l[index])
                elif line.startswith('#Representative frames:'):
                    rep_structures = map(int, line.split()[2:])
                    for index in results['pose_idx']:
                        results['is_representative'].append(index in rep_structures)

        dataset = pd.DataFrame(results)
        dataset = dataset.sort_values(['pose_idx'])
        dataset = dataset[features]
        dataset.to_csv('results.csv', index=False)

if __name__ == '__main__':
    DockAnalysis().run() 
