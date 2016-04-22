import sys
import os
import argparse
import shutil
import numpy as np
import subprocess
import glob
import time
from amber import clustr

class DockAnalysis(object):

    def initialize(self, args):

        self.dirs = []
        for dir in args.dirs:
            posedir = dir+'/poses'
            if not os.path.isdir(dir):
                raise ValueError('Directory %s not found'%dir)
            elif not os.path.isdir(posedir):
                raise ValueError('Folder %s/poses not found'%posedir)
            else:
                self.dirs.append(os.path.abspath(dir))
        self.ndirs = len(self.dirs)
        self.get_files_location(args)

    def get_files_location(self, args):

        self.files_l = []
        self.files_r = []
        self.instances_l = []

        self.scores = {}

        for jdx, dir in enumerate(self.dirs):
            posedir = dir + '/poses'
            # get location of poses and receptor files
            with open(posedir+'/info.dat') as inff:
                inff.next()
                inff.next()
                for line in inff:
                     instance, nposes, firstidx, site = line.split()
                     firstidx = int(firstidx)
                     nposes = int(nposes)
                     keep_instance = True
                     if args.instances and instance not in args.instances:
                         keep_instance = False
                     if keep_instance:
                         poses_idxs = range(firstidx,firstidx+nposes)
                         for idx in poses_idxs:
                             self.files_l.append(os.path.abspath(posedir+'/lig-%s.mol2'%idx))
                             self.files_r.append(os.path.abspath(dir+'/rec.pdb'))
                             self.instances_l.append(instance)
                         # get scores (if available)
                         rescordir = dir + '/rescoring'
                         if os.path.isdir(rescordir):
                             for filename in glob.glob(rescordir+'/*.score'):
                                 program = os.path.splitext(os.path.basename(filename))[0]
                                 scores = np.loadtxt(filename).tolist()
                                 if program not in self.scores:
                                     self.scores[program] = []
                                 for idx in poses_idxs:
                                     self.scores[program].append(scores[idx-1])

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
            help='RMSD cutoff for clustering. Default: 2.0 ')

        parser.add_argument('-instances',
            type=str,
            dest='instances',
            nargs='+',
            help='Choose instances to be used for consensus docking')

        return parser

    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()

        tcpu1 = time.time()
        print "Starting docking analysis..."

        # initialize parameters
        self.initialize(args)

        # performs clustering
        clustr.do_clustering(self.files_r, self.files_l, cutoff=args.cutoff)
        self.extract_results('clustr/info.dat')

        tcpu2 = time.time()
        print "Analysis procedure done. Total time needed: %i s" %(tcpu2-tcpu1)
   
    def extract_results(self, filename):

        avg_score = []
        heterg = []
        population = []

        ff = open(filename)
        for line in ff:
            # if line does not start with #
            if not line.startswith('#'):
                # indices = numbers of the poses involved in the current cluster
                indices = [i for i, x in enumerate(line.strip()) if x == 'X']
                population.append(len(indices)*100./self.nposes)
                instances = np.array(self.instances_l)[indices].tolist()
                avg = 0
                for key, value in self.scores.items():
                    scores = np.array(value)[indices]
                    avg += np.mean(scores)/np.amin(np.array(value))
                avg /= len(self.scores)
                avg_score.append(avg)
                # compute heterogeneity factor for the current cluster
                heterg.append(len(set(instances))*100./self.ninstances)
            elif line.startswith('#Representative frames:'):
                poses_clustrs_idxs = map(int, line[23:].split())
        ff.close()

        heterg = np.array(heterg)
        population = np.array(population)

        # Option 1
        best_score = 0
        most_populated_clustrs = np.where(population>2)
        heterg_max = np.amax(heterg)
        for idx in most_populated_clustrs[0]:
            if heterg[idx] == heterg_max:
                if avg_score[idx] > best_score:
                    best_cluster_idx = idx
                    best_score = avg_score[idx]

        best_pose_idx = poses_clustrs_idxs[best_cluster_idx]
        print self.files_l[best_pose_idx-1]
        print self.files_r[best_pose_idx-1]
