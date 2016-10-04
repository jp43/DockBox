import os
import argparse
import shutil
import numpy as np
import subprocess
import glob
import time

import pandas as pd
from amber import clustr

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
        self.instances_l = []
        self.ranks = []

        self.scores = {}
        nposes_per_instance = {}

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
                             self.sites.append(int(site))
                             if not instance in self.instances_l:
                                 nposes_per_instance[instance] = 1
                             else:
                                 nposes_per_instance[instance] += 1
                             self.ranks.append(nposes_per_instance[instance])
                             self.instances_l.append(instance)
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
            clustr.do_clustering(self.files_r, self.files_l, cutoff=args.cutoff)
        self.extract_results('clustr/info.dat', args)

        tcpu2 = time.time()
        print "Analysis procedure done. Total time needed: %i s" %(tcpu2-tcpu1)
   
    def extract_results(self, filename, args):

        #heterg = []
        #population = []
        #poses = []

        #instances = []
        #instances_no_duplicity = []
        features = ['pose_idx', 'cluster_idx', 'heterogeneity', 'population', 'instance']

        results = {}
        for column in features:
            results[column] = []

        cluster_idx = 0
        with open(filename, 'r') as ff:
            for line in ff:
                if not line.startswith('#'):
                    indices = [i for i, x in enumerate(line.strip()) if x == 'X']
                    population = len(indices)
                    instances_cluster = np.array(self.instances_l)[indices].tolist()
                    nprograms = len(list(set(instances_cluster)))
                    cluster_idx += 1 
                    for index in indices:
                        results['pose_idx'].append(index+1)
                        results['cluster_idx'].append(cluster_idx)
                        results['heterogeneity'].append(nprograms*100./self.ninstances)
                        results['population'].append(population)
                        results['instance'].append(self.instances_l[index])

        dataset = pd.DataFrame(results)
        dataset = dataset.sort_values(['pose_idx'])
        dataset = dataset[features]
        dataset.to_csv('results.csv', index=False)

        #ff = open(filename)
        #for line in ff:
        #    # if line does not start with #
        #    if not line.startswith('#'):
        #        # indices = numbers of the poses involved in the current cluster
        #        indices = [i for i, x in enumerate(line.strip()) if x == 'X']
        #        population.append(len(indices))
        #        poses.append([idx + 1 for idx in indices])
        #        # compute heterogeneity factor for the current cluster
        #        instances_cluster = np.array(self.instances_l)[indices].tolist()
        #        instances.append(instances_cluster)
        #        instances_no_duplicity.append(list(set(instances_cluster)))
        #        heterg.append(len(set(instances_cluster))*100./self.ninstances)
        #    elif line.startswith('#Representative frames:'):
        #        rep_poses_idxs = map(int, line[23:].split())
        #ff.close()

        #rep_poses = [self.files_l[idx-1] for idx in rep_poses_idxs]

        #heterg = np.array(heterg)
        #population = np.array(population, dtype=int)
        #nclusters = heterg.shape[0]

        #rank_by_rank_score = []
        #for idx in range(nclusters):
        #    rbr = 0
        #    known_instances = []
        #    for idx_pose in poses[idx]:
        #        jdx = idx_pose-1
        #        if not self.instances_l[jdx] in known_instances:
        #            rbr += self.ranks[jdx]
        #            known_instances.append(self.instances_l[jdx])
        #    rank_by_rank_score.append(rbr*1./len(known_instances))


        #kdx = 0
        #with open('anlz.out', 'w') as ff:
        #    for idx in range(nclusters):
        #        if len(instances_no_duplicity[idx]) >= args.min_overlap:
        #            kdx += 1
        #            ff.write("Cluster #%i \n"%kdx)
        #            ff.write("Poses predicted by " + ', '.join(instances_no_duplicity[idx]) + '\n')
        #            ff.write("Population: %9.2f (%i/%i) \n"%(population[idx]*100./self.nposes,population[idx],self.nposes))
        #            ff.write("Representative pose: %s\n"%os.path.basename(rep_poses[idx]))
        #            ff.write("Poses: %s\n"%(', '.join(map(str, poses[idx]))))
        #            ff.write("Rank-by-rank score: %9.1f\n\n"%(rank_by_rank_score[idx]))

        #resultdir = 'poses-filtered'
        #shutil.rmtree(resultdir, ignore_errors=True)
        #os.mkdir(resultdir)

        ## copy receptor file
        #shutil.copyfile(self.files_r[0], resultdir+'/rec.pdb')

        ## order poses per binding site
        #sites_rep_poses = [self.sites[idx-1] for idx in rep_poses_idxs]
        #rep_poses_site_ordered = np.argsort(sites_rep_poses)

        #kdx = 0
        #summary = ''
        #cursite = 0 
        #nposes = [] # number of poses involved for each binding site
        #for idx in rep_poses_site_ordered:
        #    # update binding site
        #    if cursite != sites_rep_poses[idx]:
        #        cursite = sites_rep_poses[idx]
        #        nposes.append(kdx+1)
        #    if len(instances_no_duplicity[idx]) >= args.min_overlap:
        #        shutil.copyfile(rep_poses[idx], resultdir+'/lig-%s.mol2'%(kdx+1))
        #        summary += '%30s      %10s       %10s       %10s\n'%(','.join(instances_no_duplicity[idx]), 1, kdx+1, sites_rep_poses[idx])
        #        kdx += 1
        #nposes.append(kdx+1)

        ## write files containing the number of poses
        ## generated by each software
        #with open(resultdir+'/info.dat', 'w') as ff:
        #    ff.write('   '.join(map(str,nposes))+'\n')
        #    ff.write('#                     program           nposes           firstidx             site\n')
        #    ff.write(summary)

