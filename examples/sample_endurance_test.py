#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('../isoclique')
import isolated_cliques as ics

import math
import time
import numpy as np
from scipy import stats

from sklearn.datasets import make_blobs
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.metrics import f1_score, adjusted_rand_score, adjusted_mutual_info_score

from warnings import warn
from argparse import ArgumentParser


def parser():
    usage = 'Usage: python {} [--verbose] [-N <INT>] [-C <INT>]'.format(__file__) 
    argparser = ArgumentParser(usage=usage)
    # params for control the test.
    argparser.add_argument('-v', '--verbose',
                           action='store_true',
                           help='show verbose message')
    argparser.add_argument('-t', '--trial',type=int,
                           default=1,
                           help='determine how many time the test must be executed.')
    argparser.add_argument('-s', '--stop',
                           action='store_true',
                           help='stop if a result differ from that obtained from blute forth enumeration.')

    # params for Isolated Clique Enumeration
    argparser.add_argument('--logk',
                           action='store_true',
                           help='use "f(k) = c*log(k)" as isolation factor')
    argparser.add_argument('-c', '--isolation-factor',type=float,
                           default=2,
                           help='determin "c" for c-isolation. This is used as a coefficient of log(k) when -logk option is given.')    

    # params for graph generation
    argparser.add_argument('-N', '--num-nodes',type=int,
                           default=36,
                           help='determin # of vertices in a graph.')                         
    argparser.add_argument('-d', '--dim',type=int,
                           default=32,
                           help='determin dimension of sample distributed space.')    
    argparser.add_argument('-C', '--num-communities',type=int,
                           default=10,
                           help='determin # of communities that generate samples in the space.')
    argparser.add_argument('-o', '--outlier-rate',type=float,
                           default=0.2,
                           help='determin rate of isolated vertices (against # of nodes) that do not belong to any communities.')
    argparser.add_argument('-g', '--gamma',type=float,
                           default=-1.0,
                           help='determin gamma for rbf kernel (see func. generate_community_graph). [default 1/N]')

    args = argparser.parse_args()
    # check args
    any_warns = False
    if args.isolation_factor < 0:
        warn("isolation factor must be '>0'.")
        any_warns = True    
    if args.outlier_rate > 1 or 0 > args.outlier_rate:
        warn("outlier_rate must be in range [0,1].")
        any_warns = True
    num_outliers = int(args.num_nodes * args.outlier_rate)

    if args.trial < 1:
        warn("trials must be more than 0")
        any_warns = True

    if args.dim < 1:
        warn("dim must be more than 0")
        any_warns = True

    if args.num_communities + num_outliers > args.num_nodes:
        warn("# of nodes 'N' must be larger or equal to (# of communities + outlier_rate*N)")
        any_warns = True
    
    if any_warns:
        sys.exit()

    if args.gamma <= 0:
        args.gamma = 1.0/args.num_nodes
    return argparser.parse_args()



def generate_community_graph(num_nodes,dim,num_communities,outlier_rate,gamma):
    num_outliers = int(num_nodes*outlier_rate)
    X_out, l = make_blobs(n_samples = num_outliers,
                          n_features = dim,
                          centers=num_outliers,
                          center_box=(-1.0,1.0),
                          shuffle=False)
                            
 #   num_nodes_per_comm = num_nodes//num_communities
    X,labels = make_blobs(n_samples= num_nodes - num_outliers,
                          n_features = dim,
                          centers=num_communities,
                          cluster_std=0.2,
                          center_box=(-1.0,1.0),
                          shuffle=True)
    # concatenate outliers and inliers that belongs to a community.
    X = np.r_[X_out,X]    
    # outlier samples have label '0', otherwise >0.
    labels = np.r_[np.zeros(len(l)),labels+1]

    AffinityMat = rbf_kernel(X,gamma=gamma)
    threshold = stats.scoreatpercentile(AffinityMat, 95) # to use arbitrary point value.
    #threshold = np.median(AffinityMat)
    E = (AffinityMat < threshold).astype(int)
    
    for i in range(num_nodes):
        E[i,i] = 0
    return E, labels

def is_same(cliques1,cliques2):
    c1 = sorted([tuple(sorted(c)) for c in cliques1],key=lambda x:len(x))
    c2 = sorted([tuple(sorted(c)) for c in cliques2],key=lambda x:len(x))
    print("C1!")
    print(c1)
    print("C2!")
    print(c2)
    return c1==c2

def test(args):
    E,labels_gt = generate_community_graph(args.num_nodes,
                                 args.dim,
                                 args.num_communities,
                                 args.outlier_rate,
                                 args.gamma)
    if args.verbose:
        print("[Graph]")
        print(E)
        print("[labels (community)]")
        print(labels_gt)
        
    log = {}
    start = time.time()
    ic_graph = ics.IsolatedCliques(E,
                                   edge_list_format='mat')    

    
    elapsed_time = time.time()-start
    log['time for sorting'] = elapsed_time
    if args.verbose:
        print("time for sorting: ",elapsed_time, "sec.")

    if args.logk:
        def callback(k):
            return args.isolation_factor*math.log(k)
        isolation_factor = -1
    else:
        callback = None
        isolation_factor = args.isolation_factor        

    # DO EFFICIENT SEARCH
    start = time.time()
    pivots, iso_cliques = ic_graph.enumerate(callback=callback,isolation_factor = isolation_factor)
    elapsed_time = time.time()-start    
    log['time for enumeration [EFFICIENT WAY]:'] = elapsed_time
    if args.verbose:
        print("time for enumeration",elapsed_time, "sec.")

    if True or args.verbose:
        print("[obtained cliques]")
        for clique in iso_cliques:            
            print(clique)
            (isolation,deg_ave,deg_min) = ic_graph._evaluate_subgraph(clique)
            print(isolation,deg_ave,deg_min)

    # DO BLUTE FORCE SEARCH
    start = time.time()
    iso_cliques_blute = ic_graph.enumerate_blute(callback=callback,isolation_factor = isolation_factor)
    elapsed_time = time.time()-start    
    log['time for enumeration (in blute force)'] = elapsed_time
    if args.verbose:
        print("time for enumeration [BLUTE FORCE]: ",elapsed_time, "sec.")

    if args.verbose:
        print("[obtaine cliques by (blute force]")
        for clique in iso_cliques_blute:            
            print(clique)
    
    # CHECK THE RESULTS
    log['is the valid result?'] = is_same(iso_cliques, iso_cliques_blute)
    if args.verbose:
        print("Is the results same?: ", log['is the valid result?'])
        
    # EVALUATE as CLUSTERING & OUTLIER DETECTION
    communities = ics.choose_largest(iso_cliques,
                                     args.num_communities,
                                     skip_overlap=False)
    if args.verbose:
        print("communities:",communities)

    labels_est = np.zeros(len(labels_gt),dtype='int32')
    for l, clique in enumerate(communities):
        # offset for outliers
        print(clique)
        for v in clique:
            labels_est[v] = l+1           

    print(labels_est)    

    score = adjusted_rand_score(labels_gt,labels_est)
    log['scores'] = {}
    log['scores']['Adjusted RAND index'] = score
    score= adjusted_mutual_info_score(labels_gt,labels_est)
    log['scores']['Adjusted Mutual Info'] = score
    labels_gt_bin = labels_gt > 0
    labels_est_bin = labels_est > 0
    score= f1_score(labels_gt_bin,labels_est_bin)
    log['scores']['F1 Measure'] = score

    return log
    


def main(args):

    for i in range(args.trial):
        if args.verbose:
            print("start %d th trial."%i)


        result = test(args)
        if args.verbose:
            print(result)



if __name__ == '__main__':
    
    args = parser()
    print(args)
    main(args)

    sys.exit()
    start = time.time()
    ic_graph = ics.IsolatedCliques(E)
    elapsed_time = time.time()-start
    print("%.5f sec. elapsed for graph sorting."%elapsed_time)

    nodes = ic_graph.nodes()
    edges = ic_graph.edges()
    for v,neigh in zip(nodes,edges):
        print(v,": ", neigh)

    isolation_factor = 2
    def callback(k):
        return isolation_factor*math.log(k)



    start = time.time()
#    pivots, iso_cliques = ic_graph.enumerate(isolation_factor=isolation_factor)
    pivots, iso_cliques = ic_graph.enumerate(callback=callback)
    elapsed_time = time.time()-start
    print("%.5f sec. elapsed for enumeration."%elapsed_time)
    

    print("Isolated Cliques")
    for pivot, ic in zip(pivots,iso_cliques):
        stats = ic_graph.evaluate_subgraph(ic)
        print("Pivot: ",pivot, " => [",ic,"]") # ic_graph.decode(ic,1)

#    _ics = ic_graph.enumerate_blute(isolation_factor=isolation_factor, at_most=-1)
    _ics = ic_graph.enumerate_blute(callback=callback, at_most=-1)
    for ic in _ics:
        stats = ic_graph.evaluate_subgraph(ic)
        print(ic) # ic_graph.decode(ic,1)

    

    sys.exit()
