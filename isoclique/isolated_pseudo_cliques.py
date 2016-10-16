#!/usr/bin/env python
# coding: utf-8

import math
from warnings import warn
import itertools as iters
from isolated_cliques import IsolatedCliques


'''
    cliques: a list of cliques (e.g. [[0,2,6],[1,3,5],...])
    k: the function returns k-largest cliques. When k<=0, len(cliques) is set as k.
       k should be '>0' as long as skip_overlap == False.
    skip_overlap: skip overlapping smaller cliques in the selection.
'''
class IsolatedPseudoCliques(IsolatedCliques):
    def __init__(self, \
                 edges, \
                 edge_list_format='list',\
                 labels='auto', \
                 encode_method='list',\
                 debug_mode=True):
        super(IsolatedPseudoCliques,self).__init__(edges,\
                                             edge_list_format=edge_list_format,\
                                             labels=labels, \
                                             encode_method=encode_method,\
                                             debug_mode=debug_mode,
                                             )

        
    def enumerate_blute(self,isolation_factor=1.0,callback=None,at_most=-1):
        ics = self._enumerate_blute(isolation_factor,callback,at_most)
        return self.decode(ics,2)
        
    def _enumerate_blute(self,isolation_factor=1.0,
                         callback=None, at_most=-1,
                         verbose=False):
        # CAUTION: This may be veryyyyy slow.
        if at_most<0:
            at_most = self.n_nodes
        ics = []
        if verbose:
            print("DO BLUTE SEARCH")
        for N in range(at_most,1,-1):
            if callback:
                isolation_factor = callback(N)
            cand_N = [c for c in range(self.n_nodes) if self._degrees[c]>=N-1]            

            if verbose:
                print(N, ": ", cand_N, " -> search start...")
            inf_ave = N-math.log(N)
            inf_min = N/math.log(N)
            for cand in iters.combinations(cand_N,N):
                # skip a vertex without edges.
                if N<2 and len(self._adjacency_list[cand[0]])==0:
                    continue
                (isolation,deg_ave,deg_min) = self._evaluate_subgraph(cand)
                if not (deg_ave >= inf_ave
                        and deg_min >= inf_min
                        and isolation < isolation_factor):
                    continue
                cand = set(cand)
                if not self._is_maximal_blute(cand,ics):
                    continue
                
                ics.append(list(cand))
                if verbose:
                    print(cand)
                    print(isolation,deg_ave,deg_min)
                
            if verbose:
                print("done")

        if verbose:
            print("BLUTE SEARCH END")
        return [list(c) for c in ics]

    def _is_maximal_blute(self,clique, found_cliques):
        for ref in found_cliques:
            if clique.issubset(ref):
                return False
        return True
    
    def enumerate(self,isolation_factor=1.0, callback=None):
        pivots, ics = self._enumerate(isolation_factor, callback=callback)
        self.decode(pivots,1)
        self.decode(ics,2)
        return self.decode(pivots,1), self.decode(ics,2)

    def _enumerate(self,isolation_factor=1.0, callback=None):
        ipcs = []
        pivots = []
        pivots_uniq = set([])
        for v in range(self.n_nodes):            
            neighbors = self._adjacency_list[v]
            Ss = self._enumerate_trim(v,neighbors,pivots_uniq,isolation_factor)
            if len(Ss)==0:
                continue
            for cand in Ss:

                if not callback:
                    c = isolation_factor
                else:
                    c = callback(len(cand))
                N = len(cand)
                if N < 2:
                    continue
                inf_ave = N-math.log(N)
                inf_min = N/math.log(N)
                (isolation,deg_ave,deg_min) = self._evaluate_subgraph(cand)
                if not (deg_ave >= inf_ave
                        and deg_min >= inf_min
                        and isolation < c):
                    continue
                cand = set(cand)                
                if not self._is_maximal_blute(cand,ipcs):
                    continue
                pivots.append(v)
                ipcs.append(list(cand))
                #pivots_uniq.add(v)
                
        return pivots,ipcs

    def _enumerate_trim(self, idx, neigh, pivots, c):
        N = int(self._degrees[idx])
        if N<=0:
            return []    
        N_start = max(1,int(math.log(N)-1))
        if N<=0:
            return []
        cand = []
        for k in range(N,N_start,-1):
            out_edges_sup = c * math.log(k)
            out_edges = N - k
            if out_edges_sup < out_edges:
                break
            for x in iters.combinations(neigh,k):
                x = set(x)
                if len(x.intersection(pivots))>0:
                    continue
                cand.append([idx] + list(x))
        return cand
                    
