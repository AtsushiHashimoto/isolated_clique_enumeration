#!/usr/bin/env python
# coding: utf-8

import math
from warnings import warn
import itertools as iters
from .adjacency_list import AdjacencyList, comm_seq,flatten


'''
    cliques: a list of cliques (e.g. [[0,2,6],[1,3,5],...])
    k: the function returns k-largest cliques. When k<=0, len(cliques) is set as k.
       k should be '>0' as long as skip_overlap == False.
    skip_overlap: skip overlapping smaller cliques in the selection.
'''
def choose_largest(cliques, k=-1, skip_overlap=False):
    cliques = sorted(cliques,key=lambda x: len(x))
    cliques.reverse()
    if skip_overlap:
        return cliques[:k] 

    if k<=0:
        k = len(cliques)
    covered = set([])
    largests = []
    for clique in cliques:
        if any([v in covered for v in clique]):
            continue
        covered = covered.union(clique)
        largests.append(clique)
        if len(largests)==k:
            break
    return largests


class IsolatedCliques(AdjacencyList):
    def __init__(self, \
                 edges, \
                 edge_list_format='list',\
                 labels='auto', \
                 encode_method='list',\
                 debug_mode=True):
        super(IsolatedCliques,self).__init__(edges,\
                                             edge_list_format=edge_list_format,\
                                             labels=labels, \
                                             encode_method=encode_method,\
                                             debug_mode=debug_mode,
                                             do_sort=True,\
                                             )


    def evaluate_subgraph(self,S,assume_clique=False):
        V = self.encode(S,1)
        if assume_clique:
            return self._evaluate_subgraph_assume_clique(V), None, None
        return self._evaluate_subgraph(V)
        
    def enumerate_blute(self,isolation_factor=1.0,callback=None,at_most=-1):
        ics = self._enumerate_blute(isolation_factor,callback,at_most)
        return self.decode(ics,2)
    
    def _enumerate_blute(self,isolation_factor=1.0,callback=None,at_most=-1,verbose=False):
        # CAUTION: This may be veryyyyy slow.
        if at_most<0:
            at_most = self.n_nodes
        ics = []
        if verbose:
            print("DO BLUTE SEARCH")
        for N in range(at_most,0,-1):
            if callback:
                isolation_factor = callback(N)
            cand_N = [c for c in range(self.n_nodes) if self._degrees[c]>=N-1]            

            if verbose:
                print(N, ": ", cand_N, " -> search start...")

            for cand in iters.combinations(cand_N,N):
                # skip a vertex without edges.
                cand = list(cand)
                if N==1 and len(self._adjacency_list[cand[0]])==0:
                    continue
                (isolation,deg_ave,deg_min) = self._evaluate_subgraph(cand)

                if isolation >= isolation_factor:                    
                    continue
                if deg_min < N-1:
                    continue
                cand = set(cand)
                if not self._is_maximal(cand,ics):
                    continue
                print(isolation, "<-", cand)
                ics.append(cand)                

            if verbose:
                print("done")

        if verbose:
            print("BLUTE SEARCH END")
        print(ics)
        return ics

    def _is_maximal(self,clique, found_cliques):
        for ref in found_cliques:
            if clique.issubset(ref):
                return False
        return True
    
    def enumerate(self,isolation_factor=1.0, callback=None):
        pivots, ics = self._enumerate(isolation_factor,callback)
        self.decode(pivots,1)
        self.decode(ics,2)
        return self.decode(pivots,1), self.decode(ics,2)

    def _enumerate(self,isolation_factor=1.0, callback=None):
        
        if isolation_factor==1.0 and not callback:
            #pivots, ics = self._enumerate1() # <= bug?
            pivots, ics = self._enumerate_gen(isolation_factor,callback)
        else:
            pivots, ics = self._enumerate_gen(isolation_factor,callback)

        # decode clique members 
        return pivots,ics
        '''
        iso_cliques = []
        for ic in ics:
            iso_cliques.append([self.lut[i] for i in ic])
        return iso_cliques
        '''

#    def _enumerate_gen(self,c,c_floor):        
    def _enumerate_gen(self,isolation_factor,callback):        
        #pivots_idx = []
        pivots_ics = {}
        #iso_cliques = []
        for idx,pivot_entry in enumerate(zip(self._degrees,self._adjacency_list)):
            # this is not mentioned in the original paper,
            # but calculation for a vertex with no edges is non-sense.
            if self._degrees[idx]==0:
                continue
            if callback:
                isolation_factor = callback(self._degrees[idx])
            c_floor = self.my_floor(isolation_factor)
            C = self._pivot_trim(idx,*pivot_entry,isolation_factor, c_floor)
            if C == False:               
                continue
            Qs = self._pivot_enum(idx,*pivot_entry,C,isolation_factor,c_floor)
            if len(Qs)==0:
                continue
            Qs = self._pivot_scr(
                idx,
                Qs,
                pivots_ics,
                isolation_factor,
                pivot_entry[1])
            if len(Qs)==0:
                continue            
            pivots_ics[idx] = Qs

        return pivots_ics.keys(), flatten(pivots_ics.values())

    '''
    Utility member functions.
    '''
    def _pivot_scr(self,idx,Qs,pivots_ics,c,neigh):
        # compare all IsoCliques whose pivot is in N_ (N_=neigh_
        neigh_ = set([n for n in neigh if self._degrees[n]<=self._degrees[idx]])


        #Qs = [q for q in Qs if self._evaluate_subgraph_assume_clique(q)<c]        
        Qs = [q for q in Qs if self._evaluate_subgraph(q)[0]<c]
        for pivot_j in neigh_.intersection(pivots_ics.keys()):
            for q in Qs:
                q_set = set(q)
                for q_j in pivots_ics[pivot_j]:
                    if q_set.issubset(q_j):
                        Qs.remove(q)
                        break
        return list(Qs)
        
    def _pivot_enum(self, idx, deg, neigh,C,c,c_floor):
#        neigh_ = [n for n in neigh if n<idx]
        c_dash = c_floor - (len(neigh)-len(C))
        G_C = self.subgraph(C,do_sort=False,use_in_global=False)        
        G_C.complement(is_dense=True)

        # if there is no edges in complement graph, return it as a clique.
        if G_C.was_clique:
            Qs = [C]
        else:
            # enumerate vertex covers: now it is in a blute force search.
            VCs = G_C.enumerate_vertex_covers(at_most = c_dash, candidates=C.copy().remove(idx))
            Qs = [[v for v in C if v not in vc] for vc in VCs]            
        return [q for q in Qs if len(q)>0]
        
    def _pivot_trim(self, idx, deg, neigh, c, c_floor):
        neigh_ = [n for n in neigh if self._degrees[n]<self._degrees[idx]]
        
        if len(neigh_)>c_floor:
            return False
        C = [idx]+[x for x in neigh if x not in neigh_]
        
        k = len(neigh)+1
        k_dash = len(C)
        num_neigh_ = len(neigh_)
        C_h = []
        E_C_h = set([])

        #print(idx,": ",C)

        removed = []
        # test (a): d(v_{i_j}) must be less than (c+1)k'-1
        for v in C:
            n_deg= self._degrees[v]
            # v  >= (c+1)*k_dash-1 => False
            if n_deg>=(c+1)*k_dash-1:
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor,num_neigh_)
                if not ret:
                    return False
        #print(idx,": ",C)

        for h,v in enumerate(C):
            if h==0:
                continue
            
            n_deg = self._degrees[v]
            n_neigh = self._adjacency_list[v]
            # test(b): v has more than or equal to c*k_dash outgoing edges => False
            outer_neigh = [x for x in n_neigh if x not in C]
            #print("outer_neigh: ", outer_neigh)
            n_out_edges = len(outer_neigh)
            if n_out_edges>=c*k_dash:
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor,num_neigh_)
                if not ret:
                    return False
                continue
            
            # test (c): v must have at least k-c_floor adjacent vertices in C => False
            if len(n_neigh)-n_out_edges<(k-c_floor):
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor,num_neigh_)
                if not ret:
                    return False
                continue
            #print(h,v)

            # test (d): for h = 1,2,...,k_dash, |E(C_h,neigh_|>=c*(c+1)*h => False
            C_h.append(v)
            E_C_h = E_C_h.union(outer_neigh)
            #print("E(V-C_h): ",E_C_h)
            if len(E_C_h)>= c*(c+1)*(h+1):
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor,num_neigh_)
                if not ret:
                    return False
                continue
        #print(idx,": ",C)
        return C


     
    def _enumerate1(self):
        # case: isolation_factor=1
        pivots_idx = []
        pivots = []
        for idx, pivot_entry in enumerate(zip(self._degrees,self._adjacency_list)):
            if self._degrees[idx]==0:
                # ignore a vertex with no edges.
                continue
            if not self._one_pivot_test_a(idx,*pivot_entry,pivots_idx):
                continue
            if not self._one_pivot_test_b(idx,*pivot_entry):
                continue
            if not self._one_pivot_test_c(idx,*pivot_entry):
                continue
            pivots_idx.append(idx)
            pivots.append([idx]+pivot_entry[1])
                
        return pivots_idx, pivots
               
    def _one_pivot_test_a(self, idx, deg, neigh, pivots):        
        neigh_ = [n for n in neigh if self._degrees[n]<self._degrees[idx]]
        if len(neigh_)>0:
            return False
        return True
        '''
        # test        
        if len(comm_seq(neigh_,pivots))==0:
            return True
        return False
        '''

    def _one_pivot_test_b(self,idx,deg,neigh):
        for n_idx in neigh:
            if self._degrees[n_idx] > 2*deg-2:
                return False
        return True
    def _one_pivot_test_c(self,idx,deg,neigh):
        C_h = [idx]
        E_C_h = set([]) 
        for h, n_idx in enumerate(neigh):
            n_neigh = self._adjacency_list[n_idx]
            n_neigh.remove(idx)
            # subset test
            #n_iter = iter(neigh)
            #if not all(v in n_iter for v in n_neigh):                
            if not set(neigh).issubset(n_neigh):
                return False
            C_h.append(n_idx)
            E_C_h = E_C_h.union([x for x in n_neigh if x not in neigh])
            if len(E_C_h)>h:
                return False
        return True
            

    def _check_c_validity(self,c):
        if c<=0:
            raise Exception('negative value is invalid')
                                                          
    def _evaluate_subgraph_assume_clique(self,V):
        k = len(V)
        n_out_edges = 0.0
        for v in V:
            n_out_edges += max(0,self._degrees[v]-(k-1))
        return n_out_edges/k

    def _evaluate_subgraph(self,S):
        if self.debug_mode and not self._debug_check_order(S,False):
            warn("DEBUG: Clique S is not sorted.")
#            warn("V => ["+" ".join(map(str,S)) + "]")
            
        k = len(S)
        min_deg = k
        ave_deg = 0
        n_out_edges = 0
        for v in S:
            neigh = self._adjacency_list[v]
            inner_neigh = comm_seq(S,neigh)
            n_in = len(inner_neigh)
            n_out_edges += len(neigh)-n_in            
            min_deg = min(min_deg,n_in)
            ave_deg += n_in
        
        return n_out_edges/k, ave_deg/k, min_deg


    def my_floor(self, flt):
        tmp = int(flt)
        tmp2 = math.floor(flt)
        if tmp==tmp2:
            return tmp-1
        return tmp2

    def remove_check(self, removed,v,c_floor, num_neigh_):
        removed.append(v)
        if len(removed)+num_neigh_>c_floor:
            return False,removed
        return True,removed
