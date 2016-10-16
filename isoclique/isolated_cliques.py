#!/usr/bin/env python
# coding: utf-8

import math
from warnings import warn
import itertools as iters
from adjacency_list import AdjacencyList, comm_seq


'''
    cliques: a list of cliques (e.g. [[0,2,6],[1,3,5],...])
    k: the function returns k-largest cliques. When k<=0, len(cliques) is set as k.
       k should be '>0' as long as skip_overlap == False.
    skip_overlap: skip overlapping smaller cliques in the selection.
'''
def choose_largest(cliques, k=-1, skip_overlap=False):
    cliques = sorted(cliques,key=lambda x: len(x))
    if skip_overlap:
        return cliques[:k] 

    if k<=0:
        k = len(cliques)
    covered = {}
    largests = []
    for clique in cliques:
        if any([v in covered for v in clique]):
            continue
        covered.add(tuple(clique))
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
        self.pivot_entries = list(zip(self._degrees, self._adjacency_list))


    def evaluate_subgraph(self,V,assume_clique=False):
        V = self.encode(V,1)
        if assume_clique:
            return self._evaluate_subgraph_assume_clique(V), None, None
        return self._evaluate_subgraph(V)
        
    def enumerate_blute(self,isolation_factor=1.0,callback=None,at_most=-1):
        ics = self._enumerate_blute(isolation_factor,callback,at_most)
        return self.decode(ics,2)
        
    def _enumerate_blute(self,isolation_factor=1.0,callback=None,at_most=-1):
        # CAUTION: This may be veryyyyy slow.
        if at_most<0:
            at_most = self.n_nodes
        ics = []
        print("DO BLUTE SEARCH")
        print(self._degrees)
        for N in range(at_most,0,-1):
            if callback:
                isolation_factor = callback(N)
            cand_N = [c for c in range(self.n_nodes) if self._degrees[c]>=N-1]            
            #print(N, ": ", cand_N, " -> search start...")
            print(N)
            for cand in iters.combinations(cand_N,N):
                (isolation,deg_ave,deg_min) = self._evaluate_subgraph(cand)
                if deg_ave == N-1 and isolation < isolation_factor:
                    #print(cand, ": ", (isolation,deg_ave,deg_min))
                    ics.append(list(cand))
            #print("done")
        print("BLUTE SEARCH END")
        return ics
       
    def enumerate(self,isolation_factor=1.0, callback=None):
        pivots, ics = self._enumerate(isolation_factor,callback)
        return self.decode(pivots,1),self.decode(ics,2)

    def _enumerate(self,isolation_factor=1.0, callback=None):
        if isolation_factor==1.0 and not callback:
            pivots, ics = self._enumerate1()
        else:
            '''
            c_floor = self.my_floor(isolation_factor)
            pivots, ics = self._enumerate_gen(isolation_factor,c_floor)
            '''
            pivots, ics = self._enumerate_gen(isolation_factor,callback)

        # decode clique members 
        return pivots,ics
        '''
        iso_cliques = []
        for ic in ics:
            iso_cliques.append([self.lut[i] for i in ic])
        return iso_cliques
        '''
    
    def _enumerate1(self):
        # case: isolation_factor=1
        pivots_idx = []
        pivots = []
        for idx, pivot_entry in enumerate(self.pivot_entries):            
            if not self._one_pivot_test_a(idx,*pivot_entry,pivots_idx):
                continue
            if not self._one_pivot_test_b(idx,*pivot_entry):
                continue
            if not self._one_pivot_test_c(idx,*pivot_entry):
                continue
            pivots_idx.append(idx)
            pivots.append((idx,[idx]+pivot_entry[1]))
                
        return pivots_idx, pivots
        
#    def _enumerate_gen(self,c,c_floor):        
    def _enumerate_gen(self,isolation_factor,callback):        
        pivots_idx = []
        pivots_ref = []
        pivots = []     
        for idx,pivot_entry in enumerate(self.pivot_entries):
            if callback:
                isolation_factor = callback(self._degrees[idx])
            c_floor = self.my_floor(isolation_factor)

            C = self._pivot_trim(idx,*pivot_entry,isolation_factor, c_floor)
            if not C:               
                continue
            Qs = self._pivot_enum(idx,*pivot_entry,C,isolation_factor,c_floor)
            if len(Qs)==0:
                continue
            Qs = self._pivot_scr(idx,Qs,pivots_idx,pivots,isolation_factor,pivot_entry[1])
            if len(Qs)>0:
                pivots_idx.append(idx)
                pivots_ref += [idx] * len(Qs)
                pivots += Qs

        return pivots_ref, pivots

    '''
    Utility member functions.
    '''
    def _pivot_scr(self,idx,Qs,pivots_idx,pivots,c,neigh):
        neigh_ = [n for n in neigh if n<idx]
        if len(set(pivots).intersection(neigh_))>0:
            return []

        Qs = set([tuple(q) for q in Qs])
        Qs = [q for q in Qs if self._evaluate_subgraph_assume_clique(q)<c]
        
        
        return list(Qs)
        
    def _pivot_enum(self, idx, deg, neigh,C,c,c_floor):
#        neigh_ = [n for n in neigh if n<idx]
        c_dd = c_floor - (len(neigh)-len(C))
        G_C = self.subgraph(C,do_sort=False,use_in_global=False)
        
        G_C.complement(is_dense=True)
        if G_C.was_clique:
            Qs = [C]
        else:
            VCs = G_C.enumerate_vertex_covers(at_most = c_dd, candidates=C.copy().remove(idx))
            Qs = [[v for v in C if v not in vc] for vc in VCs]            
        return Qs
        
    def _pivot_trim(self, idx, deg, neigh, c, c_floor):
        neigh_ = [n for n in neigh if n<idx]
        
        if len(neigh_)>c_floor:
            return False
        C = [idx]+[x for x in neigh if x not in neigh_]
        
        k = len(neigh)
        k_dash = len(C)
        C_h = []
        E_C_h = set([])

        #print(idx,": ",C)

        removed = []
        for v in C:
            n_deg= self.pivot_entries[v][0]
            # v  >= (c+1)*k_dash-1 => False
            if n_deg>=(c+1)*k_dash-1:
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor)
                if not ret:
                    return False
        #print(idx,": ",C)
                
        for h,v in enumerate(C):
            n_deg,n_neigh = self.pivot_entries[v]
            #print(h,v)
            # v has more than or equal to c*k_dash outgoing edges => False
            outer_neigh = [x for x in n_neigh if x not in C]
            #print("outer_neigh: ", outer_neigh)
            n_out_edges = len(outer_neigh)
            if n_out_edges>=c*k_dash:
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor)
                if not ret:
                    return False
                continue
            
            # v has less than k-c_floor adjacent vertices in C => False
            if len(n_neigh)-n_out_edges<(k-c_floor):
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor)
                if not ret:
                    return False
                continue
            #print(h,v)

            # for h = 1,2,...,k_dash, |E(C_h,neigh_|>=c*(c+1)*h => False
            C_h.append(v)
            E_C_h = E_C_h.union(outer_neigh)
            #print("E(V-C_h): ",E_C_h)
            if len(E_C_h)>= c*(c+1)*(h+1):
                C.remove(v)
                ret,removed = self.remove_check(removed,v,c_floor)
                if not ret:
                    return False
                continue
        #print(idx,": ",C)
        return C
        
    def _one_pivot_test_a(self, idx, deg, neigh, pivots):
        neigh_ = [n for n in neigh if n<idx]
        if not comm_seq(neigh_,pivots):
            return True
        return False
    def _one_pivot_test_b(self,idx,deg,neigh):
        print(self.pivot_entries[idx])
        for n_idx in neigh:
            n_deg = self.pivot_entries[n_idx][0]
            if n_deg > 2*deg-2:
                return False
        return True
    def _one_pivot_test_c(self,idx,deg,neigh):
        C_h = [idx]
        E_C_h = set([]) 
        for h, n_idx in enumerate(neigh):
            n_neigh = self.pivot_entries[n_idx][1]
            n_neigh.remove(idx)
            # subset test
            n_iter = iter(neigh)
            if not all(v in n_iter for v in n_neigh):
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
        n_out_edges = 0
        for v in V:
            n_out_edges += self.pivot_entries[v][0]-(k-1)
        return n_out_edges/k

    def _evaluate_subgraph(self,V):
        if self.debug_mode and not self._debug_check_order(V,False):
            warn("DEBUG: Clique V is not sorted.")
#            warn("V => ["+" ".join(map(str,V)) + "]")
            
        k = len(V)
        min_deg = k
        ave_deg = 0
        n_out_edges = 0
        for v in V:
            deg, neigh = self.pivot_entries[v]
            inner_neigh = comm_seq(V,neigh)
            n_out_edges += len(neigh)-len(inner_neigh)            
            min_deg = min(min_deg,len(inner_neigh))
            ave_deg += len(inner_neigh)
        return n_out_edges/k, ave_deg/k, min_deg


    def my_floor(self, flt):
        tmp = int(flt)
        tmp2 = math.floor(flt)
        if tmp==tmp2:
            return tmp-1
        return tmp2

    def remove_check(self, removed,v,c_floor):
        removed.append(v)
        if len(removed)>c_floor:
            return False,removed
        return True,removed
