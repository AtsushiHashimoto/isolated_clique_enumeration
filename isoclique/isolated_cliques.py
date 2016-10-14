#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
from warnings import warn
import itertools as iters

class IsolatedCliques():
    def __init__(self, edge_list,node_range='auto',debug_mode=True):
        self.debug_mode = debug_mode
        if node_range != 'auto':
            self.node_range = node_range
        else:
            self.node_range = list(set([item for sublist in edge_list for item in sublist]))
        self.n_nodes = len(self.node_range)

        self.lut_inv = np.zeros(max(self.node_range)+1,dtype='int32')
        self.lut = np.zeros(len(self.node_range),dtype='int32')
        for idx,orig_idx in enumerate(self.node_range):
            self.lut[idx]=orig_idx
            self.lut_inv[orig_idx]=idx
        degrees, self.lut, self.lut_inv = self._sort_nodes(edge_list)
        self.pivot_entries = list(zip(degrees, self._get_neighbors(edge_list,degrees)))

        if self.debug_mode and not self._debug_check_order(range(self.n_nodes)):
            warn("DEBUG: graph is not correctly sorted")
        if self.debug_mode and not self._debug_check_edge_symmetry(range(self.n_nodes)):
            warn("DEBUG: graph edges are not in symmetry.")
            

            

        #for idx,(deg,neigh) in enumerate(self.pivot_entries):
        #    print("orig_id:",self.lut[idx]," degree: ", deg, " neighbors:", [self.lut[idx] for idx in neigh])

    def evaluate_subgraph(self,V,assume_clique=False):
        if assume_clique:
            return self._evaluate_subgraph_assume_clique(V), None, None
        return self._evaluate_subgraph(V)
            
    def enumerate(self,isolation_factor=1.0, callback=None):
        if isolation_factor==1.0:
            ics = self._enumerate1()
        else:
            c_floor = my_floor(isolation_factor)
            ics = self._enumerate_gen(isolation_factor,c_floor)

        return ics
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
                
        return pivots
        
    def _enumerate_gen(self,c,c_floor):
        pivots_idx = []
        pivots = []
        for idx,pivot_entry in enumerate(self.pivot_entries):
            C = self._pivot_trim(idx,*pivot_entry,c, c_floor)
            if not C:               
                continue
            C = self._pivot_enum(idx,*pivot_entry,C,c,c_floor)
            pivots_idx.append(idx)
            pivots.append(C)

        return pivots

    '''
    Utility member functions.
    '''
    def _pivot_enum(self, idx, deg, neigh,C,c,c_floor):
        neigh_ = [n for n in neigh if n<idx]
        c_dd = c_floor - (len(neigh)-len(C))
        C_comp,C_is_clique = self.complement(C)
        if C_is_clique:
            Qs = [C]
        else:
            Qs_comp = vertex_cover(C,C_comp, at_most = c_dd, V_search=C.copy().remove(idx))
            Qs = [[v for v in C if v not in q] for q in Qs_comp]
        return Qs

    def complement(self, C):
        nodes = [[] for c in C]
        is_clique = True
        for i,j in iters.combinations(C,2):
            if j not in self.pivot_entries[i][1]:
                nodes[i].append(j)
                nodes[j].append(i)
                is_clique = False
        return nodes,is_clique
        
    def _pivot_trim(self, idx, deg, neigh, c, c_floor):
        neigh_ = [n for n in neigh if n<idx]
        
        if len(neigh_)>c_floor:
            return False
        C = [idx]+[x for x in neigh if x not in neigh_]
        
        k = len(neigh)
        k_dash = len(C)
        C_h = []
        E_C_h = set([])


        removed = []
        for v in C:
            n_deg= self.pivot_entries[v][0]
            # v  >= (c+1)*k_dash-1 => False
            if n_deg>=(c+1)*k_dash-1:
                C.remove(v)
                ret,removed = remove_check(removed,v,c_floor)
                if not ret:
                    return False
                
        for h,v in enumerate(C):
            n_deg,n_neigh = self.pivot_entries[v]
            #print(h,v)
            # v has more than or equal to c*k_dash outgoing edges => False
            outer_neigh = [x for x in n_neigh if x not in C]
            #print("outer_neigh: ", outer_neigh)
            n_out_edges = len(outer_neigh)
            if n_out_edges>=c*k_dash:
                C.remove(v)
                ret,removed = remove_check(removed,v,c_floor)
                if not ret:
                    return False
                continue
            
            # v has less than k-c_floor adjacent vertices in C => False
            if len(n_neigh)-n_out_edges<(k-c_floor):
                C.remove(v)
                ret,removed = remove_check(removed,v,c_floor)
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
                ret,removed = remove_check(removed,v,c_floor)
                if not ret:
                    return False
                continue
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
            
    def _sort_nodes(self, edge_list):
        degrees = np.zeros(self.n_nodes)
        for edge in edge_list:
            degrees[self.lut_inv[edge[0]]] += 1
            degrees[self.lut_inv[edge[1]]] += 1
        temp_lut = np.zeros(self.n_nodes,dtype='int32')
        temp_lut_inv = np.zeros(self.n_nodes,dtype='int32')
        _degrees = np.zeros(self.n_nodes)
        for rank,(deg,idx) in enumerate(sorted(zip(degrees,range(self.n_nodes)))):
            temp_lut[rank] = self.lut[idx]
            temp_lut_inv[idx] = self.lut_inv[rank]
            _degrees[rank] = deg            
        return _degrees,temp_lut,temp_lut_inv
   
    def _get_neighbors(self, edge_list,degrees):
        neighbors = [[] for i in range(self.n_nodes)]
        for edge in edge_list:
            idx0 = self.lut_inv[edge[0]]
            idx1 = self.lut_inv[edge[1]]
            neighbors[idx0].append(idx1)
            neighbors[idx1].append(idx0)
        
        #[print(idx,": ", "deg=",deg,", ",neighbors[idx]) for idx,(deg,neigh) in enumerate(zip(degree,neighbors))]           
        # sort neighbors
        neighbors = [sorted(neigh,key=lambda x:degrees[x]) for neigh in neighbors]        
        return neighbors

    def evaluate_subgraph_assume_clique(self,V):
        k = len(V)
        n_out_edges = 0
        for v in V:
            n_out_edges += self.pivot_entries[v][0]-(k-1)
        return n_out_edges/k

    def evaluate_subgraph(self,V):
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
        print(k)
        return n_out_edges/k, ave_deg/k, min_deg

    def _debug_check_order(self,V,deep=False):
        if not is_sorted([self.pivot_entries[x][0] for x in V]):
            warn("in _debug_check_order(self,V): pivots V is not sorted in their degree.")
            return False
        if not deep:
            return True
        
        for v in V:
            if not is_sorted([self.pivot_entries[x][0] for x in self.pivot_entries[v][1]]):
                warn("in _debug_check_order(self,V): neighbors of %d is not sorted in their degree."%v)
                return False
        return True
    def _debug_check_edge_symmetry(self,V):
        for v in V:
            for w in self.pivot_entries[v][1]:
                if v not in self.pivot_entries[w][1]:
                    return False
        return True
            
def is_sorted(l):
    if len(l)<2:
        return True
    return all(l[i] <= l[i+1] for i in range(len(l)-1))

def comm_seq(arr_1, arr_2):
    if len(arr_1) == 0 or len(arr_2) == 0:
        return []

    m = len(arr_1) - 1
    n = len(arr_2) - 1

    if arr_1[m] == arr_2[n]:
        return comm_seq(arr_1[:-1], arr_2[:-1]) + [arr_1[m]]

    elif arr_1[m] < arr_2[n]:
        return comm_seq(arr_1, arr_2[:-1])

    elif arr_1[m] > arr_2[n]:
        return comm_seq(arr_1[:-1], arr_2)

def my_floor(flt):
    tmp = int(flt)
    tmp2 = math.floor(flt)
    if tmp==tmp2:
        return tmp-1
    return tmp2

def remove_check(removed,v,c_floor):
    removed.append(v)
    if len(removed)>c_floor:
        return False,removed
    return True,removed

def vertex_cover(V,G, at_most=-1,V_search=None):
    VCs = []
    if not V_search:
        V_search = V
    if at_most<0:
        at_most = len(V_search)

    for k in range(at_most,0,-1):
        for vc_cand in iters.combinations(V_search,k):
            covered_nodes = set(vc_cand+sum([G[x] for x in vc_cand]))
            if len(covered_nodes)<len(V):
                continue
            VCs.append(vc_cand)

    return VCs
        