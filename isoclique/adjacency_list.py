#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
from warnings import warn
import itertools as iters

class LabelEncoder():
    '''
    encode_method = array | hash | None
    '''
    def __init__(self, labels, encode_method = 'array'):
        if not encode_method or not labels:
            self.encode_lut = None
            self.decode_lut = None
        else:
            self.decode_lut = [l for l in labels]  
            if encode_method == 'array':
                self.encode_lut = [0]*len(self.decode_lut)
                for l,idx in self.decode_lut:
                    self.encode_lut[l] = idx
            else:
                self.encode_lut = {key:val for val,key in enumerate(labels)}
            
    def encode(self,val, dim=1):
        if not self.encode_lut:
            return val
        return self._encode(val,dim)
    def _encode(self,val,dim):
        if dim==0:
            return self.encode_lut(val)
        for i,v in enumerate(val):
            array[i] = self.encode(v,dim-1)
        return array

    def decode(self,val, dim=1):
        if not self.decode_lut:
            return val
        return self._decode()

    def _decode(self,val,dim):
        if dim==0:
            return self.decode_lut(val)
        for i,v in enumerate(val):
            array[i] = self.decode(v,dim-1)
        return array

# class to represent any graph/subgraph by list of node with their neighbors
class AdjacencyList(ListEncoder):
    '''
    edge_list_format = 'auto' | 'list' | 'matrix' | 'neighbors'
    '''
    def __init__(self, edges, edge_list_format='list', node_range='auto',encode_method='array',debug_mode=True):
        self.debug_mode = debug_mode
        if node_range != 'auto':
            labels = node_range
        else:
            labels = list(set( \
                [item for sublist in edge_list for item in sublist] \
                ))
        if labels == list(range(len(labels))):
            super(AdjacencyList,self).__init__()
        else:
            super(AdjacencyList,self).__init__(labels,encode_method)

        self.n_nodes = len(self.node_range)
        if edge_list_format=='list':
            self.adjacency_list = self._get_neighbors_from_list(edges)
        elif edge_list_format=='mat':
           self.adjacency_list = self._get_neighbors_from_mat(edges)

        if self.debug_mode and \
           not self._debug_check_edge_symmetry(range(self.n_nodes)):
            warn("DEBUG: graph edges are not in symmetry.")

        # count d(v)
        self.degrees = np.array((self.n_nodes))
        [self.degrees[i] = len(self.adjacency_list[i]) \
         for i in range(self.degrees)]

        # sort adjacency list
        self._sort_nodes(edge_list)
        if self.debug_mode and not self._debug_check_order(range(self.n_nodes)):
            warn("DEBUG: graph is not correctly sorted")

            
    def edges(self, edge_list_format='neighbors'):
        if edge_list_format=='neighbors':
            return self.adjacency_list
        return self._conv(edge_list_format)
    def _conv(edge_list_format):
        if edge_list_format=='matrix':
            mat = np.zeros((self.n_nodes,self.n_nodes),dtype='int32')
            for i in range(self.n_nodes):
                for j in self.adjacency_list[i]:
                    mat[i][j] += 1
            return mat
        edges = [None] * 
        for i in in range(self.n_nodes):
            for j in self.adjacency_list[i]:
                edges.append
            

    def complement(self, C):
        nodes = [[] for c in C]
        is_clique = True
        for i,j in iters.combinations(C,2):
            if j not in self.pivot_entries[i][1]:
                nodes[i].append(j)
                nodes[j].append(i)
                is_clique = False
        return nodes,is_clique
                   
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
        
