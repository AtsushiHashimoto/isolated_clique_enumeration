#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
from warnings import warn
import itertools as iters

class LabelEncoder:
    '''
    encode_method = array | hash | None
    '''
    def __init__(self, labels, encode_method = 'array'):
        self.encode_method = encode_method
        if not encode_method or not labels:
            self.encode_lut = list(range(len(labels)))
            self.decode_lut = list(range(len(labels)))
        else:
            self.decode_lut = [None] * len(labels)
            for idx,l in enumerate(labels):
                self.decode_lut[idx] = l
            self._gen_encode_lut()
        self.labels = self.decode_lut                                
    
    def encode(self,val, dim=1):
        if not self.encode_method:
            return val
        return self._encode(val,dim)
    def decode(self, val, dim=1):
        if not self.encode_method:
            return val
        return self._decode(val,dim)

    def _encode(self,val,dim):
        if dim==0:
            return self.encode_lut[val]
        val = list(val)
        for i,v in enumerate(val):
            val[i] = self._encode(v,dim-1)
        return val
    def _decode(self,val,dim):
        if dim==0:
            return self.decode_lut[val]
        val = list(val)
        for i,v in enumerate(val):            
            val[i] = self._decode(v,dim-1)
        return val

    def _gen_encode_lut(self):
        if self.encode_method == 'array':
            self.encode_lut = np.zeros(max(self.decode_lut)+1,dtype='int32')
            for idx,l in enumerate(self.decode_lut):
                self.encode_lut[l] = idx
        else: # if encode_method == 'hash':
            self.encode_lut = {}
            for idx,l in enumerate(self.decode_lut):
                self.encode_lut[l] = idx

    
    def update_codes(self, _decode_lut, _encode_lut):
        if None == self.encode_method:
            print("empty enocde")
            self.encode_method = 'array'
            
        n_labels = len(self.labels)
        if len(_decode_lut)!=n_labels:
            warn("Invaild decode Look-up-table.")
            return
        if len(_encode_lut)!=n_labels:
            warn("Invaild encode Look-up-table.")
            return

        dl = self.decode_lut.copy()
        for idx,label in enumerate(self.decode_lut):
           dl[_encode_lut[idx]] = label
        self.decode_lut = dl
        self._gen_encode_lut()
        self.labels = self.decode_lut


# class to represent any graph/subgraph by list of node with their neighbors
class AdjacencyList(LabelEncoder):
    '''
    edge_list_format = 'auto' | 'list' | 'matrix' | 'neighbors'
    '''
    def __init__(self, edges, edge_list_format='list', _labels='auto',encode_method='array',debug_mode=True, do_sort=True):
        self.debug_mode = debug_mode
        if _labels != 'auto':
            labels = _labels
        else:
            labels = list(set( \
                [item for sublist in edges for item in sublist] \
                ))

        if labels == list(range(len(_labels))):
            super(AdjacencyList,self).__init__(labels,encode_method=None)
        else:
            super(AdjacencyList,self).__init__(labels,encode_method=encode_method)

        self.n_nodes = len(labels)
        _edges = self.encode(edges,2)

        if edge_list_format=='list':
            self.adjacency_list = self._get_neighbors_from_list(_edges)
        elif edge_list_format=='mat':
           self.adjacency_list = self._get_neighbors_from_mat(_edges)

        if self.debug_mode and \
           not self._debug_check_edge_symmetry(range(self.n_nodes)):
            warn("DEBUG: graph edges are not in symmetry.")

        if not do_sort:
            return
    
        # count d(v)
        self.degrees = np.zeros(self.n_nodes)
        for i in range(self.n_nodes):
            self.degrees[i]=len(self.adjacency_list[i])
        self.n_edges = np.sum(self.degrees) / 2

        # sort adjacency list
        self._sort_nodes()
        if self.debug_mode and not self._debug_check_order(range(self.n_nodes),deep=True):
            warn("DEBUG: graph is not correctly sorted")

    def _sort_nodes(self):
        deg_idx = list(zip(self.degrees,range(self.n_nodes)))
        
        encode_lut = np.zeros(self.n_nodes,dtype='int32')
        decode_lut = np.zeros(self.n_nodes,dtype='int32')
        _degrees = np.zeros(self.n_nodes)
        _neighbors = [[] for i in range(self.n_nodes)]
        for rank,(deg,idx) in enumerate(sorted(deg_idx)):
            decode_lut[rank] = idx
            encode_lut[idx] = rank
            _degrees[rank] = deg
            _neighbors[rank] = sorted(self.adjacency_list[idx],key=lambda x:self.degrees[x])

        self.degrees = _degrees

        # temporary set encode_lut to encode current adjacency list
        self.encode_lut = encode_lut
        self.adjacency_list = self.encode(_neighbors)

        # update lookup tables.
        self.update_codes(decode_lut,encode_lut)
   
    def _get_neighbors_from_list(self, edge_list):
        neighbors = [[] for i in range(self.n_nodes)]
        if self.debug_mode:
            for edge in edge_list:
                if self.debug_mode and len(edge)!=2:
                    warn("DEBUG: edge_list violate its format Mx2 at %s."%str(edge))

        for edge in edge_list:
            if self.debug_mode and len(edge)!=2:
                warn("DEBUG: edge_list is not Mx2")
            neighbors[edge[0]].append(edge[1])
            neighbors[edge[1]].append(edge[0])
        return neighbors
    def _get_neighbors_from_mat(self,edge_mat):
        if self.debug_mode:
            if len(edge_mat) != self.n_nodes:
                warn("DEBUG: edge matrix height %s is invalid. It must be %s"%(len(edge_mat),self.n_nodes))
            for arr in edge_mat:
                if len(arr) != self.n_nodes:
                    warn("DEBUG: an array in edge matrix has invlid array size %s. It must be %s"%(len(arr),self.n_nodes))
                        
        neighbors = [[] for i in range(self.n_nodes)]
        for i,j in iters.combinations(self.n_nodes,2):
            neighbors[i].append[j]
            neighbors[j].append[i]
        return neighbors

            
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

        e = 0
        edges = [None] * self.n_edges * 2
        for i in range(self.n_nodes):
            for j in self.adjacency_list[i]:
                if i > j:
                    edges[e] = [j,i]
                else:
                    edges[e] = [i,j]
                e += 1
        return list(set(edges))

    def complement(self, C):
        nodes = [[] for c in C]
        is_clique = True
        for i,j in iters.combinations(C,2):
            if j not in self.pivot_entries[i][1]:
                nodes[i].append(j)
                nodes[j].append(i)
                is_clique = False
        return nodes,is_clique
                   

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
            for w in self.adjacency_list[v]:
                if v not in self.adjacency_list[w]:
                    return False
        return True

    def _debug_check_order(self,V,deep=False):
        if not is_sorted(self.degrees):
            warn("in _debug_check_order(self,V): pivots V is not sorted in their degree.")
            return False
        if not deep:
            return True
        
        for v in V:
            if not is_sorted([self.degrees[x] for x in self.adjacency_list[v]]):
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
        
