#!/usr/bin/env python
# coding: utf-8

import numpy as np
from warnings import warn
import itertools as iters

class LabelEncoder:
    '''
    encode_method = 'list' | 'hash'
    '''
    def __init__(self, labels, encode_method = 'list',parent=None):
        self.encode_method = encode_method

        if encode_method not in ['list','hash']:
            message = "Unknown encode_method: %s"%str(encode_method)
            warn(message)
        
        decode_lut = [None] * len(labels)
 
        for idx,l in enumerate(labels):
            decode_lut[idx] = l
        self.n_labels = len(labels)
        self.set_decode_lut(decode_lut)
        if parent:
            self.set_parent(parent)
        else:
            self.parent = None
    def set_parent(self,parent):
        self.parent = parent
        self.reflect_parent()
        
    def reflect_parent(self):
        labels = [self.parent.decode_lut[self.decode_lut[idx]] for idx in range(self.n_labels)]
        self.set_decode_lut(labels)
    
    def encode(self,val, depth):
        if not self.encode_method:
            return val
        return self._encode(val,depth)
    def decode(self, val, depth):
        if not self.encode_method:
            return val
        return self._decode(val,depth)
    
    def _encode(self,val,depth):
        if depth==0:
            return self.encode_lut[val]
        val = list(val)
        for i,v in enumerate(val):
            val[i] = self._encode(v,depth-1)
        return val
    def _decode(self,val,depth):
        if depth==0:
            return self.decode_lut[val]
        val = list(val)
        for i,v in enumerate(val):            
            val[i] = self._decode(v,depth-1)
        return val
        
    def set_decode_lut(self,decode_lut):
        if len(decode_lut)!=self.n_labels:
            warn("Invaild decode Look-up-table.")
            return
        self.decode_lut = decode_lut
        self._gen_encode_lut()
        self._labels = self.decode_lut

    def _gen_encode_lut(self):
        self.encode_lut = self._gen_encode_lut_(self.decode_lut,self.encode_method)        
        
    def _gen_encode_lut_(self,decode_lut,encode_method):
        if encode_method == 'hash':
            encode_lut = {}
            for idx,l in enumerate(decode_lut):
                encode_lut[l] = idx
        elif encode_method == 'list':
            # check untouched label with initial value '-1'.
            encode_lut = np.zeros(max(decode_lut)+1,dtype='int32') - 1
            for idx,l in enumerate(decode_lut):
                encode_lut[l] = idx
        else:
            message = "Unknown encode_method: %s"%str(encode_method)
            warn(message)
        return encode_lut
    


# class to represent any graph/subgraph by list of node with their neighbors
class _AdjacencyList(LabelEncoder):
    '''
    edge_list_format = 'auto' | 'list' | 'mat' | 'neighbors'
    '''
    def __init__(self, edges, edge_list_format='list',\
                 labels='auto',encode_method='list',\
                 debug_mode=True, do_sort=True):
        self.debug_mode = debug_mode
        
        
        # Graph State Flags
        self.is_sorted = False

        if labels != 'auto':
            labels_ = labels
            self.n_nodes = len(labels_)
        elif edge_list_format == 'mat':
            self.n_nodes = len(edges)
            labels_ = list(range(self.n_nodes))
        else:           
            labels_ = list(set( \
                [item for sublist in edges for item in sublist] \
                ))
            self.n_nodes = len(labels_)
        super(_AdjacencyList,self).__init__(labels_,encode_method=encode_method)

        #print("labels: ",labels_)
        self.set_decode_lut(labels_)
        #print("encode_lut from above labels: ",self.encode_lut)

        _edges = self.encode(edges,2)
        

        if edge_list_format=='list':
            self._adjacency_list = self._get_neighbors_from_list(_edges)
        elif edge_list_format=='mat':
            self._adjacency_list = self._get_neighbors_from_mat(_edges)
        else:
            #print(_edges)
            self._adjacency_list = _edges #self.reorder(_edges)


        if self.debug_mode and not self._debug_check_self_loop():
            warn("DEBUG: graph has self loop.")            
        if self.debug_mode and \
           not self._debug_check_edge_symmetry(range(self.n_nodes)):
            warn("DEBUG: graph edges are not in symmetry.")

    
        # count d(v)
        self._degrees = np.zeros(self.n_nodes,dtype='int32')
        for i in range(self.n_nodes):
            self._degrees[i]=len(self._adjacency_list[i])
        self.n_edges = int(np.sum(self._degrees) / 2)


        if not do_sort:
            return
        # sort adjacency list
        self._sort_nodes()
        
        if self.debug_mode and not self._debug_check_order(range(self.n_nodes),deep=True):
            warn("DEBUG: graph is not correctly sorted")

        self.is_sorted = True

        
        
    def _sort_nodes(self):
        outer_label_set = self.decode_lut

        #self._print("BEFORE adj. list sort")

        # sort adjacency list with degree together
        deg_dec_adjs = sorted(zip(self._degrees,outer_label_set,self._adjacency_list))
        self._degrees = [x[0] for x in deg_dec_adjs]
        self._adjacency_list = [x[2] for x in deg_dec_adjs]

        abs_decode_lut = [x[1] for x in deg_dec_adjs]
        abs_encode_lut = self._gen_encode_lut_(abs_decode_lut, encode_method=self.encode_method)
        #print("abs_decode_lut : ",abs_decode_lut)
        #print("abs_encode_lut: ",abs_encode_lut)
        #print("self.decode_lut: ",self.decode_lut)
        #print("self.encode_lut: ",self.encode_lut)
        

        ## self._print("AFTER adj. list sort") #<- mearningless because of broken degrees


        # encode _adjacency_list and degrees
        for idx,neigh in enumerate(self._adjacency_list):
            '''
            print(idx,": ",neigh)
            print("(pre encode)  ",idx,": ",self._adjacency_list[idx])
            for x in neigh:                
                print(x, "->", self.decode_lut[x], " -> ", abs_encode_lut[self.decode_lut[x]])
            '''
            self._adjacency_list[idx] = [abs_encode_lut[self.decode_lut[x]] for x in neigh] 
            #print("(post encode) ",idx,": ",self._adjacency_list[idx])
        #self._print("AFTER adj. list sorting & all elem. encoding")
    
        # sort each elem in adjacency_list
        self._adjacency_list = [sorted(neigh,key=lambda x:self._degrees[x]) for neigh in self._adjacency_list]        
            
        #self._print("AFTER elem sort")


        # set new decode/encode/labels
        self.set_decode_lut(abs_decode_lut)

    def _subgraph(self, S, do_sort=False):
        adj_list = [[] for v in S]            
        for idx,v in enumerate(S):
            adj_list[idx] = list(set(self._adjacency_list[v]).intersection(S))               
        return adj_list
        
   
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
        for i,j in iters.combinations(range(self.n_nodes),2):
            if edge_mat[i][j]:
                neighbors[i].append(j)
                neighbors[j].append(i)
        return neighbors

            
    def _edges(self, edge_list_format='neighbors'):
        if edge_list_format=='neighbors':
            return self._adjacency_list
        return self.__conv(edge_list_format)
        
    def __conv(self,edge_list_format):
        if edge_list_format=='matrix':
            mat = np.zeros((self.n_nodes,self.n_nodes),dtype='int32')
            for i in range(self.n_nodes):
                for j in self._adjacency_list[i]:
                    mat[i][j] += 1
            return mat

        e = 0
        edges = []
        for i in range(self.n_nodes):
            for j in self._adjacency_list[i]:
                if i > j:
                    edges.append((j,i))
                else:
                    edges.append((i,j))
        return list(set(edges))


    def complement(self,is_dense=True):
        if is_dense:
            return self.__complement_dense_graph()
        warn("AdjcencyList.complement() for general case is not implemented.")
        
    def __complement_dense_graph(self):
        c_adjacency_list = [[] for v in self._labels]
        self.was_clique = True
        for i,j in iters.combinations(range(self.n_nodes),2):
            if j not in self._adjacency_list[i]:
                c_adjacency_list[i].append(j)
                c_adjacency_list[j].append(i)
                self.was_clique = False
        self._adjacency_list = c_adjacency_list
        
    def _debug_check_self_loop(self):
        return all([idx not in neigh for idx,neigh in enumerate(self._adjacency_list)])         

    def _debug_check_edge_symmetry(self,V):
        for v in V:
            for w in self._adjacency_list[v]:
                if v not in self._adjacency_list[w]:
                    return False
        return True

    def _debug_check_order(self,V,deep=False):
        if not is_sorted(self._degrees):
            warn("in _debug_check_order(self,V): pivots V is not sorted in their degree.")
            return False
        if not deep:
            return True
        
        for v in V:
            if not is_sorted([self._degrees[x] for x in self._adjacency_list[v]]):
                warn(("in _debug_check_order(self,V): neighbors of %d is not sorted in their degree."%v))
                print("INNER LABEL")
                for i,neigh in enumerate(self._adjacency_list):
                    print(i, ": ", neigh)
                print("OUTER LABEL")
                for i,neigh in zip(self.nodes(),self.edges()):
                    print(i, ": ", neigh)
                return False
        return True
    def _is_vertex_cover(self,S,edges=None):
        if not edges:
            edges = self.__conv(edge_list_format='list')
        
        for e in edges:
            if e[0] not in S:
                return False
            if e[1] not in S:
                return False
        return True

    def _enumerate_vertex_covers(self, at_most = -1, candidates=None):
        if at_most<0:
            at_most = self.n_nodes
        if not candidates:
            candidates = range(self.n_nodes)
        VCs = []
        edges = self.__conv(edge_list_format='list')
        for size in range(at_most,0,-1):
            for vc in iters.combinations(candidates,size):
                vc = list(vc)
                if self._is_vertex_cover(vc,edges):
                    VCs.append(vc)                
        return VCs
    
    def _print(self,message=""):
        if len(message)>0:
            print(message)
        for v,neigh in enumerate(self._adjacency_list):
            print(v,": ",neigh, " degree=> ("+str(self._degrees[v])+"): ",", ".join(["("+str(self._degrees[n])+")" for n in neigh]))

    def is_clique(self):
        for neigh in self._adjacency_list:
            if len(set(neigh))!=self.n_nodes-1:
                return False
        return True

'''
A Graph that have edges as Neighbors for each graph.
A simple graph is expected.
'''
class AdjacencyList(_AdjacencyList):
    def __init__(self, edges, edge_list_format='list',\
                 labels='auto',encode_method='list',\
                 debug_mode=True, do_sort=True,):
        super(AdjacencyList,self).__init__(edges,edge_list_format=edge_list_format,\
                labels=labels,encode_method=encode_method,\
                debug_mode=debug_mode,do_sort=do_sort,
                )
    def edges(self, edge_list_format='neighbors'):
        return self.decode(self._edges(edge_list_format),2)
        
    def nodes(self):
        return self._labels
    def vertices(self):
        return self._labels     
        
    def subgraph(self, S, do_sort=False,use_in_global=True):
        if use_in_global:
            S_ = self.encode(S,1)
        else:
            S_ = S
        adj_list = self._subgraph(S_, False)
        
        G_S = AdjacencyList(adj_list, edge_list_format='neighbors', \
                             labels=S_,\
                             encode_method='list',debug_mode=self.debug_mode,
                             do_sort=do_sort)
        if use_in_global:
            G_S.set_parent(self)
        return G_S
    def enumerate_vertex_covers(self, at_most = -1, candidates=None):
        if candidates:
            candidates = self.encode(candidates,2)
        VCs = self._enumerate_vertex_covers(at_most,candidates)
        return self._decode(VCs,2)
    def is_vertex_cover(self, S, edges=None):
        return self._is_vertex_cover(self._encode(S,1),edges)
        


# get a mapping look-up-table from a set A to another set B        
def get_map_lut(A,B):
    if len(A)!=len(B):
        warn("get_map_lut(A,B): A and B has different length!")
    f = [None] * (max(A)+1)
    for a,b in zip(A,B):
        f[a] = b
    return f
          
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

def flatten(l):
    return [item for sublist in l for item in sublist]
