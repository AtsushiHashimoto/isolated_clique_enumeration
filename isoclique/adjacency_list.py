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
    edge_list_format = 'auto' | 'list' | 'matrix' | 'neighbors'
    '''
    def __init__(self, edges, edge_list_format='list',\
                 labels='auto',encode_method='list',\
                 debug_mode=True, do_sort=True):
        self.debug_mode = debug_mode
        
        
        # Graph State Flags
        self.is_sorted = False
        self.is_clique = None

        if labels != 'auto':
            labels_ = labels
        else:
            labels_ = list(set( \
                [item for sublist in edges for item in sublist] \
                ))

        if labels_ == list(range(len(labels))):
            super(_AdjacencyList,self).__init__(labels_,encode_method=None)
        else:
            super(_AdjacencyList,self).__init__(labels_,encode_method=encode_method)

        print("labels: ",labels_)
        self.set_decode_lut(labels_)
        print("encode_lut from above labels: ",self.encode_lut)

        self.n_nodes = len(labels_)
        _edges = self.encode(edges,2)
        

        if edge_list_format=='list':
            self._adjacency_list = self._get_neighbors_from_list(_edges)
        elif edge_list_format=='mat':
            self._adjacency_list = self._get_neighbors_from_mat(_edges)
        else:
            print(_edges)
            self._adjacency_list = _edges #self.reorder(_edges)

        for idx,neigh in enumerate(self._adjacency_list):
            print(idx,": ",neigh)

        if self.debug_mode and not self._debug_check_self_loop():
            warn("DEBUG: graph has self loop.")            
        if self.debug_mode and \
           not self._debug_check_edge_symmetry(range(self.n_nodes)):
            warn("DEBUG: graph edges are not in symmetry.")

    
        # count d(v)
        self._degrees = np.zeros(self.n_nodes,dtype='int32')
        for i in range(self.n_nodes):
            self._degrees[i]=len(self._adjacency_list[i])
        self.n_edges = np.sum(self._degrees) / 2


        if not do_sort:
            return
        # sort adjacency list
        self._sort_nodes()
        
        if self.debug_mode and not self._debug_check_order(range(self.n_nodes),deep=True):
            warn("DEBUG: graph is not correctly sorted")

        self.is_sorted = True

    def _print(self,message=""):
        if len(message)>0:
            print(message)
        for v,neigh in enumerate(self._adjacency_list):
            print(v,": ",neigh, " degree=> ("+str(self._degrees[v])+"): ",", ".join(["("+str(self._degrees[n])+")" for n in neigh]))
        
        
    def _sort_nodes(self):
        outer_label_set = sorted(self.nodes())
        abs_decode_lut = [None] * len(self.decode_lut)
        for idx,l in sorted(enumerate(outer_label_set),key=lambda x:self._degrees[x[0]]):
            abs_decode_lut[idx] = l
 #       abs_decode_lut = [temp[i] for i in sorted(range(self.n_nodes),key=lambda x:self._degrees[x])]
        print("abs_decode_lut: ", abs_decode_lut)
        abs_encode_lut = self._gen_encode_lut_(abs_decode_lut, encode_method='list')
        print("abs_encode_lut: ",abs_encode_lut)


        # sort each elem in adjacency_list
        self._adjacency_list = [sorted(neigh,key=lambda x:self._degrees[x]) for neigh in self._adjacency_list]        

        # sort adjacency list with degree together
        deg_adjs = sorted(zip(self._degrees,self._adjacency_list))
        self._degrees = [x[0] for x in deg_adjs]
        self._adjacency_list = [x[1] for x in deg_adjs]

        # encode _adjacency_list and degrees
        for idx,neigh in enumerate(self._adjacency_list):
            #print(idx,": ",neigh)
            #print("(pre encode)  ",idx,": ",self._adjacency_list[idx])
            #for x in neigh:                
            #    print(x, "->", self.decode_lut[x], " -> ", abs_encode_lut[self.decode_lut[x]])
            self._adjacency_list[idx] = [abs_encode_lut[self.decode_lut[x]] for x in neigh] 
            #print("(post encode) ",idx,": ",self._adjacency_list[idx])
    
            

        self._print("After ENCODE")


        print("self.decode_lut (must be [4,3,2,1]):",self.decode_lut)
        print("self.encode_lut (must be [-1,3,2,1,0]):",self.encode_lut)
        
        print("abs_decode_lut (must be [1,2,3,4]):",abs_decode_lut)
        print("abs_encode_lut (must be [-1,0,1,2,3]):",abs_encode_lut)


        # set new decode/encode/labels
        self.set_decode_lut(abs_decode_lut)

    def _subgraph(self, S, do_sort=False):
        adj_list = [[] for v in S]
        print("S = ", S)
        print("graph G: INNER")
        for v,neigh in enumerate(self._adjacency_list):
            print(v, ": ", neigh)
        
            
        print("subgraph G(",S,"): INNER")
        for idx,v in enumerate(S):
            adj_list[idx] = list(set(self._adjacency_list[v]).intersection(S))
            print(v, ": ",adj_list[idx],"=",S," and ",self._adjacency_list[v])
               
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
        for i,j in iters.combinations(self.n_nodes,2):
            neighbors[i].append[j]
            neighbors[j].append[i]
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
        edges = [None] * self.n_edges * 2
        for i in range(self.n_nodes):
            for j in self._adjacency_list[i]:
                if i > j:
                    edges[e] = [j,i]
                else:
                    edges[e] = [i,j]
                e += 1
        return list(set(edges))


    def complement(self,is_dense=True):
        if is_dense:
            return self.__complement_dense_graph()
        warn("AdjcencyList.complement() for general case is not implemented.")
        
    def __complement_dense_graph(self):
        c_adjacency_list = [[] for v in self._labels]
        self.is_clique = True
        for i,j in iters.combinations(range(self.n_nodes),2):
            if j not in self._adjacency_list[i]:
                c_adjacency_list[i].append(j)
                c_adjacency_list[j].append(i)
                self.is_clique = False
        self._adjacency_list = c_adjacency_list
        self._sort_nodes()
        
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
        
    def subgraph(self, S, do_sort=False):
        S_ = self.encode(S,1)
        print("encode lut: ",self.encode_lut)
        print("outer index: ",S)
        print("inner index: ",S_)
        adj_list = self._subgraph(S_, False)

        print("###########")
        for idx,neigh in zip(S_,adj_list):
            print(idx, ": ", neigh)
        
        G_S = AdjacencyList(adj_list, edge_list_format='neighbors', \
                             labels=S_,\
                             encode_method='list',debug_mode=self.debug_mode,
                             do_sort=do_sort)
        G_S.set_parent(self)
        return G_S

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
        


