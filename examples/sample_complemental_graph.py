#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('../isoclique')
import isoclique as ic


import random
import itertools as iters
'''
def generate_random_color():
    return '#{:X}{:X}{:X}'.format(*[random.randint(16, 255) for _ in range(3)])
def generate_random_color_list(num):
    colors = [None]*num
    for i in range(num):
        colors[i] = generate_random_color()
    return colors
'''

if __name__ == '__main__':

    encode_method='list'
    n_nodes = 10
    density = 0.4
    n_edges = int((n_nodes*(n_nodes-1)/2)*density)
    V = list(range(n_nodes))
    E = list(iters.combinations(V,2))
    random.shuffle(E)
    E = E[:n_edges]
    print(V)
    print(E)  

    random.shuffle(V)
    n_subgraph = int(n_nodes/2)
    S = V[:n_subgraph]

    random.shuffle(V) 

    print("indicated label order: ",V)
    G = ic.AdjacencyList(E,labels=V,do_sort=True,encode_method=encode_method)
    nodes = G.nodes()
    edges = G.edges()
    print("G(V,E):")
    for v,neigh in zip(nodes,edges):
        print(v,": ", neigh)
        
    #sys.exit()
        

    subgraph = G.subgraph(S,do_sort=False)
    nodes = subgraph.nodes()
    edges = subgraph.edges()
    print("\nG(S,E(S)):")
    for v,neigh in zip(nodes,edges):
        print(v,": ", neigh)
    
    
    print("complementary subgraph")    
    subgraph.complement()
    # complement calculation judges is 
    print("Is subgraph G(S) a clique?: ",subgraph.was_clique)
    nodes = subgraph.nodes()
    edges = subgraph.edges()
    for v,neigh in zip(nodes,edges):
        print(v,": ", neigh)





