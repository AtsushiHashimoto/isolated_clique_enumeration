#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('../isoclique')
import isoclique as ic
import networkx as nx
import random

import time

def generate_random_color():
    return '#{:X}{:X}{:X}'.format(*[random.randint(16, 255) for _ in range(3)])
def generate_random_color_list(num):
    colors = [None]*num
    for i in range(num):
        colors[i] = generate_random_color()
    return colors

if __name__ == '__main__':

    E = nx.karate_club_graph().edges()
    '''
    encode_method=None
    E = [[x,y] for x,y in E]
    '''
    '''
    # generate a graph with node id offset.
    # edge must be a list for irregular node IDs.
    # (e.g. (x,y)->[x+5,y+5])
    E = [[x+5,y+5] for x,y in E]
    encode_method = 'array'
    '''
    #'''
    # for string node labels
    E = [[e1,e2] for e1,e2 in E]
    for i in range(len(E)):
        if E[i][0]==1:
            E[i][0] = 'one'
        elif E[i][1]==1:
            E[i][1] = 'one'
    encode_method = 'hash'
    #'''
            
    #print(E)
    print("WITHOUT_SORT")
    start = time.time()
    print("start to constust AdjacencyList.")
    sorted_graph = ic.AdjacencyList(E,do_sort=False,encode_method=encode_method)
    elapsed_time = time.time()-start
    print("%.5f sec. elapsed for graph sorting."%elapsed_time)

    nodes = sorted_graph.nodes()
    neighbors = sorted_graph.edges()

    for v,neigh in zip(nodes,neighbors):
        print(v,": ", neigh)

    print("WITH_SORT")
    start = time.time()
    sorted_graph = ic.AdjacencyList(E,do_sort=True,encode_method=encode_method)
    elapsed_time = time.time()-start
    print("%.5f sec. elapsed for graph sorting."%elapsed_time)

    nodes = sorted_graph.nodes()
    neighbors = sorted_graph.edges()
    for v,neigh in zip(nodes,neighbors):
        print(v,": ", neigh)

    Gs = sorted_graph.subgraph([4,5,6,7,8])
    print(Gs.edges(edge_list_format='list'))
    Gs.complement()
    print(Gs.edges(edge_list_format='list'))
    

    sys.exit()


