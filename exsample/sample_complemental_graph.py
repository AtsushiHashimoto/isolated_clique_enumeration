#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('../isoclique')
import adjacency_list as al
#import isolated_cliques as ic

'''
import random


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
    V = [0,1,2,4]
    E = [[0,4],[0,2],[1,4],[2,4]]
    S = [1,2,4]
    indicated_label_order = [4,2,1,0]    
    
    V = [v+1 for v in V]
    E = [[e[0]+1,e[1]+1] for e in E]
    S = [v+1 for v in S]
    indicated_label_order = [5,3,2,1]    
    print(E)  

    print("indicated label order: ",indicated_label_order)
    G = al.AdjacencyList(E,labels=indicated_label_order,do_sort=True,encode_method=encode_method)
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





