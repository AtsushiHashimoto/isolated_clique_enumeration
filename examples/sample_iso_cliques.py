#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('../isoclique')
import isoclique as ic
import networkx as nx
import matplotlib.pyplot as plt
import random

import math
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


    start = time.time()
    ic_graph = ic.IsolatedCliques(E)
    elapsed_time = time.time()-start
    print("%.5f sec. elapsed for graph sorting."%elapsed_time)

    nodes = ic_graph.nodes()
    edges = ic_graph.edges()
    for v,neigh in zip(nodes,edges):
        print(v,": ", neigh)

    isolation_factor = 2
    def callback(k):
        return isolation_factor*math.log(k)



    start = time.time()
#    pivots, iso_cliques = ic_graph.enumerate(isolation_factor=isolation_factor)
    pivots, iso_cliques = ic_graph.enumerate(callback=callback)
    elapsed_time = time.time()-start
    print("%.5f sec. elapsed for enumeration."%elapsed_time)
    

    print("Isolated Cliques")
    for pivot, ic in zip(pivots,iso_cliques):
        stats = ic_graph.evaluate_subgraph(ic)
        print("Pivot: ",pivot, " => [",ic,"]") # ic_graph.decode(ic,1)

#    _ics = ic_graph.enumerate_blute(isolation_factor=isolation_factor, at_most=-1)
    _ics = ic_graph.enumerate_blute(callback=callback, at_most=-1)
    for ic in _ics:
        stats = ic_graph.evaluate_subgraph(ic)
        print(ic) # ic_graph.decode(ic,1)

    

    sys.exit()


    # drawing
    rand_colors = generate_random_color_list(len(cliques))
    pos=nx.spring_layout(G) # positions for all nodes
    node_list = set(G.nodes())
    edge_list = set(G.edges())
    for i in range(len(cliques)):
        H = G.subgraph(cliques[i])        
        nx.draw_networkx_nodes(H,pos,
                       nodelist=cliques[i],
                       node_color=rand_colors[i])
        print(H.edges())
        nx.draw_networkx_edges(H,pos,
                               edge_list=H.edges(),
                               edge_color=rand_colors[i],
                               width=4)
        node_list = node_list - set(cliques[i])
        edge_list = edge_list - set(H.edges())
        
    nx.draw_networkx_nodes(H,pos,nodelist=node_list,node_color="#808080")
    nx.draw_networkx_edges(H,pos,edgelist=edge_list,edge_color="#808080")
    plt.show()
    
