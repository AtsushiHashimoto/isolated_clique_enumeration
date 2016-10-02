#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import isolated_cliques as ic
import networkx as nx
import matplotlib.pyplot as plt
import random

def generate_random_color():
    return '#{:X}{:X}{:X}'.format(*[random.randint(16, 255) for _ in range(3)])
def generate_random_color_list(num):
    colors = [None]*num
    for i in range(num):
        colors[i] = generate_random_color()
    return colors

if __name__ == '__main__':
    G = nx.karate_club_graph()
    E = G.edges()
    print(E)

    #cliques = [[1,2,3],[5,6,7]]
    cliques = ic.find_isolated_cliques(G.number_of_nodes(),G.edges,c)

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
    
