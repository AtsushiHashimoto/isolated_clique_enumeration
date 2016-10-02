#!/usr/bin/env python
# coding: utf-8

import numpy as np


def check_matrix_shape(gmat):
    shape = gmat.shape
    if shape[0]==0:
        raise Exception('empty gmat')
    if shape[0]!=shape[1]:
        raise Exception('non-square gmat')
    return shape[0]

def check_c_validity(c):
    if c<=0:
        raise Exception('negative value is invalid')
        
def sort_descending_order_by_pivot_degree(isocliques):

def get_degrees(gmat):
    
def test_b(clique,degrees)    

def find_isolated_cliques(num_of_nodes, edge_list, c):
    gmat_size = check_matrix_shape(gmat)
    check_c_validity(c)
    
    # initialize vertex labels
    labels = [0]*gmat_size

    # get degree of each vertex
    degrees = get_degrees(gmat)

    # find pivots
    isocliques = {}
    for i in range(gmat_size):
        isoclique = get_neighbor(gmat,i)
        isoclique.append(i)
        if !is_pivots(gmat,c,isoclique):
            continue
        if !test_b(isoclique,degrees):
            continue
        if !clique_test(isoclique,gmat):
            continue
        isocliques[i] = isoclique

    sort_descending_order_by_pivot_degree(isocliques)
    cnum = 0
    for i,clique in isocliques.items():
        if labels[i]>0:
            continue
        cnum = cnum+1
        for j in clique:
            labels[j] = cnum

    return cnum, labels

if __name__ == '__main__':
    gmat = np.matrix(
        [[ 0, 1, 0, 1, 1, 0],
         [ 1, 0, 1, 1, 0, 1],
         [ 0, 1, 0, 0, 1, 1],
         [ 1, 1, 0, 0, 1, 1],
         [ 1, 0, 1, 1, 0, 1],
         [ 0, 1, 1, 1, 1, 0]])
    c = 0.8
    cnum, labels = find_isolated_cliques(gmat,c)

    print("Number of isolated cliques with c=%f: %d" % (c,cnum) )
    print("clique labels for each vertex (0 is non-clique vertices)")
    print(labels)
    
