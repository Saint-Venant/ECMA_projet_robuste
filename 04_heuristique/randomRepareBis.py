# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 21:50:14 2020

@author: thibaut
"""

import numpy as np
import os
from time import time
from transformation_fichiers import load_data
from classes import make_yki, make_x, make_edges
from SP2 import SP2k, SP1
from graphes_perf import affiche_diag_perf
from matplotlib import pyplot as plt


def init_yki(n, K):
    yki = [n*[0] for k in range(K)]
    for i in range(n):
        k = np.random.randint(K)
        yki[k][i] = 1
    return yki

def compute_SP2k(yki, w_v, W_v, W, n):
    values_SP2k = []
    valMax = 0
    indexMax = -1
    for k in range(K):
        val = SP2k(yki,k,w_v,W_v,W,n)
        values_SP2k.append(val)
        if val > valMax:
            valMax = val
            indexMax = k
    return values_SP2k, valMax, indexMax

def getFeasibleSolution(n,L,W,K,B,w_v,W_v,lh):
    yki = init_yki(n, K)

    # compute the values of all SP2-k sub-problems and get the max
    values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
    print(values_SP2k)
    feasible = (valMax <= B)

    maxIter = 1000
    iteration = 0
    while not(feasible) and (iteration < maxIter):
        # change one vertex to a random cluster
        vertices = [i for i in range(n) if yki[kMax][i] == 1]
        assert(len(vertices) > 0)
        np.random.shuffle(vertices)

        acceptSwitch = False
        for i_change in vertices:
            for new_k in range(K):
                yki[kMax][i_change] = 0
                yki[new_k][i_change] = 1

                _, newValMax, _ = compute_SP2k(yki, w_v, W_v, W, n)

                if (newValMax < valMax):
                    acceptSwitch = True
                    break
                else:
                    yki[kMax][i_change] = 1
                    yki[new_k][i_change] = 0
            if acceptSwitch:
                break

        if not(acceptSwitch):
            # random change
            i_change = vertices[0]
            new_k = np.random.randint(K)
            yki[kMax][i_change] = 0
            yki[new_k][i_change] = 1

        values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
        print(values_SP2k)
        feasible = (valMax <= B)
        iteration += 1

    if feasible:
        print("Found feasible solution\n\n")
    else:
        print("No feasible solution found\n\n")
    return yki

def V1(yki, n, K):
    yyki = [x.copy() for x in yki]
    i = np.random.randint(n)
    for k in range(K):
        yyki[k][i] = 0
    k = np.random.randint(K)
    yyki[k][i] = 1
    return yyki

def isFeasible(yki, n, K, B, w_v, W_v, W):
    values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
    return valMax <= B

def descent(yki,edges,lij,lh,L,n,K,B,w_v,W_v,W):
    x=make_x(yki,edges,K)
    value = SP1(x, edges, lij, L, n, lh)
    for p in range(1000):
        yyki = V1(yki,n,K)
        xx = make_x(yyki,edges,K)
        if isFeasible(yyki,n,K,B,w_v,W_v,W):
            valueBis = SP1(xx, edges, lij, L, n, lh)
            if valueBis <= value:
                yki = yyki
    return yki
        
        
if __name__ == '__main__':
    np.random.seed(2)
    # load data
    #instanceFile = '..\\modified_data\\22_ulysses_6.dat'
    instanceFile = '..\\modified_data\\52_berlin_9.dat'
    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path=instanceFile)
    edges=make_edges(n)
    
    yki = getFeasibleSolution(n,L,W,K,B,w_v,W_v,lh)
    x=make_x(yki,edges,K)

    value = SP1(x, edges, lij, L, n, lh)
    print("value : ", value)

    yki = descent(yki,edges,lij,lh,L,n,K,B,w_v,W_v,W)
    x=make_x(yki,edges,K)
    value = SP1(x, edges, lij, L, n, lh)
    print("value after descent : ", value)
