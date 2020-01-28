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
    yki = [n*[1]]
    for k in range(1, K):
        yki.append(n*[0])
    return yki

#def compute_SP2k(yki, w_v, W_v, W, n):
#    values_SP2k = []
#    valMax = 0
#    indexMax = -1
#    for k in range(K):
#        val = SP2k(yki,k,w_v,W_v,W,n)
#        values_SP2k.append(val)
#        if val > valMax:
#            valMax = val
#            indexMax = k
#    return values_SP2k, valMax, indexMax
def calcule_SP2k(yki, w_v, W_v, W, n,K):
    values_SP2k = []
    valMax = 0
    valMin =0
    indexMax = -1
    indexMin = -1
    for k in range(K):
        val = SP2k(yki,k,w_v,W_v,W,n)
        values_SP2k.append(val)
        if val > valMax:
            valMax = val
            indexMax = k
    indexMin=np.argmin(np.array(values_SP2k))
    valMin=values_SP2k[indexMin]
    return values_SP2k, valMax, indexMax, valMin, indexMin


def heuristique(chemin_fichier):
    #attention les fichiers 80_gr_3.dat 80_gr_6.dat ne marche pas
    np.random.seed(1)
    # load data
    #n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path='modified_data\\10_ulysses_6.dat')
    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path=chemin_fichier)
    edges=make_edges(n)
    
    #yki = init_yki(n, K)
    yki= make_yki(K,n)
    x=make_x(yki,edges,K)

    # compute the values of all SP2-k sub-problems and get the max
    values_SP2k, valMax, kMax, valMin, kMin  = calcule_SP2k(yki, w_v, W_v, W, n,K)
    #print('initialisatoion :\n',values_SP2k)
    feasible = (valMax <= B)

    maxIter = 1000
    iteration = 0
    while not(feasible) and (iteration < maxIter):
        # change one vertex to a random cluster
        vertices = [i for i in range(n) if yki[kMax][i] == 1]
        assert(len(vertices) > 0)
        np.random.shuffle(vertices)
        
        i_change=vertices[np.random.randint(len(vertices))]
        yki[kMax][i_change] = 0
        yki[kMin][i_change] = 1



#
#        acceptSwitch = False
#        for i_change in vertices:
#            for new_k in range(K):
#                yki[kMax][i_change] = 0
#                yki[new_k][i_change] = 1
#
#                _, newValMax, _, _ ,_ = calcule_SP2k(yki, w_v, W_v, W, n)
#
#                if (newValMax < valMax):
#                    acceptSwitch = True
#                    break
#                else:
#                    yki[kMax][i_change] = 1
#                    yki[new_k][i_change] = 0
#            if acceptSwitch:
#                break
#
#        if not(acceptSwitch):
#            # random change
#            i_change = vertices[0]
#            new_k = np.random.randint(K)
#            yki[kMax][i_change] = 0
#            yki[new_k][i_change] = 1
#        
        
        
        

        values_SP2k, valMax, kMax, valMin, kMin  = calcule_SP2k(yki, w_v, W_v, W, n,K)
        #print(values_SP2k)
        feasible = (valMax <= B)
        iteration += 1
    print("_____________________________________________\n",chemin_fichier)
    if feasible:
        print("\n Found feasible solution:\n",values_SP2k)
    else:
        print("No feasible solution found",maxIter)

def resultat_heuristique():
    times_list=[]
    racine='C:\\Users\\thibaut\\Documents\\Ponts_3A\\CNAM\\ECMA\\projet_ECMA_2019_2020\\modified_data' 
    for _, _, files in os.walk(racine):
        for file in files:
            t_debut=time()
            heuristique(chemin_fichier='modified_data\\'+file)
            t_fin=time()
            times_list.append(t_fin-t_debut)
    affiche_diag_perf(times_list)
    plt.title('Diagramme de performances - heuristique')
    plt.xlabel('Temps (s)')
    plt.ylabel("Nombre d'instances résolues par l'heuristique")
    plt.show()
    return

if __name__ == '__main__':
    resultat_heuristique()
    

            
     
        
        
        
        
        
        
        
        
        
        
        
        
        
#if __name__ == '__main__':
#    np.random.seed(2)
#    # load data
#    #n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path='modified_data\\10_ulysses_3.dat')
#    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path='modified_data\\532_att_6.dat')
#    edges=make_edges(n)
#    
#    yki = init_yki(n, K)
#    x=make_x(yki,edges,K)
#
#    # compute the values of all SP2-k sub-problems and get the max
#    values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
#    print(values_SP2k)
#    feasible = (valMax <= B)
#
#    maxIter = 1000
#    iteration = 0
#    while not(feasible) and (iteration < maxIter):
#        # change one vertex to a random cluster
#        vertices = [i for i in range(n) if yki[kMax][i] == 1]
#        assert(len(vertices) > 0)
#        np.random.shuffle(vertices)
#
#        acceptSwitch = False
#        for i_change in vertices:
#            for new_k in range(K):
#                yki[kMax][i_change] = 0
#                yki[new_k][i_change] = 1
#
#                _, newValMax, _ = compute_SP2k(yki, w_v, W_v, W, n)
#
#                if (newValMax < valMax):
#                    acceptSwitch = True
#                    break
#                else:
#                    yki[kMax][i_change] = 1
#                    yki[new_k][i_change] = 0
#            if acceptSwitch:
#                break
#
#        if not(acceptSwitch):
#            # random change
#            i_change = vertices[0]
#            new_k = np.random.randint(K)
#            yki[kMax][i_change] = 0
#            yki[new_k][i_change] = 1
#
#        values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
#        print(values_SP2k)
#        feasible = (valMax <= B)
#        iteration += 1
#
#    if feasible:
#        print("Found feasible solution")
#    else:
#        print("No feasible solution found")