# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 15:01:00 2020

@author: thibaut
credits: certaines portion de code ont été écrites par Google ORTools
"""
import numpy as np
from ortools.linear_solver import pywraplp
from transformation_fichiers import load_data
from classes import make_yki, make_x,make_edges

def SP2k(yki,k,w_v,W_v,W,n):
    # Instantiate a Glop solver, naming it LinearExample.
    solver = pywraplp.Solver('SP2solver',pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)

    variables=[[]for k in range(n) ]
    objective = solver.Objective()
    for i in range(0, n):
        variables[i] = solver.NumVar(w_v[i], w_v[i]*(1+W_v[i]) , 'w2_v'+str(i))
        objective.SetCoefficient(variables[i], yki[k][i])
#    # Create the constraints, one per nutrient.
#  constraints = [0] * len(nutrients)
#  for i in range(0, len(nutrients)):
#    constraints[i] = solver.Constraint(nutrients[i][1], solver.infinity())
#    for j in range(0, len(data)):
#      constraints[i].SetCoefficient(food[j], data[j][i+3])      
    contrainte=solver.Constraint(-solver.infinity(),W+n)
    for i in range(len(variables)):
        contrainte.SetCoefficient(variables[i],1/(w_v[i]))
    
    objective.SetMaximization()
    solver.Solve()
    
    valeur_objectif=0
    for i in range(n):
#        print(variables[i],variables[i].solution_value())
        valeur_objectif+=variables[i].solution_value()*yki[k][i]
#    print(valeur_objectif)
    return(valeur_objectif)
  
    
    
    
    
def SP1(x,edges,lij,L,w_v,W_v,W,n):
    # Instantiate a Glop solver, naming it LinearExample.
    solver = pywraplp.Solver('SP2solver',pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)

    variables=[[]for k in range(len(x)) ]
    objective = solver.Objective()
    for p in range(len(x)):
        variables[p] = solver.NumVar(lij[edges[p].i][edges[p].j], lij[edges[p].i][edges[p].j]+3*(lh[edges[p].i]+lh[edges[p].j]) , 'l1_ij'+str(p))
        #variables[p] = solver.NumVar(borne inf, borne sup , nom)
        objective.SetCoefficient(variables[p], x[p])
#    # Create the constraints, one per nutrient.
#  constraints = [0] * len(nutrients)
#  for i in range(0, len(nutrients)):
#    constraints[i] = solver.Constraint(nutrients[i][1], solver.infinity())
#    for j in range(0, len(data)):
#      constraints[i].SetCoefficient(food[j], data[j][i+3])
    S=0
    for e in edges:
        S+=lij[e.i][e.j]/(lh[e.i]+lh[e.j])
    contrainte=solver.Constraint(-solver.infinity(),L+S)
    for i in range(len(variables)):
        contrainte.SetCoefficient(variables[i],1/(lh[edges[p].i]+lh[edges[p].j]))
    
    objective.SetMaximization()
    solver.Solve()
    
    valeur_objectif=0
    for p in range(len(variables)):
        print(variables[p],variables[p].solution_value())
        valeur_objectif+=variables[p].solution_value()*x[p]
    print(valeur_objectif,L+S)
    return(valeur_objectif)

    
def prochain(tab,index,n,K):
    if(tab[index]<K):
        tab[index]+=1
    else:
        tab[index]=0
        prochain(tab,index+1,n,K)
    

def naif(n,K):
    indicationY=[0 for k in range(n)]
    index=0
    while(indicationY!=[K for k in range(n)]):
        prochain(indicationY,index,n,K)
        #print(indicationY)
        yki=make_yki(K,n,indicationY)
        reponse=0
        for k in range(K):
            if(SP2k(yki,k,w_v,W_v,W,n)<B):
                reponse+=1
        if(reponse==K):
            print('reponse admissible', yki)
        
            
            
        
    


if __name__=="__main__":
    '''les variables du pb sont des variables globales'''
    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path='modified_data\\10_ulysses_3.dat')
    #yki= make_yki(K,n)
    yki=[[0, 0, 0, 1, 0, 1, 1, 1, 0, 0],
         [1, 1, 1, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 1, 0, 0, 0, 1, 0]]#c'est une sol optimale
#    SP2k(yki,1,w_v,W_v,W,n)
    #naif(n,K)
    edges=make_edges(n)
    x=make_x(yki,edges,K)
    #print('test',edges[0].j)
    SP1(x,edges,lij,L,w_v,W_v,W,n)
    