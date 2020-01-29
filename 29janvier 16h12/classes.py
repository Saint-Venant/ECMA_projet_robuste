# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 13:08:20 2020

@author: thibaut
"""

import numpy as np
from transformation_fichiers import load_data



class Edge(object):
    def __init__(self,i=None,j=None):
        self.i=i
        self.j=j
    def __str__(self):
        return('i='+str(self.i)+';j='+str(self.j))
        
    pass

'''la liste des edges'''
def make_edges(n):
    edges=[]
    for i in range(n):
        for j in range(i+1,n):
            edges.append(Edge(i=i,j=j))
    return(edges)

def make_yki(K,n, indication=None):
    np.random.seed(123)
    yki=np.zeros(shape=(K,n))
    #print('indication',indication)
    for i in range (n):
        if(indication==None):
            k=np.random.randint(low=0,high=K)
        else:
            k=indication[i]-1
        yki[k][i]=1
    return(yki)

def make_x(yki,edges,K):
    x=[]
    for e in edges:
        memeCluster=0
        for k in range(K):
            memeCluster+=max(0,yki[k][e.i]+yki[k][e.j]-1)
        x.append( memeCluster )
    return(x)





if __name__=="__main__":
    '''les variables du pb sont des variables globales'''
    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path='modified_data\\10_ulysses_3.dat')
    edges=make_edges(n)
    for e in edges:
        print(e)