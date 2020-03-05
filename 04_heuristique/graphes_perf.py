# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:25:26 2020

@author: thibaut
"""

import os
from matplotlib import pyplot as plt

def recup_temps1(path):
    temps_fichiers=[]
    for root, dirs, files in os.walk(path):
        for file in files:
            fichier = open(path+'\\'+file, 'r')     
            for ligne in fichier:
                if(ligne[0]=='c'):
                    temps_fichiers.append(int(ligne.split(' ')[-1].rstrip('\n')))
    return(temps_fichiers)
    
def recup_temps_dat(path):
    temps_fichiers=[]
    for root, dirs, files in os.walk(path):
        for file in files:
            fichier = open(path+'\\'+file, 'r')     
            for ligne in fichier:
                if(ligne[0]=='s' and ligne[1]=='o'):
                    temps_fichiers.append(  int( float(ligne.split(' ')[-1].rstrip('\n')))       )
    return(temps_fichiers)
    
    
def nombres_sous(tab,nbre):
    compteur=0
    for elt in tab:
        if(elt<=nbre):
            compteur+=1
    return(compteur)
        
def affiche_diag_perf(temps_fichers):
    X=[t for t in range(int(max(temps_fichers)*1.05))]
    Y=[]
    for x in X:
        Y.append(nombres_sous(temps_fichers,x))
    plt.plot(X,Y)
    
        
    return

if __name__=='__main__':
    temps1=recup_temps1(path='simulation1')
    temps2=recup_temps_dat(path='simulation2')
    affiche_diag_perf(temps1)
    affiche_diag_perf(temps2)
    
    plt.title('Diagramme de performances')
    plt.xlabel('Temps (s)')
    plt.ylabel("Nombre d'instances résolues")
    plt.savefig("Diagramme de perf - comparaisons",dpi=150)
    plt.show()