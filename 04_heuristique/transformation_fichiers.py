# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 16:08:12 2020

@author: thibaut
"""

import os
import numpy as np


def load_data(path='modified_data\\10_ulysses_3.dat'):
    
    fichier = open(path, 'r')     
    for ligne in fichier:
        if(ligne[0]=='n'):
            n=int(ligne.rstrip(';\n').split(' ')[-1])
        if(ligne[0]=='L'):
            L=int(ligne.rstrip(';\n').split(' ')[-1])
        if(ligne[0]=='W' and ligne[1]==' '):
            W=int(ligne.rstrip(';\n').split(' ')[-1])
        if(ligne[0]=='K'):
            K=int(ligne.rstrip(';\n').split(' ')[-1])
        if(ligne[0]=='B'):
            B=int(ligne.rstrip(';\n').split(' ')[-1])
        if(ligne[0]=='w'):
            lecture_w=(ligne.rstrip(';\n').split('=')[-1]).replace(' ','').replace('[','').replace(']','')
            w_v=[int(x) for x in lecture_w.split(',')]
        if(ligne[0]=='W' and ligne[1]=='_'):
            lecture_W=(ligne.rstrip(';\n').split('=')[-1]).replace(' ','').replace('[','').replace(']','')
            W_v=[float(x) for x in lecture_W.split(',')]
        if(ligne[0]=='l'):
            lecture_lh=(ligne.rstrip(';\n').split('=')[-1]).replace(' ','').replace('[','').replace(']','')
            lh=[int(x) for x in lecture_lh.split(',')]
        if(ligne[0]=='c'):
            coordonnees=[]
            lecture_coord=(ligne.rstrip(';\n').split('=')[-1]).split('[')#.replace(' ','').replace('[','').replace(']','')
            for iter in range(2,len(lecture_coord)):
                coordonnees.append([float(lecture_coord[iter].rstrip('],').split(',')[0]),float(lecture_coord[iter].rstrip('],').split(',')[1])])
            '''coordonnees contient bien les coordonnées comme dans le .dat'''
            lij=np.zeros(shape=(n,n))
            for i in range(n):
                for j in range(n):
                    lij[i][j]=np.sqrt((coordonnees[i][0]-coordonnees[j][0])**2+(coordonnees[i][1]-coordonnees[j][1])**2)
            return(n,L,W,K,B,w_v,W_v,lh,coordonnees,lij)            
        


def transform(input_path,output_path,file):
    chemin_acces=input_path+'\\'+file
    chemin_sortie=output_path+'\\'+file.split('.')[0]+'.dat'
    fichier = open(chemin_acces, 'r') 
    ecriture=open(chemin_sortie,'w')
    coord='coordinates = ['
    
    for ligne in fichier:
        if(ligne[0]in['n','L','W','K','B','w','W','l']):
            #print(ligne.rstrip('\n')+';\n')
            ecriture.write(ligne.rstrip('\n')+';\n')
            pass
        if(ligne[0].isdigit()):
            #print(ligne.split(' '))
            coord=coord+'['+ligne.split(' ')[0]+','+ligne.split(' ')[1]+'],'
    
    coord=coord.rstrip(',')
    coord=coord+'];'
    ecriture.write(coord)
    #print(coord)
    print(chemin_sortie)
    
    


if __name__ == "__main__":
    input_path='C:\\Users\\thibaut\\Documents\\Ponts_3A\\CNAM\\ECMA\\projet_ECMA_2019_2020\\data'
    output_path='C:\\Users\\thibaut\\Documents\\Ponts_3A\\CNAM\\ECMA\\projet_ECMA_2019_2020\\modified_data'
    
    for root, dirs, files in os.walk(input_path):
        for file in files:
            transform(input_path,output_path,file)