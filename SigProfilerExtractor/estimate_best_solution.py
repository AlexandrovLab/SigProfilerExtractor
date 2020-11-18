#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 16:01:59 2020

@author: mishugeb
"""


import pandas as pd
import numpy as np
import os
from SigProfilerExtractor import subroutines as sub


def estimate_solution(base_csvfile="All_solutions_stat.csv", 
                    All_solution="All_Solutions", 
                    genomes="Samples.txt", 
                    output="results", 
                    title="Selection_Plot",
                    stability=0.8, 
                    min_stability=0.2, 
                    combined_stability=1.0,
                    statistics=True,
                    select=None):
    
    
    
    base_csvfile = pd.read_csv(base_csvfile, sep=",", index_col=0)
    signatures=list(base_csvfile.index)
    genomes=pd.read_csv(genomes, sep="\t", index_col=0)
    colnames=genomes.columns
    genomes=np.array(genomes)
    all_similarities_list=[]
    layer_directory=output
    
    if genomes.shape[0]==78:
        mtype="DBS78"
    elif genomes.shape[0]==83:
        mtype="ID83"
    elif genomes.shape[0]==48:
        mtype="CNV48"
    else:
        mtype="SBS"+str(genomes.shape[0])
    
    #prepare the csvfile
    csvfile=np.zeros([len(signatures),4])
    csvfile=csvfile
    for i in range(len(signatures)):
        if base_csvfile.shape[1]!=3:
            signatures[i]=signatures[i].rstrip("*")
            fnorm=base_csvfile.iloc[i,4].rstrip("%")
            csvfile[i,[1,3]]=base_csvfile.iloc[i,[1,0]]
        elif base_csvfile.shape[1]==3:
            signatures[i]=str(signatures[i])
            fnorm=base_csvfile.iloc[i,1].astype(float)
            fnorm=fnorm*100
            csvfile[i,[1,3]]=base_csvfile.iloc[i,[0,2]]
        csvfile[i,0]=signatures[i]    
        csvfile[i,2]=fnorm
        w=pd.read_csv(All_solution+"/"+mtype+"_"+signatures[i]+"_Signatures/Signatures/"+mtype+"_S"+signatures[i]+"_Signatures.txt", sep="\t", index_col=0)
        h=pd.read_csv(All_solution+"/"+mtype+"_"+signatures[i]+"_Signatures/Activities/"+mtype+"_S"+signatures[i]+"_NMF_Activities.txt", sep="\t", index_col=0)
        w = np.array(w)
        h = np.array(h).T
        est_genomes=np.dot(w,h)
        all_similarities, cosine_similarities = sub.calculate_similarities(genomes, est_genomes, colnames)
        all_similarities_list.append(all_similarities)
        
    csvfile=pd.DataFrame(csvfile)
    csvfile.columns=["Total Signatures",  "Stability",  "Matrix Frobenius%",  "avgStability"]
    
    try:
        if not os.path.exists(layer_directory):
            os.makedirs(layer_directory)
    except: 
        print ("The {} folder could not be created".format("output"))
        
    solution, all_stats= sub.stabVsRError(csvfile, 
                                          layer_directory, 
                                          title, all_similarities_list, 
                                          input_type="dataframe",  
                                          stability=stability, 
                                          min_stability=min_stability, 
                                          combined_stability=combined_stability, 
                                          mtype=mtype,
                                          statistics=statistics,
                                          select=select)
    
    all_stats.insert(1, 'Stability (Avg Silhouette)', csvfile["avgStability"]) 
    all_stats=all_stats.set_index(["Signatures"])
    all_stats.to_csv(layer_directory+"/All_solutions_stat.csv", sep = ",")
    #print("\nSelected Solution: ", solution)
    
   
        
    return(solution)



