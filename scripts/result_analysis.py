#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 11:50:55 2018

@author: mishugeb
"""
import numpy as np
from numpy import linalg as LA
import pickle
import argparse
#import saved results from hard disk



parser = argparse.ArgumentParser(description='Extraction results from sigprofiler analysis')
    
parser.add_argument('results', type =str,
                    help= 'Name of the file containing the results from sigprofiler analysis')


args = parser.parse_args()


f = open('output/'+args.results, 'rb')
result = pickle.load(f)
f.close()

signatures = list()
norm = list()
stb = list()


for i in result:
    genome= i[0]
    processAvg= (i[1])
    exposureAvg= (i[2])
    processStabityAvg= (i[5])
    
    norm.append(LA.norm(genome-np.dot(processAvg, exposureAvg), 'fro'))
    stb.append(processStabityAvg)
    signatures.append(i[-1])
 
    
fh = open("output/results_stat.csv", "w")   
fh.write("Number of signature, Reconstruction Error, Process stability\n") 
for i, j, k in zip(signatures, norm, stb):
    print ('The reconstruction error is {} and the process stability is {} for {} signatures'.format(j, k, i))
    fh.write('{}, {}, {}\n'.format(i, j, k))
    


fh.close()