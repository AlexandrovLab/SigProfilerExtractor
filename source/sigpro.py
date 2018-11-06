#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 11:04:39 2018

@author: mishugeb
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 13:39:29 2018

@author: mishugeb
"""
import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 
import sys
import argparse
import scipy.io
import numpy as np
import pandas as pd
import scipy.sparse as spr
import nimfa
import time
from sklearn import metrics
import pickle
from functools import partial
from numpy import linalg as LA
import time
import subroutines as sub
sys.path.append('../SigProfilerMatrixGenerator/scripts/')
import sigProfilerMatrixGeneratorFunc as datadump 
import sigProfilerPlotting as plot
import string   
import shutil
#Take the external inputs

parser = argparse.ArgumentParser(description='Extraction of mutational signatures from Cancer genomes')
    
parser.add_argument('-t','--input_type', type=str, 
                    help= 'The input type. There are three available input types: "vcf", "text", "matobj". The "vcf" type input will load the mutational catalogue from\
                    a varriant caller data. As a reminder the user has to create a project folder and place that to correct location (please see the readme file for the details\
                    about creating and location). The "text" input will load the matutational catalogue from a plain text file delimited by tab. The user has to place the text\
                    file in the input folder beforehand. Finally, the "matobj" type input will load the mutational catalogue from the matlab object file and the user has to\
                    place the matlab object file in the input folder beforehand as well')

parser.add_argument('-p','--project', type =str,
                    help= 'Name of the project file. This argument is mandatory for the "vcf" type input')

parser.add_argument('-r','--refgen', type =str,
                    help= 'Name of the reference genome. This argument is mandatory for the "vcf" type input')
parser.add_argument('-i','--inputfile', type=str, help='The name of the input file. This argument is mandatory for "text" or "matobj" type input' )

parser.add_argument('-s','--startprocess', type=int, 
                    help= 'The minimum number of processes to be extracted')

parser.add_argument('-e','--endprocess', type=int, 
                    help= 'The maximum number of processes to be extracted.')

parser.add_argument('-n','--n_iterations', type=int, 
                    help= 'The number of iterations to be executed. This parameter is valid only for the "vcf" input.')
 
parser.add_argument('-c','--cpu', type=int, 
                    help= 'The number of cpu to be executed for parallel computation. The default value will use the maximum number of the available cpus')
parser.add_argument ('-o','--output', type =str, 
                     help= 'Posisional argument. Name of the output directory.')
  
parser.add_argument('-l',"--layer", help='Optional parameter that set if the signatures will extrated in a hierarchical manner', action='store_true')
parser.add_argument("--exome", help='Optional parameter instructs script to create the catalogues using only the exome regions. Whole genome context by default. This parameter is valid only for the "vcf" input.', action='store_true')
parser.add_argument("--indel", help='Optional parameter instructs script to create the catalogue for limited INDELs. This parameter is valid only for the "vcf" input.', action='store_true')
parser.add_argument("--extended_indel", help='Optional parameter instructs script to create the catalogue for extended INDELs. This parameter is valid only for the "vcf" input.', action='store_true')

parser.add_argument('-m','--mtypes', type =str,
                    help= 'The context of mutations and is optional. This is valid when the input type is "vcf". User should pass the inteded mutation types among to be analyzed\
                    separeted by coma "," with no space. The sigporfiler engine will analyze the specific mutation\
                    types those are passed to this argument. The valid mutation type are  6, 12, 96, 1536, 192, 3072 and  DINUC.\
                    For example, if the user wants analyze mutation type 96, 192 and DINUC, that person should pass\
                    "--mtypes 96,192,DINUC" as in the argument. If the argument is not used, all the mutations will\
                    be analyzed')

args=parser.parse_args()

 


   

   
    
################################ take the inputs from the mandatory arguments ####################################
if args.input_type:
    input_type = args.input_type
else:
    raise Exception("Please provide an among ('vcf', 'text' or 'matobj'). Use --help to get more help.")
    

if args.output:
    out_put = args.output
else:
    raise Exception("Please provide the name of the output file. Use --help to get more help.")
    

   ################################ take the inputs from the general optional arguments ####################################
if args.startprocess:
    startProcess=args.startprocess
else:
    startProcess=1

if args.endprocess:
    endProcess=args.endprocess
else:
    endProcess=2

if args.n_iterations:
    totalIterations=args.n_iterations
else:
    totalIterations=3

if args.cpu:
    cpu = args.cpu
else:
    cpu = -1
    
if args.layer:
    hierarchi = True
else:
    hierarchi = False
   
if input_type=="text":
    ################################### For text input files ######################################################
    if args.inputfile:
        text_file = args.inputfile
        title = "" # set the title for plotting 
    else:
        raise Exception("Please provide an input file. Use --help to get more help.")
        
    data = pd.read_csv("../input/"+text_file, sep="\t").iloc[:,:-1]
    genomes = data.iloc[:,1:]
    genomes = np.array(genomes)
    
    
    #Contruct the indeces of the matrix
    #setting index and columns names of processAvg and exposureAvg
    index = data.iloc[:,0]
    colnames  = data.columns[1:]
    
    #creating list of mutational type to sync with the vcf type input
    mtypes = [str(genomes.shape[0])]
    os.chdir("../")  # get out of the source directory
###############################################################################################################

###########################################################################################################################################################################################
elif input_type=="matobj":
    ################################# For matlab input files #######################################################
    if args.inputfile:
        mat_file = args.inputfile
        title = "" # set the title for plotting 
    else:
        raise Exception("Please provide an input file. Use --help to get more help.")
        
    
    mat = scipy.io.loadmat('../input/'+ mat_file)
    mat = sub.extract_input(mat)
    genomes = mat[1]
    
    
  
    
    #Contruct the indeces of the matrix
    #setting index and columns names of processAvg and exposureAvg
    index1 = mat[3]
    index2 = mat[4]
    index = []
    for i, j in zip(index1, index2):
        index.append(i[0]+"["+j+"]"+i[2])
    colnames = np.array(pd.Series(mat[2]))
    index = np.array(pd.Series(index))
    
    #creating list of mutational type to sync with the vcf type input
    mtypes = [str(genomes.shape[0])]
    
    os.chdir("../")  # get out of the source directory
    #################################################################################################################
    
    
    
elif input_type=="vcf":
    ################################# For vcf input files #######################################################
    if args.project:
        project = args.project
        title = project # set the title for plotting 
    else:
        raise Exception("Please provide the project name. Use --help to get more help.")
       
    
    if args.refgen:
        refgen = args.refgen
    else:
        raise Exception("Please provide the refence genome name. Use --help to get more help.")
        
    
    if args.exome:
        exome = True
    else:
        exome = False

    if args.extended_indel:
        indel = True
    else: 
        indel = False

    if args.indel:
        indel = True
        limited_indel = True
        mlist = ["96", "DINUC", "INDEL"] 
    else:
        limited_indel = False
        mlist = ["96", "DINUC"] 
        
        
    os.chdir("../SigProfilerMatrixGenerator/scripts")
    #data = datadump.sigProfilerMatrixGeneratorFunc ("project2", "GRCh37") 
    data = datadump.sigProfilerMatrixGeneratorFunc (project, refgen, exome= exome, indel=limited_indel, indel_extended=indel, bed_file=None) 
    #create a data to broadcast
    
    
    
# =============================================================================
#     mlist = []
#     for m in data:
#         mlist.append(m)
# =============================================================================
        
     
    if args.mtypes:
        
        mtypes = args.mtypes.split(",")
        if any(x not in mlist for x in mtypes):
             raise Exception("Please pass valid mutation types seperated by comma with no space. Also please use the uppercase characters")
            
             
    else:
         mtypes = mlist
    #print (mtypes)
    #change working directory 
    
    os.chdir("../../")
    
       
###########################################################################################################################################################################################                  
for m in mtypes:
    
    # Determine the types of mutation which will be needed for exporting and copying the files
    if not (m=="DINUC"or m=="INDEL"):
        mutation_type = "SNV"
        
    else:
        mutation_type = m
    
    if input_type=="vcf":
        genomes = data[m]
        index = genomes.index.values
        colnames  = genomes.columns
        
        
        
    #create output directories to store all the results 
    output = out_put+"/"+m
    
    
   
    est_genomes = np.zeros([1,1])
    listofsignatures=[]
    H_iteration = 1 
    
    # While loop starts here
    while genomes.shape[1]>10 and est_genomes.shape[1]!=genomes.shape[1]:
        genomes = np.array(genomes)
        information =[] 
        if hierarchi is True:
            layer_directory = output+"/L"+str(H_iteration)
        elif hierarchi is False:
            layer_directory = output
            
        try:
            if not os.path.exists(layer_directory):
                os.makedirs(layer_directory)
                #os.makedirs(output+"/pickle_objects")
                #os.makedirs(output+"/All solutions")
            
         
        
    
            
        except: 
            print ("The {} folder could not be created".format("output"))
        
        
        fh = open(layer_directory+"/results_stat.csv", "w")   
        fh.write("Number of signature, Reconstruction Error, Process stability\n") 
        fh.close()
        # The following for loop operates to extract data from each number of signature
        for i in range(startProcess,endProcess+1):
                
            processAvg, \
            exposureAvg, \
            processStd, \
            exposureStd, \
            avgSilhouetteCoefficients, \
            clusterSilhouetteCoefficients, \
            finalgenomeErrors, \
            finalgenomesReconstructed, \
            finalWall, \
            finalHall, \
            processes = sub.decipher_signatures(genomes= genomes, \
                                                i = i, \
                                                totalIterations=totalIterations, \
                                                cpu=cpu, \
                                                mut_context=m) 
            
            
            ####################################################################### add sparsity in the exposureAvg #################################################################
            
            stic = time.time() 
            for s in range(exposureAvg.shape[1]):
                #print (s)
                exposureAvg[:,s] = sub.remove_all_single_signatures(processAvg, exposureAvg[:,s], genomes[:,s])
                #print ("Optimization for Sample {} is completed".format(s+1))
                #print ("\n\n\n\n")
            stoc = time.time()
            print ("Optimization time is {} seconds".format(stoc-stic))
                
            ##########################################################################################################################################################################
            # store the resutls of the loop            
            loopResults = [genomes, processAvg, exposureAvg, processStd, exposureStd, avgSilhouetteCoefficients, clusterSilhouetteCoefficients, finalgenomeErrors, finalgenomesReconstructed, finalWall, finalHall, processes]    
            information.append([processAvg, exposureAvg]) #Will be used during hierarchical approach
            
            ################################# Export the results ###########################################################    
            sub.export_information(loopResults, m, layer_directory, index, colnames)
            
            
            
        ################################################################################################################
        ########################################## Plot Stabiltity vs Reconstruction Error #############################        
        ################################################################################################################    
        
        solution = sub.stabVsRError(layer_directory+"/results_stat.csv", layer_directory, title)
        
        if os.path.exists(layer_directory+"/Selected solution"):
            shutil.rmtree(layer_directory+"/Selected solution") 
        # Copy the best solution the "selected solution" folder
        solutionFolderFrom= layer_directory+"/All solutions/"+str(solution)+" "+ mutation_type+ " Signature"
        solutionFolderTo = layer_directory+"/Selected solution/"+str(solution)+" "+ mutation_type+ " Signature"
        shutil.copytree(solutionFolderFrom, solutionFolderTo)
    
# =============================================================================
#         # get the solution for this hierarchy
#         for i in range(len(avgSilhouetteCoefficients)):
#             if avgSilhouetteCoefficients[-(i+1)]>0.85:
#                 solution = len(avgSilhouetteCoefficients)-(i)
#                 break
#             else:
#                 solution = 0
# =============================================================================
        if hierarchi is True:
            processAvg = information[solution-1][0]
            exposureAvg = information[solution-1][1]
            #del information
            
            est_genomes = np.dot(processAvg, exposureAvg) 
            
            low_similarity_idx = []
            for i in range(genomes.shape[1]):
                similarity = sub.cos_sim(genomes[:,i], est_genomes[:,i])
                #print (similarity)
                if similarity < 0.95:
                    low_similarity_idx.append(i)
            
            
            if len(low_similarity_idx)==0:   
                low_similarity_idx = []
            #print(low_similarity_idx)
            
            
            listofsignatures.append(processAvg) 
            genomes = genomes[:,low_similarity_idx]
            colnames=colnames[low_similarity_idx]
            H_iteration = H_iteration + 1
        
        elif hierarchi is False:
            break
        
        
        
    
       