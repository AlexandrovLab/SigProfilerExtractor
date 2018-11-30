#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 22:06:07 2018

@author: mishugeb
"""
import os
from mpi4py import MPI
import time 
import subroutines as sub
import numpy as np
import scipy.io
import pandas as pd
import pickle
from numpy import linalg as LA
import sys
import argparse
from random import shuffle
sys.path.append('../SigProfilerMatrixGenerator/scripts/')
import sigProfilerMatrixGeneratorFunc as datadump 
import string   
import shutil
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sigProfilerPlotting as plot
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


# take the inputs to master node and then broad cast them to the worker nodes
if rank == 0:

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
     
    
    parser.add_argument ('-o','--output', type =str, 
                         help= 'Posisional argument. Name of the output directory.')
      
    
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
    parser.add_argument("--snv",help='Optional parameter instructs script to create the catalogue for SNVs. This parameter is valid only for the "vcf" input.', action='store_true')
    
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
    
    
    
    
    
    if input_type=="text":
        ################################### For text input files ######################################################
        if args.inputfile:
            text_file = args.inputfile
            title = "" # set the title for plotting 
        else:
            raise Exception("Please provide an input file. Use --help to get more help.")
            
        data = pd.read_csv("../input/"+text_file, sep="\t").iloc[:,:]
        data=data.dropna(axis=1, inplace=False)
        genomes = data.iloc[:,1:]
        genomes = np.array(genomes)
        
        
        #Contruct the indeces of the matrix
        #setting index and columns names of processAvg and exposureAvg
        index = data.iloc[:,0]
        colnames  = data.columns[1:]
        
        
        #create a data to broadcast
        broadcast_data=genomes
        #creating list of mutational type to sync with the vcf type input
        mtypes = [str(genomes.shape[0])]
        os.chdir("../")  # get out of the source directory
        broadcast = [broadcast_data, input_type, mtypes, startProcess, endProcess]
        ###############################################################################################################
    
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
        
        
        #create a data to broadcast
        broadcast_data=genomes
        #Contruct the indeces of the matrix
        #setting index and columns names of processAvg and exposureAvg
        index1 = mat[3]
        index2 = mat[4]
        index = []
        for i, j in zip(index1, index2):
            index.append(i[0]+"["+j+"]"+i[2])
        colnames = pd.Series(mat[2])
        index = pd.Series(index)
        
        #creating list of mutational type to sync with the vcf type input
        mtypes = [str(genomes.shape[0])]
        
        os.chdir("../")  # get out of the source directory
        broadcast = [broadcast_data, input_type, mtypes, startProcess, endProcess]
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
        
        if args.snv:
            SNV = True
        else:
            SNV = False
            
        
            
        if args.exome:
            exome = True
        else:
            exome = False
    
        if args.extended_indel:
            indel = True
        else: 
            indel = False
    
        if args.indel and args.snv:
            indel = True
            limited_indel = True
            mlist = ["96", "DINUC", "INDEL"] 
            #print(mlist)
        elif args.snv:
            limited_indel = False
            mlist = ["96", "DINUC"] 
            #print(mlist)
        elif args.indel:
            indel = True
            limited_indel = True
            mlist = ["INDEL"] 
            #print(mlist)
        else:
            raise Exception("Please pass either --snv or --indel or both arguments to get a valide result")
            
        os.chdir("../SigProfilerMatrixGenerator/scripts")
        #data = datadump.sigProfilerMatrixGeneratorFunc ("project2", "GRCh37") 
        data = datadump.sigProfilerMatrixGeneratorFunc(project, refgen, exome=exome, SNVs=SNV,indel=limited_indel, indel_extended=indel, bed_file=None, chrom_based=False, plot=False, gs=False)
        
        #print (data)
        #create a data to broadcast
        
        broadcast_data=data
        
# =============================================================================
#         mlist = []
#         for m in data:
#             mlist.append(m)
# =============================================================================
            
         
      
            
        if args.mtypes:
            mkeys = data.keys()
            mtypes = args.mtypes.split(",")
            if any(x not in mkeys for x in mtypes):
                raise Exception("Please pass valid mutation types seperated by comma with no space. Also please use the uppercase characters")
                
                 
        else:
             mtypes = mlist
        #print (mtypes)
        #change working directory 
        os.chdir("../../")
        #print(os.getcwd())
        broadcast = [broadcast_data, input_type, mtypes, startProcess, endProcess]
        #print (broadcast)
        
        

else:
    broadcast = None


broadcast = comm.bcast(broadcast, root=0)






################################################# Start the main computations from the master node ##############################################

if rank == 0:
    for m in mtypes:
        
        if input_type=="vcf":
            genomes = data[m]
            #print (genomes)
            index = genomes.index.values
            colnames  = genomes.columns
            #print ("index is okay", index)
            #print ("colnames is okay", colnames)
            genomes = np.array(genomes)
    
        #create output directories to store all the results 
        output = out_put+"/"+m
        try:
            if not os.path.exists(output):
                    os.makedirs(output)
                    #os.makedirs(output+"/pickle_objects")
                    os.makedirs(output+"/All solutions")
                    #os.makedirs(output+"/Selected solution")
                    
            #Variables for the final output for all signatures in the csv file
            signatures = list() # will store the signature numbers
            norm = list() #will store the reconstruction errors from each number of signatures
            stb = list()  #will store the process stabilities for each number of signatures 
            fh = open(output+"/results_stat.csv", "w")   
            fh.write("Number of signature, Reconstruction Error, Process stability\n") 
            fh.close()
            
                    
        except: 
            print ("The {} folder could not be created".format("output"))
            
        
        
        # The following for loop operates to extract data from each number of signature
        for i in range(startProcess,endProcess+1):
            genomes = np.array(genomes) # Convert the data from dataframe to array to use the remove_single_signature function
            # The initial values accumute the results for each number of 
            totalMutationTypes = genomes.shape[0];
            totalGenomes = genomes.shape[1];
            totalProcesses = i
            Wall = np.zeros((totalMutationTypes, totalProcesses * totalIterations));
            #print (Wall.shape)
            Hall = np.zeros((totalProcesses * totalIterations, totalGenomes));
            genomeErrors = np.zeros((totalMutationTypes, totalGenomes, totalIterations));
            genomesReconstructed = np.zeros((totalMutationTypes, totalGenomes, totalIterations))
            print ("Extracting signature {} for mutation type {}".format(i, m))
            
            
            
    ##############################################################################################################################################################################         
    ############################################################# start of message passing interface for this loop ###############################################################  
    ##############################################################################################################################################################################         
            tic = time.time()
            results = []  #this list accumulate the data received from different workers
            iteration = 0
            process = 1
            firstBatch = 1
            
            # Send the first batch of processes to the nodes
            while firstBatch< min(size, totalIterations):
                comm.send(1, dest=process)
                #print ("Sending signal to process",process)
                iteration += 1
                process += 1
                firstBatch+=1
                
            
            
            # Wait for the data to come back
            received_processes = 0
            iteration_id = 1
            while received_processes< totalIterations :
                results.append(comm.recv(source=MPI.ANY_SOURCE, tag=1))
                process = comm.recv(source=MPI.ANY_SOURCE, tag=2)
                print ("Iteration {} completed for signature {}".format(iteration_id,  i))
                received_processes += 1
                iteration_id += 1
                
                if iteration<totalIterations :
                    comm.send(1, dest=process)
                    #print ("Sending on signal to process",process)
                    iteration += 1
                    
            # Send the shutdown signal
            for process in range(1,size):
                comm.send(0, dest=process, tag=1)
                    
                
            toc = time.time()
            print ("Time taken to collect {} iterations for {} signatures is {} seconds".format(totalIterations , i, round(toc-tic, 2)))
    ##############################################################################################################################################################################       
    ######################################################### end of message passing interface for this loop #####################################################################      
    ##############################################################################################################################################################################        
            
            avgSilhouetteCoefficients = -1.1
            clusterSilhouetteCoefficients = [0]
            processclust=[0]
            exposerclust=[0]
            finalWall=[0]
            finalHall = [0]
            finalgenomeErrors=[0]
            finalgenomesReconstructed = [0]
            
            for k in range(25):
                shuffle(results)
                Wall = np.zeros((totalMutationTypes, totalProcesses * totalIterations));
                #print (Wall.shape)
                Hall = np.zeros((totalProcesses * totalIterations, totalGenomes));
                genomeErrors = np.zeros((totalMutationTypes, totalGenomes, totalIterations));
                genomesReconstructed = np.zeros((totalMutationTypes, totalGenomes, totalIterations))
                
                processCount=0
                for j in range(len(results)):
                    W = results[j][0]
                    H = results[j][1]
                    genomeErrors[:, :, j] = genomes -  np.dot(W,H);
                    genomesReconstructed[:, :, j] = np.dot(W,H);
                    #print ("W", W.shape)
                    Wall[ :, processCount : (processCount + totalProcesses) ] = W;
                    Hall[ processCount : (processCount + totalProcesses), : ] = H;
                    processCount = processCount + totalProcesses;
                #print (Wall.shape, Hall.shape)
                
                
                W= np.array_split(Wall, totalIterations, axis=1)
                H= np.array_split(Hall, totalIterations, axis=0)
                   
                
                loop_processclust, loop_exposerclust, loop_avgSilhouetteCoefficients, loop_clusterSilhouetteCoefficients= sub.find_clusters_v1(W, H)
                
                #print ("stability", loop_avgSilhouetteCoefficients)
               
                if loop_avgSilhouetteCoefficients>avgSilhouetteCoefficients:
                    avgSilhouetteCoefficients=loop_avgSilhouetteCoefficients
                    clusterSilhouetteCoefficients = loop_clusterSilhouetteCoefficients
                    processclust = loop_processclust
                    exposerclust = loop_exposerclust
                    finalWall = Wall
                    finalHall = Hall
                    finalgenomeErrors = genomeErrors 
                    finalgenomesReconstructed = genomesReconstructed
                
            
            # computing the avg and std of the processes and exposures:
            processes=i
            processAvg = np.zeros((genomes.shape[0], processes))
            exposureAvg = np.zeros((processes, genomes.shape[1]))
            processStd = np.zeros((genomes.shape[0], processes))
            exposureStd = np.zeros((processes, genomes.shape[1]))
            
            for j in range(0, processes):
                processAvg[:,j]=np.mean(processclust[j], axis=1)
                processStd[:,j] = np.std(processclust[j], axis=1)
                exposureAvg[j,:] = np.mean(exposerclust[j], axis=1)
                exposureStd[j,:] = np.std(exposerclust[j], axis=1)
                
                
             ####################################################################### add sparsity in the exposureAvg #################################################################
            stic = time.time() 
            for s in range(exposureAvg.shape[1]):
                #print (i)
                exposureAvg[:,s] = sub.remove_all_single_signatures(processAvg, exposureAvg[:,s], genomes[:,s])
                #print ("\n\n\n\n")
            stoc = time.time()
            print ("Optimization time is {} seconds".format(stoc-stic))
        ##########################################################################################################################################################################    
                
            loopResults =[genomes, processAvg, exposureAvg, processStd, exposureStd, avgSilhouetteCoefficients, clusterSilhouetteCoefficients, genomeErrors, genomesReconstructed, Wall, Hall, processes]
           
            
           
            
    #################################################################### The result exporting part ################################################################    
    
            # Determine the mutation type:
            if not (m=="DINUC"or m=="INDEL"):
                mutation_type = "SNV"
                
            else:
                mutation_type = m
            # Create the neccessary directories
            subdirectory = output+"/All solutions/"+str(i)+" "+ mutation_type+ " Signature"
            if not os.path.exists(subdirectory):
                os.makedirs(subdirectory)
        
        
        
            #Export the loopResults as pickle objects
        
            resultname = "signature"+str(i)


    
            #Export the loopResults as pickle objects
            
# =============================================================================
#             resultname = "signature"+str(i)
#             
#             f = open(output+"/pickle_objects/"+resultname, 'wb')
#         
#             pickle.dump(loopResults, f)
#             f.close()
#             
# =============================================================================
       
            #preparing the column and row indeces for the Average processes and exposures:  
            listOfSignatures = []
            letters = list(string.ascii_uppercase)
            letters.extend([i+b for i in letters for b in letters])
            letters = letters[0:i]
            
            for j,l in zip(range(i),letters)  :
                listOfSignatures.append("Signature "+l)
            listOfSignatures = np.array(listOfSignatures)
            
        
                
            
            #Extract the genomes, processAVG, processStabityAvg
            genome= loopResults[0]
            #print ("genomes are ok", genome)
            processAvg= (loopResults[1])
            exposureAvg= (loopResults[2])
            processStabityAvg= (loopResults[5])
       
            
            # Calculating and listing the reconstruction error, process stability and signares to make a csv file at the end
            reconstruction_error = LA.norm(genome-np.dot(processAvg, exposureAvg), 'fro')/LA.norm(genome, 'fro')
            norm.append(reconstruction_error) 
            stb.append(processStabityAvg)
            signatures.append(loopResults[-1])
            
            
            # Preparing the results to export as textfiles for each signature
            
            #First exporting the Average of the processes
            processAvg= pd.DataFrame(processAvg)
            processes = processAvg.set_index(index)
            processes.columns = listOfSignatures
            processes = processes.rename_axis("signatures", axis="columns")
            #print("process are ok", processes)
            processes.to_csv(subdirectory+"/processes.txt", "\t", index_label=[processes.columns.name])  
            
            #Second exporting the Average of the exposures
            exposureAvg = pd.DataFrame(exposureAvg)
            exposures = exposureAvg.set_index(listOfSignatures)
            exposures.columns = colnames
            #print("exposures are ok", exposures)
            exposures.to_csv(subdirectory+"/exposures.txt", "\t", index_label=[exposures.columns.name]) 
            
            
           
            fh = open(output+"/results_stat.csv", "a") 
            print ('The reconstruction error is {} and the process stability is {} for {} signatures\n\n'.format(reconstruction_error, round(processStabityAvg,4), i))
            fh.write('{}, {}, {}\n'.format(i, reconstruction_error, processStabityAvg))
            fh.close()
            
            
                    
           ########################################### PLOT THE SIGNATURES ################################################
            if m=="DINUC" or m=="78":
                plot.plotDBS(subdirectory+"/processes.txt", subdirectory+"/Signature_plot" , "", "78", True)
            elif m=="INDEL" or m=="83":
                plot.plotID(subdirectory+"/processes.txt", subdirectory+"/Signature_plot" , "", "94", True)
            else:
                plot.plotSBS(subdirectory+"/processes.txt", subdirectory+"/Signature_plot", "", m, True)
                
            
            
            processAvg = pd.read_csv(subdirectory+"/processes.txt", sep="\t", index_col=0)
            exposureAvg = pd.read_csv(subdirectory+"/exposures.txt", sep="\t", index_col=0)
            probability = sub.probabilities(processAvg, exposureAvg)
            probability=probability.set_index("Sample")
            probability.to_csv(subdirectory+"/probabilities.txt", "\t") 
            
        ################################################################################################################
        ########################################## Plot Stabiltity vs Reconstruction Error #############################        
        ################################################################################################################    
        
        solution = sub.stabVsRError(output+"/results_stat.csv", output, title)
        
        if os.path.exists(output+"/Selected solution"):
            shutil.rmtree(output+"/Selected solution") 
        # Copy the best solution the "selected solution" folder
        solutionFolderFrom= output+"/All solutions/"+str(solution)+" "+ mutation_type+ " Signature"
        solutionFolderTo = output+"/Selected solution/"+str(solution)+" "+ mutation_type+ " Signature"
        shutil.copytree(solutionFolderFrom, solutionFolderTo)
        

else:
    
    data = broadcast[0]
    input_type= broadcast[1]
    mtypes = broadcast[2]
    startProcess = broadcast[3]
    endProcess = broadcast[4]
    
    
    
    for m in mtypes:
        
        if input_type=="vcf":
            genomes = data[m]
            genomes = np.array(genomes) # Convert the data from dataframe to array to use the remove_single_signature function
        else:
            genomes = data
            genomes = np.array(genomes) # Convert the data from dataframe to array to use the remove_single_signature function
        for i in range(startProcess,endProcess+1):
            while True:
                start = comm.recv(source=0)
                #print ("The signal is", start)
                if start == 0:
                    break
                W, H= sub.nmf(genomes=genomes, totalProcesses =i)
                comm.send([W,H], dest=0, tag=1)
                comm.send(rank, dest=0, tag=2)
                #print ("sending data from rank {} for an iteration of signature {}".format(rank, i))
                #comm.barrier()
  
