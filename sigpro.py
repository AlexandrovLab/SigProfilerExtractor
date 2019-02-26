#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 13:39:29 2018

@author: S M Ashiqul Islam (Mishu)


    ##########################################
    SigProfilerExtractor (``sigproextractor``)
    ##########################################
    
    SigProfilerExtractor allows de novo extraction of mutational signatures from data 
    generated in a matrix format. The tool identifies the number of operative mutational 
    signatures, their activities in each sample, and the probability for each signature to 
    cause a specific mutation type in a cancer sample. The tool makes use of SigProfilerMatrixGenerator 
    and SigProfilerPlotting. 

"""
import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import scipy.io
import numpy as np
import pandas as pd
import time
from sigproextractor import subroutines as sub
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as datadump   
import shutil
import multiprocessing as mp
import sigproextractor as cosmic



def importdata(datatype="matobj"):
    
    """
    Imports the path of example data.
    
    parameters
    ----------
    
    datatype: A string. Type of data. The type of data should be one of the following:
            - "vcf": used for vcf format data.
            - "text": used for text format data.
            - "matobj": used for matlab object format data.
    
    
    Returns:
    -------

    The path of the example data.

    Example: 
    -------
    >>> from sigproextractor import sigpro as sig
    >>> data = sig.importdata("text")
    
    This "data" variable can be used as a parameter of the "project" argument of the sigProfilerExtractor function
        
    """
    
    paths = cosmic.__path__[0]
    if datatype=="matobj":
        data = paths+"/data/21_breast_WGS_substitutions.mat"
    elif datatype=="text":
        data = paths+"/data/all_mice_silvio.txt"
    elif datatype=="csv":
        data = paths+"/data/csvexample.csv"
    elif datatype=="vcf":
        directory = os.getcwd()
        dataold = paths+"/data/vcftest"
        datanew = directory+"/vcftest"
        if not os.path.exists(datanew):
            shutil.copytree(dataold , datanew)
        
        data="vcftest"
    return data


def sigProfilerExtractor(input_type, out_put, project, refgen="GRCh37", startProcess=1, endProcess=10, totalIterations=8, cpu=-1, hierarchy = False, mtype = ["default"],exome = False): 
    
    """
    Extracts mutational signatures from an array of samples.
    
    
    Parameters
    ----------
    
    input_type: A string. Type of input. The type of input should be one of the following:
            - "vcf": used for vcf format inputs.
            - "text": used for text format inputs.
            - "matobj": used for matlab object format of inputs. 
        
    out_put: A string. The name of the output folder. The output folder will be generated in the current working directory. 
            
    project: A string. Name of the input folder (in case of "vcf" type input) or the input file (in case of "text" or "matobj" type input). The project file or folder should be inside the current working directory. For the "vcf" type input,the project has to be a folder which will contain the vcf files in vcf format or text formats. The "text" or "matobj" type projects have to be a file. "matobj" projects should have .mat extension.  
            
    refgen: A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".
            
    startProcess: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1 
    
    endProcess: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 10
    
    totalIterations: A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 8
            
    cpu: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors. 
    
    hierarchy: boolean, optional. Defines if the signature will be extracted in a hierarchical fashion. The default value is "False".
    
    mtype: A list of strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is ["96", "DINUC" , "INDEL"].
            
    exome: boolean, optional. Defines if the exomes will be extracted. The default value is "False".
    
    
    
    
    Returns
    -------
    
    After sigProfilerExtractor is successfully executed, an output directory will be generated in the current working directory 
    according to the name of the  parameter of the "out_put" argument. In the "output" directory there will be subfolder 
    for each type of mutational contexts. 
    
    If the "hierarchy" parameter is false, inside of each mutational context subdirectory, there will be subdirectories named 
    "All solutions" and "Final solution". Besides the subdirectories, there will be a file named "results_stat.csv" which 
    will contain the record of the relative reconstruction error and process stability for each number of signatures. 
    Another file named stibility.pdf will contain the plot of recontruction error vs process stability. The "All solution"
    directory will contain the subdirectories for each number of signatures which will further contain the solution files 
    ("signature.txt", "exposure.txt", "probabilities.txt" and a pdf file that depicts the  proportion of the mututaions 
    for each number signatures. On the other hand, the "Final solution" directory contains two subdirectories: "De Novo Solution"
    and "Decomposed Solution". The "De Novo Solution" subdirectory will contain the solution files for the optimum number of 
    "De Novo Signatures" signatures with a dendrogram file where the samples are clustered by the de novo signatures. The "Decomposed 
    Solution" subfolder contains the records  where "De Novo Signatures" are further decomposed into the global signatures. 
    
    If the "hierarchy" parameter is true, inside of each mutational context subdirectory, there will be a subdirectory named
    "All_Solution_by_Layer" which will further contain the solutions  in the layer (L) subdirectories. Everything else will be similar to
    the previously deccribed directory structures. The structure of the result folder is synopsized below:
        
        If Hierarchy is False:
            
        -Mutational Context folder
            -All solution folder
                -Signature folder
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
            -Selected_Solution folder
                -De_Novo_Solution folder
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
                    -dendrogram plot file
                -Decomposed_Solution folder
                    -comparison with global signature.csv file
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
                    -dendrogram plot file
            -results_stat.csv file
            -stability plot pdf
            
                    
        If Hierarchy is True:
            
        -Mutational Context folder
            -All Solution by Layer folder
                -Layer folder (L)
                    -All solution folder
                        -Signature folder
                            -exposure.txt file
                            -signature.txt file
                            -probabilities.txt file
                            -signature plot pdf file
                    -L1_solution folder
                        -exposure.txt file
                        -signature.txt file
                        -probabilities.txt file
                        -signature plot pdf file
                    -results_stat.csv file
                    -stability plot pdf
            -Selected_Solution folder
                -De_Novo_Solution folder
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
                    -dendrogram plot file
                -Decomposed_Solution folder
                    -comparison with global signature.csv file
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
                    -dendrogram plot file
            -results_stat.csv file
            -stability plot pdf
    
    Examples
    --------
    
    >>> from sigproextractor import sigpro as sig
    >>> data = sig.importdata("vcf")
    >>> sig.sigProfilerExtractor("vcf", "example_output", data, startProcess=1, endProcess=3)
    
    Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
    Check the current working directory for the "example_output" folder.
    
    
    """
    ################################ take the inputs from the mandatory arguments ####################################
    input_type = input_type;
    out_put = out_put;  
    project = project
        
    
    ################################ take the inputs from the general optional arguments ####################################
    startProcess=startProcess ; 
    endProcess=endProcess;
    totalIterations=totalIterations
    cpu = cpu
    hierarchi = hierarchy
    
       
    if input_type=="text":
        
        ################################### For text input files ######################################################
        
        text_file = project
        title = "" # set the title for plotting 
            
    
            
        data = pd.read_csv(text_file, sep="\t").iloc[:,:]
        data=data.dropna(axis=1, inplace=False)
        data = data.loc[:, (data != 0).any(axis=0)]
        genomes = data.iloc[:,1:]
        genomes = np.array(genomes)
        allgenomes = genomes.copy()  # save the allgenomes for the final results 
        
        #Contruct the indeces of the matrix
        #setting index and columns names of processAvg and exposureAvg
        index = data.iloc[:,0]
        colnames  = data.columns[1:]
        allcolnames = colnames.copy() # save the allcolnames for the final results
        
        #creating list of mutational type to sync with the vcf type input
        mtypes = [str(genomes.shape[0])]
        if mtypes[0] == "78":
            mtypes = ["DBS78"]
        elif mtypes[0] == "94":
            mtypes = ["ID83"]
        
    ###############################################################################################################
    
        
    ###########################################################################################################################################################################################
    elif input_type=="csv":
    ################################# For matlab input files #######################################################
    
        filename = project
        title = "" # set the title for plotting 
        
        genomes, index, colnames, mtypes = sub.read_csv(filename)   
        allgenomes = genomes.copy()
        allcolnames = colnames.copy() 
        
        
        # Define the mtypes
        mtypes = [mtypes]
        if mtypes[0] == "78":
            mtypes = ["DBS78"]
        elif mtypes[0] == "94":
            mtypes = ["ID83"]
        
       
        
    
    #################################################################################################################
    
    
        ###########################################################################################################################################################################################
    elif input_type=="matobj":
        ################################# For matlab input files #######################################################
        
        mat_file = project
        title = "" # set the title for plotting 
        
            
        
        mat = scipy.io.loadmat(mat_file)
        mat = sub.extract_input(mat)
        genomes = mat[1]
        allgenomes = genomes.copy()  # save the allgenomes for the final results 
        
        
      
        
        #Contruct the indeces of the matrix
        #setting index and columns names of processAvg and exposureAvg
        index1 = mat[3]
        index2 = mat[4]
        index = []
        for i, j in zip(index1, index2):
            index.append(i[0]+"["+j+"]"+i[2])
        colnames = np.array(pd.Series(mat[2]))
        allcolnames = colnames.copy() # save the allcolnames for the final results
        index = np.array(pd.Series(index))
        
        #creating list of mutational type to sync with the vcf type input
        mtypes = [str(genomes.shape[0])]
        if mtypes[0] == "78":
            mtypes = ["DBS78"]
        elif mtypes[0] == "94":
            mtypes = ["ID83"]
        
        #################################################################################################################
        
        
    elif input_type=="vcf":
        ################################# For vcf input files #######################################################
        
        project = project
        title = project # set the title for plotting 
        
        refgen = refgen
        
        
        exome = exome
    
        
        
    
            
        project_name = project.split("/")[-1]
        data = datadump.SigProfilerMatrixGeneratorFunc(project_name, refgen, project, exome=exome,  bed_file=None, chrom_based=False, plot=False, gs=False)
        
        
        
    
        # Selecting the mutation types    
        if mtype != ["default"]:
            mkeys = data.keys()
            mtypes = mtype
            if any(x not in mkeys for x in mtypes):
                 raise Exception("Please pass valid mutation types seperated by comma with no space. Carefully check (using SigProfilerMatrixGenerator)"\
                                 "what mutation contexts should be generated by your VCF files. Also please use the uppercase characters")
                
                 
        else:
            if set(["96", "DINUC", "INDEL"]).issubset(data):            
                mtypes = ["96", "DINUC", "INDEL"] 
            elif set(["96", "DINUC"]).issubset(data): 
                mtypes = ["96", "DINUC"]
            elif set(["INDEL"]).issubset(data):            
                mtypes = ["INDEL"] 
        #print (mtypes)
        #change working directory 
        
        
    else:
        raise ValueError("Please provide a correct input_type. Check help for more details")
          
    ###########################################################################################################################################################################################                  
    for m in mtypes:
        
        # Determine the types of mutation which will be needed for exporting and copying the files
        if not (m=="DINUC"or m=="INDEL"):
            mutation_type = "SBS"+m
            
        else:
            if m == "DINUC":
                mutation_type = "DBS78"
            elif m== "INDEL":
                mutation_type = "ID83"
                
        
        if input_type=="vcf":
            genomes = data[m]
            genomes = genomes.loc[:, (genomes != 0).any(axis=0)]
            allgenomes = genomes.copy()  # save the allgenomes for the final results 
            index = genomes.index.values
            colnames  = genomes.columns
            allcolnames = colnames.copy() # save the allcolnames for the final results 
            
           
            
        #create output directories to store all the results 
        output = out_put+"/"+mutation_type
        
        
        
        
        est_genomes = np.zeros([1,1])
        listofsignatures=[]
        listofsignaturesSTE = []
        H_iteration = 1 
        flag = True # We need to enter into the first while loop regardless any condition
        # While loop starts here
        while flag:
            genomes = np.array(genomes)
            

            information =[] 
            if hierarchi is True:
                layer_directory = output+"/All_Solution_Layer/L"+str(H_iteration)
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
             
            similarity_dataframe = pd.DataFrame({"Samples": list(colnames)})
            
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
                
    
                # remove signatures only if the process stability is above a thresh-hold of 0.85
                if  avgSilhouetteCoefficients> -0.85:   
                    stic = time.time() 
                    pool = mp.Pool()
                    results = [pool.apply_async(sub.remove_all_single_signatures_pool, args=(x,processAvg,exposureAvg,genomes,)) for x in range(genomes.shape[1])]
                    pooloutput = [p.get() for p in results]
                    #print(results)
                    pool.close()
                    
                    for i in range(len(pooloutput)):
                        #print(results[i])
                        exposureAvg[:,i]=pooloutput[i]
                    stoc = time.time()
                    print ("Optimization time is {} seconds".format(stoc-stic))    
                
                
                ##########################################################################################################################################################################
                # store the resutls of the loop.  Here,  processStd and exposureStd are standard Errors, NOT STANDARD DEVIATIONS.           
                loopResults = [genomes, processAvg, exposureAvg, processStd, exposureStd, avgSilhouetteCoefficients, clusterSilhouetteCoefficients, finalgenomeErrors, finalgenomesReconstructed, finalWall, finalHall, processes]    
                information.append([processAvg, exposureAvg, processStd, exposureStd]) #Will be used during hierarchical approach
                
                ################################# Export the results ###########################################################    
                sub.export_information(loopResults, m, layer_directory, index, colnames)
                
                # Compute the estimated genome from the processAvg and exposureAvg
                est_genomes = np.dot(processAvg, exposureAvg) 
                
                #check the similarities between the original and estimated genome for each number of signatures 
                all_similarities = []
                for i in range(genomes.shape[1]):
                    similarity = sub.cos_sim(genomes[:,i], est_genomes[:,i])
                    similarity = round(similarity,2)
                    #print(similarity)
                    all_similarities.append(similarity)
                similarity_dataframe["Signature "+str(processes)] = all_similarities
            ################################################################################################################
            ########################################## Plot Stabiltity vs Reconstruction Error #############################        
            ################################################################################################################    
            # Print the Stabiltity vs Reconstruction Error as get the solution as well
            solution = sub.stabVsRError(layer_directory+"/results_stat.csv", layer_directory, title)
            #print ("The mutution type is %s"%(m)
            
            
            
            
            
            ################################### Hierarchical Extraction  #########################
            if hierarchi is True:
                
                # write the name of Samples and Matrix participating in each Layer.
                layer_genome = pd.DataFrame(genomes)
                layer_genome = layer_genome.set_index(index)
                layer_genome.columns = colnames
                layer_genome.to_csv(layer_directory+"/Samples_in_Layer_"+str(H_iteration)+".text", sep = "\t")
                similarity_dataframe.to_csv(layer_directory+"/Samples_Similarity_in_Layer_"+str(H_iteration)+".text", sep = "\t")
                del layer_genome
                
                sample_record = open(output+"/Samples_Selected_by_Layers.text", "a")
                sample_record.write("\nSamples participating in Layer"+str(H_iteration)+"\n"+"Total number of samples in this layer is: "+str(len(colnames))+"\n\n" )
                for sn in colnames:
                    # sn is the abbreviation of "Sample Name", used as a iterator variable
                    sample_record.write(sn+" ,\n" )
                sample_record.write("######################################################################################\n")   
                sample_record.write("######################################################################################\n")   
                sample_record.write("######################################################################################\n")   
                sample_record.write("######################################################################################\n\n\n\n\n")                        
                sample_record.close()
                    
                
                
                if os.path.exists(layer_directory+"/L"+str(H_iteration)+"_solution"):
                    shutil.rmtree(layer_directory+"/L"+str(H_iteration)+"_solution") 
                # Copy the best solution the "selected solution" folder
                solutionFolderFrom= layer_directory+"/All_solutions/"+mutation_type+"_Signature_"+str(solution)
                solutionFolderTo = layer_directory+"/L"+str(H_iteration)+"_Solution/"+mutation_type+"_Signature_"+str(solution)
                shutil.copytree(solutionFolderFrom, solutionFolderTo)
                
                # load the best processAvg, exposureAvg and processSTE based on the solution
                processAvg = information[solution-startProcess][0]
                exposureAvg = information[solution-startProcess][1]
                processSTE = information[solution-startProcess][2]
                #del information
                
                # Compute the estimated genome from the processAvg and exposureAvg
                est_genomes = np.dot(processAvg, exposureAvg) 
                
                # make the list of the samples which have similarity lower than the thresh-hold with the estimated ones
                low_similarity_idx = []
                
                for i in range(genomes.shape[1]):
                    similarity = sub.cos_sim(genomes[:,i], est_genomes[:,i])
                    
                    # The tresh-hold for hierarchy is 0.95 for now
                    if similarity < 0.90:    
                        low_similarity_idx.append(i)
                
                
                if len(low_similarity_idx)==0:   
                    low_similarity_idx = []
                #print(low_similarity_idx)
                
                # Accumulated the signatures and signaturesSTE for the final results
                listofsignatures.append(processAvg) 
                listofsignaturesSTE.append(processSTE)
                
                
                genomes = genomes[:,low_similarity_idx]
                colnames=colnames[low_similarity_idx]
                H_iteration = H_iteration + 1
                
                
                
                #########################################################################################################
                # do the necessary operations and put the outputs in the "Final Solution" folder when the while loop ends
                if genomes.shape[1]<10 or est_genomes.shape[1]==genomes.shape[1]:
                    flag = False #update the flag for the whileloop
                    
                    # create the folder for the final solution/ De Novo Solution
                    layer_directory1 = output+"/Selected_Solution/De_Novo_Solution"
                    try:
                        if not os.path.exists(layer_directory1):
                            os.makedirs(layer_directory1)
                    except: 
                        print ("The {} folder could not be created".format("output"))
            
                    
                    
                    count = 0
                    for p,q in zip(listofsignatures, listofsignaturesSTE):
                        if count==0:
                            processAvg=p
                            processSTE = q
                        else:
                            processAvg = np.hstack([processAvg, p]) 
                            processSTE = np.hstack([processSTE, q])
                        count+=1
                        
                    
                    # make de novo solution(processAvg, allgenomes, layer_directory1)
                    listOfSignatures = sub.make_letter_ids(idlenth = processAvg.shape[1])
                    sub.make_final_solution(processAvg, allgenomes, listOfSignatures, layer_directory1, m, index, allcolnames, process_std_error = processSTE)    
                    
                    try:
                        # create the folder for the final solution/ Decomposed Solution
                        layer_directory2 = output+"/Selected_Solution/Decomposed_Solution"
                        try:
                            if not os.path.exists(layer_directory2):
                                os.makedirs(layer_directory2)
                        except: 
                            print ("The {} folder could not be created".format("output"))
                
                        
                        final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2)                
                        # extract the global signatures and new signatures from the final_signatures dictionary
                        globalsigs = final_signatures["globalsigs"]
                        globalsigs = np.array(globalsigs)
                        newsigs = final_signatures["newsigs"]
                        processAvg = np.hstack([globalsigs, newsigs])  
                        allsigids = final_signatures["globalsigids"]+final_signatures["newsigids"]
                        
                        sub.make_final_solution(processAvg, allgenomes, allsigids, layer_directory2, m, index, allcolnames)
                    except:
                        print("\nWARNING!!! We apolozize we don't have a global signature database for the mutational context you provided. We have a database only for SBS96, DINUC and INDELS.\nTherefore no result for signature Decomposition is generated." )
                    
                
                    
                #######################################################################################################
            elif hierarchi is False:
                
                # write the name of Samples and Matrix participating in the experiment.
                layer_genome = pd.DataFrame(genomes)
                layer_genome = layer_genome.set_index(index)
                layer_genome.columns = colnames
                layer_genome.to_csv(output+"/Samples.text", sep = "\t")
                similarity_dataframe.to_csv(output+"/Samples_Similarity.text", sep = "\t")
                del layer_genome
                
                ################################### Decompose the new signatures into global signatures   #########################
                processAvg = information[solution-startProcess][0]
                processSTE = information[solution-startProcess][2]
                
               
                # create the folder for the final solution/ De Novo Solution
                layer_directory1 = output+"/Selected_Solution/De_Novo_Solution"
                try:
                    if not os.path.exists(layer_directory1):
                        os.makedirs(layer_directory1)
                except: 
                    print ("The {} folder could not be created".format("output"))
                
                # make de novo solution(processAvg, allgenomes, layer_directory1)
                listOfSignatures = sub.make_letter_ids(idlenth = processAvg.shape[1])
                sub.make_final_solution(processAvg, allgenomes, listOfSignatures, layer_directory1, m, index, allcolnames)    
               
                try:
                   # create the folder for the final solution/ Decomposed Solution
                    layer_directory2 = output+"/Selected_Solution/Decomposed_Solution"
                    try:
                        if not os.path.exists(layer_directory2):
                            os.makedirs(layer_directory2)
                    except: 
                        print ("The {} folder could not be created".format("output"))
                
                    
                    final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2)
                    # extract the global signatures and new signatures from the final_signatures dictionary
                    globalsigs = final_signatures["globalsigs"]
                    globalsigs = np.array(globalsigs)
                    newsigs = final_signatures["newsigs"]
                    processAvg = np.hstack([globalsigs, newsigs])  
                    allsigids = final_signatures["globalsigids"]+final_signatures["newsigids"]
                    
                    sub.make_final_solution(processAvg, genomes, allsigids, layer_directory2, m, index, colnames, process_std_error = processSTE)
                
                except:
                    print("\nWARNING!!! We apolozize we don't have a global signature database for the mutational context you provided. We have a database only for SBS96, DINUC and INDELS.\nTherefore no result for signature Decomposition is generated." )
                    shutil.rmtree(layer_directory2)
                
                
               
                break
    print("\n\n \nYour Job Is Successfully Terminated! Thank You For Using SigProfiler Extractor.\n ")
             
                

     

if __name__=="__main__":
    
    sigProfilerExtractor("text", "textfunc", "all_mice_silvio.txt", refgen="GRCh37", startProcess=1, endProcess=2, totalIterations=3, \
                         cpu=-1, hierarchy = False, mtype = ["default"],exome = False, indel_extended = False)
