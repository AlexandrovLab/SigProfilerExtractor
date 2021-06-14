#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 12:21:06 2019

@author: mishugeb
"""

from SigProfilerExtractor import subroutines as sub
import numpy as np
import pandas as pd
import SigProfilerExtractor as cosmic
import os




def decompose(signatures, activities, samples,  output, signature_database=None, nnls_add_penalty=0.05, 
              nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, 
              genome_build="GRCh37", refit_denovo_signatures=True, make_decomposition_plots=True, 
              collapse_to_SBS96=True,connected_sigs=True, verbose=False):

    
    """
    Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples.
    
    Parameters: 
        
        signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs. 
        activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
        samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
        output: A string. Path to the output folder.
        genome_build = A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
        verbose = Boolean. Prints statements. Default value is False. 
        
    Values:
        The files below will be generated in the output folder. 
        
        Cluster_of_Samples.txt
        comparison_with_global_ID_signatures.csv
        Decomposed_Solution_Activities.txt
        Decomposed_Solution_Samples_stats.txt
        Decomposed_Solution_Signatures.txt
        decomposition_logfile.txt
        dendogram.pdf
        Mutation_Probabilities.txt
        Signature_assaignment_logfile.txt
        Signature_plot[MutatutionContext]_plots_Decomposed_Solution.pdf
        
    Example:
        >>>from SigProfilerExtractor import decomposition as decomp
        >>>signatures = "path/to/dDe_Novo_Solution_Signatures.txt"
        >>>activities="path/to/De_Novo_Solution_Activities.txt"
        >>>samples="path/to/Samples.txt"
        >>>output="name or path/to/output.txt"
        decomp.decompose(signatures, activities, samples, output, genome_build="GRCh37", verbose=False)

    """
    
    processAvg = pd.read_csv(signatures, sep = "\t", index_col=0) 
    originalProcessAvg=processAvg
    exposureAvg = pd.read_csv(activities, sep = "\t", index_col = 0)
    genomes = pd.read_csv(samples, sep = "\t", index_col = 0)
    
    mutation_type = str(genomes.shape[0])
    m=mutation_type
    
    index = genomes.index
    colnames = genomes.columns
    listOfSignatures = processAvg.columns
    
    
    
    
    
    
    
    
    #creating list of mutational type to sync with the vcf type input
    if mutation_type == "78":
        mutation_context = "DBS78"
    elif mutation_type == "83":
        mutation_context = "ID83"
    elif mutation_type=="48":
        mutation_context = "CNV48"
        print("Mutation Type is: CNV")
        #paths = cosmic.__path__[0]
        #sigDatabase = pd.read_csv(paths+"/data/CNV_signatures.txt", sep="\t",index_col=0)
        #genomes=genomes.loc[sigDatabase.index]
        #processAvg=processAvg.loc[sigDatabase.index]
    else:
        mutation_context = "SBS"+mutation_type
        
    processAvg = np.array(processAvg)    
    signature_names = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype = mutation_context)
    exposureAvg.columns=signature_names   
    
    
    
    
    
    # create the folder for the final solution/ De Novo Solution
    
    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except: 
        print ("The {} folder could not be created".format("output"))
    
    # make the texts for signature plotting
    layer_directory1 = output+"/De_Novo_Solution"
    try:
        if not os.path.exists(layer_directory1):
            os.makedirs(layer_directory1)
    except: 
        print ("The {} folder could not be created".format("De_Novo_Solution"))
    
    
   
    
    
   
    listOfSignatures = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype=mutation_context)
    genomes = pd.DataFrame(genomes)
    
    
   
    denovo_exposureAvg = np.array(exposureAvg.T)
    exposureAvg = sub.make_final_solution(processAvg, genomes, listOfSignatures, layer_directory1, m, index,\
                   colnames,denovo_exposureAvg  = denovo_exposureAvg, add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, initial_remove_penalty=initial_remove_penalty, de_novo_fit_penalty=de_novo_fit_penalty, connected_sigs=connected_sigs, refit_denovo_signatures=refit_denovo_signatures)    

    
    
    layer_directory2 = output+"/Decompose_Solution"
    try:
        if not os.path.exists(layer_directory2):
            os.makedirs(layer_directory2)
    except: 
        print ("The {} folder could not be created".format("Decomposed_Solution"))
    
    
    if processAvg.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
        processAvg = pd.DataFrame(processAvg, index=index)
        processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
        genomes = genomes.groupby(genomes.index.str[1:8]).sum()
        index = genomes.index
        processAvg = np.array(processAvg)
    
    if processAvg.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
        processAvg = pd.DataFrame(processAvg, index=index)
        processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
        genomes = pd.DataFrame(genomes, index=index)
        genomes = genomes.groupby(genomes.index.str[2:9]).sum()
        index = genomes.index
        processAvg = np.array(processAvg)
            
    
    final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build,signature_database=signature_database, mutation_context=mutation_context, add_penalty=0.05, connected_sigs=connected_sigs,remove_penalty=0.01, make_decomposition_plots=make_decomposition_plots, originalProcessAvg=originalProcessAvg)    
    #final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build)
    # extract the global signatures and new signatures from the final_signatures dictionary
    globalsigs = final_signatures["globalsigs"]
    globalsigs = np.array(globalsigs)
    newsigs = final_signatures["newsigs"]
    processAvg = np.hstack([globalsigs, newsigs])  
    allsigids = final_signatures["globalsigids"]+final_signatures["newsigids"]
    attribution = final_signatures["dictionary"]
    background_sigs= final_signatures["background_sigs"]
    index = genomes.index
    colnames = genomes.columns
    
    
    
    
    result = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory2, m, index, colnames, \
                            cosmic_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg ,  background_sigs=background_sigs, verbose=verbose, genome_build=genome_build, add_penalty=nnls_add_penalty, remove_penalty=nnls_remove_penalty, initial_remove_penalty=initial_remove_penalty,connected_sigs=connected_sigs,refit_denovo_signatures=False)

    return result


