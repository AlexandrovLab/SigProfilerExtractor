#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 12:21:06 2019

@author: mishugeb
"""

from SigProfilerExtractor import subroutines as sub
import numpy as np
import pandas as pd
import os




def decompose(signatures, activities, samples, output,  genome_build="GRCh37", verbose=False):

    
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
    exposureAvg = pd.read_csv(activities, sep = "\t", index_col = 0)
    
    
    
    processAvg = np.array(processAvg)
    
    
    genomes = pd.read_csv(samples, sep = "\t", index_col = 0)
    mutation_type = str(genomes.shape[0])
    #creating list of mutational type to sync with the vcf type input
    if mutation_type == "78":
        mutation_context = "DBS78"
    elif mutation_type == "83":
        mutation_context = "ID83"
    else:
        mutation_context = "SBS"+mutation_type
        
    signature_names = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype = mutation_context)
    exposureAvg.columns=signature_names   
    
    index = genomes.index
    m=mutation_type
    layer_directory2 = output
    if not os.path.exists(layer_directory2):
        os.makedirs(layer_directory2)
    
    
    if processAvg.shape[0]==1536: #collapse the 1596 context into 96 only for the deocmposition 
        processAvg = pd.DataFrame(processAvg, index=index)
        processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
        genomes = genomes.groupby(genomes.index.str[1:8]).sum()
        index = genomes.index
        processAvg = np.array(processAvg)
       
    
    final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build, mutation_context=mutation_context)    
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
                            remove_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg , penalty=0.01, background_sigs=background_sigs, verbose=verbose, genome_build=genome_build)

    return result


