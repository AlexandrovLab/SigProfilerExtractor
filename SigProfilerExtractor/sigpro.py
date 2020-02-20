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
import scipy
import scipy.io
import sklearn
import numpy as np
import pandas as pd
import time
from SigProfilerExtractor import subroutines as sub
import SigProfilerMatrixGenerator
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as datadump   
import shutil
import multiprocessing as mp
import SigProfilerExtractor as cosmic
import platform
import datetime
import psutil
import sigProfilerPlotting 
from SigProfilerExtractor import single_sample as ss
def memory_usage():
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse1 = py.memory_info()[0]/2.**30  # memory use in GB...I think
    print('\n************** Reported Current Memory Use: '+ str(round(memoryUse1,2))+" GB *****************\n")
    #print('\n************** Reported Current Memory Use: '+ str(round(memoryUse2,2))+" GB *****************\n")






def importdata(datatype="matrix"):
    
    """
    Imports the path of example data.
    
    parameters
    ----------
    
    datatype: A string. Type of data. The type of data should be one of the following:
            - "vcf": used for vcf format data.
            - "matrix": used for text format data. This format represents the catalog of mutations seperated by tab. 
            - "matobj": used for matlab object format data.
            
    
    
    Returns:
    -------

    The path of the example data.

    Example: 
    -------
    >>> from SigProfilerExtractor import sigpro as sig
    >>> data = sig.importdata("table")
    
    This "data" variable can be used as a parameter of the "project" argument of the sigProfilerExtractor function
        
    """
    
    paths = cosmic.__path__[0]
    if datatype=="matobj":
        data = paths+"/data/21_breast_WGS_substitutions.mat"
    elif datatype=="text" or datatype=="table" or datatype=="matrix":
        data = paths+"/data/Samples.txt"
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


def sigProfilerExtractor(input_type, out_put, input_data, refgen="GRCh37", genome_build = 'GRCh37', startProcess=1, endProcess=10, totalIterations=100, init="alexandrov-lab-custom", cpu=-1,  mtype = "default",exome = False, penalty=0.05, resample = True, wall= False, gpu=False): 
    memory_usage()
    """
    Extracts mutational signatures from an array of samples.
    
    
    Parameters
    ----------
    
    input_type: A string. Type of input. The type of input should be one of the following:
            - "vcf": used for vcf format inputs.
            - "matrix": used for table format inputs using a tab seperated file.
             
        
    out_put: A string. The name of the output folder. The output folder will be generated in the current working directory. 
            
    input_data: A string. Name of the input folder (in case of "vcf" type input) or the input file (in case of "table"  type input). The project file or folder should be inside the current working directory. For the "vcf" type input,the project has to be a folder which will contain the vcf files in vcf format or text formats. The "text"type projects have to be a file.   
            
    refgen: A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".
            
    startProcess: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1 
    
    endProcess: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 10
    
    totalIterations: A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 100
    
    init: A String. The initialization algorithm for W and H matrix of NMF
    
    wall: A Boolean. If true, the Ws and Hs from all the NMF iterations are generated in the output. 
            
    cpu: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors. 
    
    mtype: A list of strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is ["96", "DINUC" , "ID"], where "96" is the SBS96 context, "DINUC"
    is the DINULEOTIDE context and ID is INDEL context. 
            
    exome: Boolean, optional. Defines if the exomes will be extracted. The default value is "False".
    
    penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the thresh-hold cutoff to asaign signatures to a sample.    
    
    resample: Boolean, optional. Default is True. If True, add poisson noise to samples by resampling.  
    
    
    Returns
    -------
    To learn about the output, please visit https://osf.io/t6j7u/wiki/home/
    
    
    Examples
    --------
    
    >>> from SigProfilerExtractor import sigpro as sig
    >>> data = sig.importdata("vcf")
    >>> sig.sigProfilerExtractor("vcf", "example_output", data, startProcess=1, endProcess=3)
    
    Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
    Check the results in the "example_output" folder.
    """
    
    if gpu == True:
        import torch
    
        if gpu and (torch.cuda.device_count() == 0):
            raise RuntimeError("GPU not available!")
    
    
    #################################### At first create the system data file ####################################
    if not os.path.exists(out_put):
        os.makedirs(out_put)
    sysdata = open(out_put+"/JOB_METADATA.txt", "w")
    sysdata.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
    sysdata.write("-------System Info-------\n")
    sysdata.write("Operating System Name: "+ platform.uname()[0]+"\n"+"Nodename: "+platform.uname()[1]+"\n"+"Release: "+platform.uname()[2]+"\n"+"Version: "+platform.uname()[3]+"\n")
    sysdata.write("\n-------Python and Package Versions------- \n")
    sysdata.write("Python Version: "+str(platform.sys.version_info.major)+"."+str(platform.sys.version_info.minor)+"."+str(platform.sys.version_info.micro)+"\n")
    sysdata.write("Sigproextractor Version: "+cosmic.__version__+"\n")
    sysdata.write("SigprofilerPlotting Version: "+sigProfilerPlotting.__version__+"\n")
    sysdata.write("SigprofilerMatrixGenerator Version: "+SigProfilerMatrixGenerator.__version__+"\n")
    sysdata.write("Pandas version: "+pd.__version__+"\n")
    sysdata.write("Numpy version: "+np.__version__+"\n")
    sysdata.write("Scipy version: "+scipy.__version__+"\n")
    sysdata.write("Scikit-learn version: "+sklearn.__version__+"\n")
    #sysdata.write("Nimfa version: "+nimfa.__version__+"\n")
    
    
    
    sysdata.write("\n-------Vital Parameters Used for the execution -------\n")
    #format the project_name first:
    project = input_data  #will use this variable as the parameter for project argument in SigprofilerMatrixGenerator
    try:
        if project[-1] != "/":
            project_name = project.split("/")[-1]   #will use this variable as the parameter for project_name argument in SigprofilerMatrixGenerator
        else:
            project_name = project.split("/")[-2]
    except:
        project_name = "Input from DataFrame"
    sysdata.write("input_type: {}\ninputdata: {}\nstartProcess: {}\nendProcess: {}\ntotalIterations: {}\ncpu: {}\nrefgen: {}\ngenome_build: {}\nmtype: {} \ninit: {}\n".format(input_type, project_name, startProcess, endProcess, totalIterations, cpu, refgen, genome_build, mtype, init))
    
    sysdata.write("\n-------Date and Time Data------- \n")
    tic = datetime.datetime.now()
    sysdata.write("Date and Clock time when the execution started: "+str(tic)+"\n")
    sysdata.close()
    
    
    
    
    
    ################################ take the inputs from the mandatory arguments ####################################
    input_type = input_type;
    out_put = out_put;  
    #project = input_data   #the variable was already set above
        
    
    ################################ take the inputs from the general optional arguments ####################################
    startProcess=startProcess ; 
    endProcess=endProcess;
    
    totalIterations=totalIterations
    cpu = cpu
    hierarchy = False #No use
    
    
   
    
    
    
    if input_type=="text" or input_type =="table" or input_type=="matrix":
        
        ################################### For text input files ######################################################
        
        text_file = project
        title = "" # set the title for plotting 
            
    
        if type(text_file)!=str:
            data=text_file
        else:
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
        elif mtypes[0] == "83":
            mtypes = ["ID83"]
        else:
            mtypes = ["SBS"+mtypes[0]]
        
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
        mtypes = [str(genomes.shape[0])]
        if mtypes[0] == "78":
            mtypes = ["DINUC"]
        elif mtypes[0] == "83":
            mtypes = ["ID"]
        
       
        
    
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
            mtypes = ["DINUC"]
        elif mtypes[0] == "83":
            mtypes = ["ID"]
        
        #################################################################################################################
        
        
    elif input_type=="vcf":
        ################################# For vcf input files #######################################################
        
        project = project
        title = project # set the title for plotting 
        
        refgen = refgen
        
        
        exome = exome
    
        
        
    
            
        #project_name = project.split("/")[-1]
        data = datadump.SigProfilerMatrixGeneratorFunc(project_name, refgen, project, exome=exome,  bed_file=None, chrom_based=False, plot=False, gs=False)
        
        
        
        
        # Selecting the mutation types    
        if mtype == ["default"]:
            if set(["96", "DINUC", "ID"]).issubset(data):  
                mtypes = ["SBS96", "DBS78", "ID83"] 
            elif set(["96", "DINUC"]).issubset(data): 
                mtypes = ["SBS96", "DBS78"]
            elif set(["ID"]).issubset(data):            
                mtypes = ["ID83"]    
        
        elif mtype == "default":
            if set(["96", "DINUC", "ID"]).issubset(data):  
                mtypes = ["SBS96", "DBS78", "ID83"] 
            elif set(["96", "DINUC"]).issubset(data): 
                mtypes = ["SBS96", "DBS78"]
            elif set(["ID"]).issubset(data):            
                mtypes = ["ID83"]
          
        else:
            #mkeys = data.keys()
            mtype = mtype.upper()
            mtype = mtype.replace(" ", "")
            mtypes = mtype.split(",")
# =============================================================================
#             if any(x not in mkeys for x in mtypes):
#                  raise Exception("Please pass valid mutation types seperated by comma with no space. Carefully check (using SigProfilerMatrixGenerator)"\
#                                  "what mutation contexts should be generated by your VCF files. Also please use the uppercase characters")
# =============================================================================
            
             
            
        
        #change working directory 
        
        #set the genome_build
        genome_build=refgen 
        
    else:
        raise ValueError("Please provide a correct input_type. Check help for more details")
        
          
    ###########################################################################################################################################################################################                  
   
    for m in mtypes:
        
       
        mutation_context = m
        
        # we need to rename the m because users input could be SBS96, SBS1536, DBS78, ID83 etc
        if m.startswith("SBS"):
            m = m[3:] #removing "SBS"
        elif m.startswith("DBS"):
            m = "DINUC"
        elif m.startswith("ID"):
            m = "ID"
        
        # Determine the types of mutation which will be needed for exporting and copying the files
        if not (m=="DINUC" or m.startswith("DBS") or m.startswith("ID")):
            
            if m.startswith("SBS"):
                mutation_type = m
            else:
                mutation_type = "SBS"+m
            
        else:
            if m == "DINUC" or m.startswith("DBS"):
                mutation_type = "DBS78"
            elif m== "ID" or m.stratswith("ID"):
                mutation_type = "ID83"
                
       
            
        if input_type=="vcf":
            
            try: 
                
                genomes = pd.DataFrame(data[m])
            except: 
                raise Exception("Please pass valid mutation types seperated by comma with no space. Carefully check (using SigProfilerMatrixGenerator)"\
                                 "what mutation contexts should be generated by your VCF files. Also please use the uppercase characters")
            
            #check if the genome is a nonzero matrix
            shape= genomes.shape
            if shape==(0,0):
                sysdata = open(out_put+"/JOB_METADATA.txt", "a")
                sysdata.write("Sample is not a nonzero matrix for the mutation context "+ m+"\n")
                print("Sample is not a nozero matrix for the mutation context "+ m)
                sysdata.close()
                continue
                
            genomes = genomes.loc[:, (genomes != 0).any(axis=0)]
            
            allgenomes = genomes.copy()  # save the allgenomes for the final results 
            index = genomes.index.values
            colnames  = genomes.columns
            allcolnames = colnames.copy() # save the allcolnames for the final results 
            
        #check if start and end processes are bigger than the number of samples
        startProcess = min(startProcess, genomes.shape[1])
        endProcess = min(endProcess, genomes.shape[1])   
        
        #in the plotting funciton "ID" is used as "INDEL"
        if m=="ID":
            m="INDEL" #for plotting     
            
        #create output directories to store all the results 
        output = out_put+"/"+mutation_type
        
        est_genomes = np.zeros([1,1])
        H_iteration = 1 
        genomes = np.array(genomes)
        information =[] 
        layer_directory = output
        try:
            if not os.path.exists(layer_directory):
                os.makedirs(layer_directory)
                #os.makedirs(output+"/pickle_objects")
                #os.makedirs(output+"/All solutions")
        except: 
            print ("The {} folder could not be created".format("output"))
        
        
        fh = open(layer_directory+"/All_solutions_stat.csv", "w")   
        fh.write("Total Signatures,Stability,Matrix Frobenius%,avgStability\n") 
        fh.close()
        # The following for loop operates to extract data from each number of signature
        
        all_similirities_list = [] #this list is going to store the dataframes of different similirieties as items
        minimum_stabilities = []
        #similarity_dataframe = pd.DataFrame({"Sample Name": list(colnames)})
        
        
        
        # set up the seeds generation same matrices for different number of signatures
        seeds = np.random.randint(0, 10000000, size=totalIterations) # set the seeds ranging from 0 to 10000000 for resampling and same seeds are used in different number of signatures
        
        # get the cutoff for normatization to handle the hypermutators 
        
        normalization_cutoff = sub.get_normalization_cutoff(genomes)
        #print("Normalization Cutoff is :", normalization_cutoff)
        
        #genomes = sub.normalize_samples(genomes, normalize=False, all_samples=False, number=30000)
        
        
        for i in range(startProcess,endProcess+1):
            current_time_start = datetime.datetime.now()
            
            #memory_usage()    
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
            converge_information, \
            reconstruction_error, \
            processes = sub.decipher_signatures(genomes= genomes, \
                                                i = i, \
                                                totalIterations=totalIterations, \
                                                cpu=cpu, \
                                                mut_context=m, \
                                                resample = resample,
                                                seeds=seeds, 
                                                init = init,
                                                normalization_cutoff=normalization_cutoff,
                                                gpu=gpu, batch_size=128)
            
            
            
            
            
            
            
            
            #denormalize the genomes and exposures
            #genomes = sub.denormalize_samples(genomes, totalMutations, normalization_value=100000)
            #exposureStd = sub.denormalize_samples(exposureStd, totalMutations, normalization_value=100000)    
            ####################################################################### add sparsity in the exposureAvg #################################################################
            
            
            # remove signatures only if the process stability is above a thresh-hold of 0.85
            if  avgSilhouetteCoefficients> -1.0:   
                stic = time.time() 
                
                #removing signatures:
# =============================================================================
#                     pool = mp.Pool()
#                     results = [pool.apply_async(sub.remove_all_single_signatures_pool, args=(x,processAvg,exposureAvg,genomes,)) for x in range(genomes.shape[1])]
#                     pooloutput = [p.get() for p in results]
#                     
#                     #print(results)
#                     pool.close()
#                     
#                     for i in range(len(pooloutput)):
#                         #print(results[i])
#                         exposureAvg[:,i]=pooloutput[i]
# =============================================================================
                    
                #refitting signatures:
                #removing signatures:
                pool = mp.Pool()
                results = [pool.apply_async(ss.fit_signatures_pool, args=(genomes,processAvg,x,)) for x in range(genomes.shape[1])]
                pooloutput = [p.get() for p in results]
                pool.close()
                                    
                for i in range(len(pooloutput)):
                    
                    exposureAvg[:,i]=pooloutput[i][0] 
                    
                stoc = time.time()
                print ("Optimization time is {} seconds".format(stoc-stic))    
                
            #report progress to the system file:
            current_time_end = datetime.datetime.now()
            sysdata = open(out_put+"/JOB_METADATA.txt", "a")
            if  hierarchy is True:
                sysdata.write("\nSignature extraction for {} completed for layer {} {} signatures for {}! TimeStamp: {}\n".format(mutation_type,  H_iteration, processes,  current_time_end-current_time_start, current_time_end))
            else:
                sysdata.write("\nSignature extraction for {} completed for {} signatures for {}! TimeStamp: {}\n".format(mutation_type,  processes,  current_time_end-current_time_start, current_time_end))
            
            #Get total mutationation for each signature in reverse order and order the signatures from high to low mutation barden
            signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
            sorted_idx = np.argsort(-signature_total_mutations)
            processAvg = np.take(processAvg, sorted_idx, axis=1)
            exposureAvg = np.take(exposureAvg, sorted_idx, axis=0)
            signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
            
            signature_stats = pd.DataFrame({"Stability": clusterSilhouetteCoefficients, "Total Mutations": signature_total_mutations})
            minimum_stabilities.append(round(np.mean(clusterSilhouetteCoefficients),2)) #here minimum stability is the average stability !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Compute the estimated genome from the processAvg and exposureAvg
            est_genomes = np.dot(processAvg, exposureAvg) 
            
            #check the similarities between the original and estimated genome for each number of signatures
            
            all_similarities, cosine_similarities = sub.calculate_similarities(genomes, est_genomes, colnames)
            #print(totalMutations)
            ##########################################################################################################################################################################
            # store the resutls of the loop.  Here,  processStd and exposureStd are standard Errors, NOT STANDARD DEVIATIONS.           
            loopResults = [genomes, processAvg, exposureAvg, processStd, exposureStd, avgSilhouetteCoefficients, clusterSilhouetteCoefficients, signature_total_mutations, all_similarities, signature_stats, reconstruction_error, finalgenomeErrors, finalgenomesReconstructed, converge_information, finalWall, finalHall,  processes]    
            information.append([processAvg, exposureAvg, processStd, exposureStd, clusterSilhouetteCoefficients, signature_total_mutations, signature_stats, all_similarities]) #Will be used during hierarchycal approach
            
            ################################# Export the results ###########################################################    
            sub.export_information(loopResults, m, layer_directory, index, colnames, wall=wall)
            
          
            
            all_similirities_list.append(all_similarities)
                #
            #similarity_dataframe["Total Signatures "+str(processes)] = cosine_similarities
            
            
            
            
        ################################################################################################################
        ########################################## Plot Stabiltity vs Reconstruction Error #############################        
        ################################################################################################################    
        # Print the Stabiltity vs Reconstruction Error as get the solution as well
        solution, all_stats = sub.stabVsRError(layer_directory+"/All_solutions_stat.csv", layer_directory, title, all_similirities_list, mutation_type)
        all_stats.insert(0, 'Stability (Avg Silhouette)', minimum_stabilities) #!!!!!!!!!!!!!!!!1 here minimum stability is avg stability
        all_stats.to_csv(layer_directory+"/All_solutions_stat.csv", sep = ",")
        
        # add more information to results_stat.csv
         
        
        #Set index for the  the Similarity Dataframe
        #similarity_dataframe = similarity_dataframe.set_index("Sample Name")
        
        #Add the total mutations of each sample
        #sample_total_mutations = list(np.sum(genomes, axis =0))
       
        #similarity_dataframe.insert(loc=0, column = "Total Mutations", value = sample_total_mutations)
        
        
        
        # write the name of Samples and Matrix participating in each Layer.
        layer_genome = pd.DataFrame(genomes)
        layer_genome = layer_genome.set_index(index)
        layer_genome.columns = colnames
        layer_genome = layer_genome.rename_axis("Mutation Types", axis="columns")
        
        
        
        
        
# =============================================================================
#                 data_stat_folder = output+"/Data_Stats"
#                 try:
#                     if not os.path.exists(data_stat_folder):
#                         os.makedirs(data_stat_folder)
#                 except: 
#                         print ("The {} folder could not be created".format("Data_Stats"))
#                 
#                 layer_genome.to_csv(data_stat_folder+"/Samples.text", sep = "\t", index_label=[layer_genome.columns.name])
#                 similarity_dataframe.to_csv(data_stat_folder+"/Similatiry_Data_All_Sigs.text", sep = "\t")
#                 del layer_genome
#                 for i in range(startProcess,endProcess+1):
#                     all_similirities_list[i-startProcess].to_csv(data_stat_folder+"/Similatiry_Data_Sig_"+str(i)+".text", sep="\t")
# =============================================================================
        # record the samples
        layer_genome.to_csv(output+"/Samples.txt", sep = "\t", index_label=[layer_genome.columns.name])
        #similarity_dataframe.to_csv(data_stat_folder+"/Similatiry_Data_All_Sigs"+str(H_iteration)+".text", sep = "\t")
        del layer_genome
        ################################### Decompose the new signatures into global signatures   #########################
        processAvg = information[solution-startProcess][0]
        processSTE = information[solution-startProcess][2]
        signature_stabilities = information[solution-startProcess][4]
        signature_total_mutations = information[solution-startProcess][5]  
        signature_stats = information[solution-startProcess][6] 
        all_similarities = information[solution-startProcess][7]
        
       
        # create the folder for the final solution/ De Novo Solution
        layer_directory1 = output+"/Suggested_Solution/De_Novo_Solution"
        try:
            if not os.path.exists(layer_directory1):
                os.makedirs(layer_directory1)
        except: 
            print ("The {} folder could not be created".format("output"))
        
        # make the texts for signature plotting
        signature_stabilities = sub.signature_plotting_text(signature_stabilities, "Stability", "float")
        signature_total_mutations = sub.signature_plotting_text(signature_total_mutations, "Total Mutations", "integer")
        # make de novo solution(processAvg, allgenomes, layer_directory1)
        
        
       
        
        
       
        listOfSignatures = sub.make_letter_ids(idlenth = processAvg.shape[1], mtype=mutation_context)
        allgenomes = pd.DataFrame(allgenomes)
        
        
        
        
        
        exposureAvg = sub.make_final_solution(processAvg, allgenomes, listOfSignatures, layer_directory1, m, index, \
                       allcolnames, process_std_error = processSTE, signature_stabilities = signature_stabilities, \
                       signature_total_mutations = signature_total_mutations, signature_stats = signature_stats, penalty=penalty)    
          
        try:
            # create the folder for the final solution/ Decomposed Solution
            
            layer_directory2 = output+"/Suggested_Solution/Decomposed_Solution"
            try:
                if not os.path.exists(layer_directory2):
                    os.makedirs(layer_directory2)
            except: 
                print ("The {} folder could not be created".format("output"))
        
            if processAvg.shape[0]==1536: #collapse the 1596 context into 96 only for the deocmposition 
                processAvg = pd.DataFrame(processAvg, index=index)
                processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
                genomes = pd.DataFrame(genomes, index=index)
                genomes = genomes.groupby(genomes.index.str[1:8]).sum()
                index = genomes.index
                processAvg = np.array(processAvg)
                genomes = np.array(genomes)
                
           
            final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build, mutation_context=mutation_context)
            
            # extract the global signatures and new signatures from the final_signatures dictionary
            globalsigs = final_signatures["globalsigs"]
            globalsigs = np.array(globalsigs)
            newsigs = final_signatures["newsigs"]
            processAvg = np.hstack([globalsigs, newsigs])  
            allsigids = final_signatures["globalsigids"]+final_signatures["newsigids"]
            attribution = final_signatures["dictionary"]
            background_sigs= final_signatures["background_sigs"]
            genomes = pd.DataFrame(genomes)
            
            
            
            exposureAvg = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory2, m, index, colnames, \
                                    remove_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg , background_sigs=background_sigs, penalty=penalty, genome_build=genome_build)
            
        except:
            print("\nWARNING!!! We apolozize we don't have a global signature database for the mutational context you provided. We have a database only for SBS96, DINUC and INDELS.\nTherefore no result for signature Decomposition is generated." )
            shutil.rmtree(layer_directory2)
                
            
           
    
    sysdata = open(out_put+"/JOB_METADATA.txt", "a")
    toc = datetime.datetime.now()
    sysdata.write("\nDate and Clock time when the execution ended: "+str(toc)+"\n")
    
    sysdata.write("-------Job Status------- \n")
    sysdata.write("Analysis of mutational signatures completed successfully! Total execution time: "+str(toc-tic)+". Results can be found in: ["+out_put+"] folder")
    sysdata.close()

    print("\n\n \nYour Job Is Successfully Completed! Thank You For Using SigProfiler Extractor.\n ")
             


     


