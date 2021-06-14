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
import multiprocessing
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

    
def record_parameters(sysdata, excecution_parameters, start_time):
            #genomes = sub.normalize_samples(genomes, normalize=False, all_samples=False, number=30000)
            sysdata.write("\n--------------EXECUTION PARAMETERS--------------\n")
            
            sysdata.write("INPUT DATA\n")
            sysdata.write("\tinput_type: {}\n".format(excecution_parameters["input_type"]))
            sysdata.write("\toutput: {}\n".format(excecution_parameters["output"]))
            sysdata.write("\tinput_data: {}\n".format(excecution_parameters["input_data"]))
            sysdata.write("\treference_genome: {}\n".format(excecution_parameters["reference_genome"]))
            sysdata.write("\tcontext_types: {}\n".format(excecution_parameters["context_type"]))
            sysdata.write("\texome: {}\n".format(excecution_parameters["exome"]))
            
            sysdata.write("NMF REPLICATES\n")
            sysdata.write("\tminimum_signatures: {}\n".format(excecution_parameters["minimum_signatures"]))
            sysdata.write("\tmaximum_signatures: {}\n".format(excecution_parameters["maximum_signatures"]))
            sysdata.write("\tNMF_replicates: {}\n".format(excecution_parameters["NMF_replicates"]))
            
            sysdata.write("NMF ENGINE\n")
            sysdata.write("\tNMF_init: {}\n".format(excecution_parameters["NMF_init"]))
            sysdata.write("\tprecision: {}\n".format(excecution_parameters["precision"]))
            sysdata.write("\tmatrix_normalization: {}\n".format(excecution_parameters["matrix_normalization"]))
            sysdata.write("\tresample: {}\n".format(excecution_parameters["resample"]))
            sysdata.write("\tseeds: {}\n".format(excecution_parameters["seeds"]))
            sysdata.write("\tmin_NMF_iterations: {}\n".format(format(excecution_parameters["min_NMF_iterations"],',d')))
            sysdata.write("\tmax_NMF_iterations: {}\n".format(format(excecution_parameters["max_NMF_iterations"], ',d')))
            sysdata.write("\tNMF_test_conv: {}\n".format(format(excecution_parameters["NMF_test_conv"],',d')))
            sysdata.write("\tNMF_tolerance: {}\n".format(excecution_parameters["NMF_tolerance"]))
            
            sysdata.write("CLUSTERING\n")
            sysdata.write("\tclustering_distance: {}\n".format(excecution_parameters["dist"]))
            
            
            sysdata.write("EXECUTION\n")
            if excecution_parameters["cpu"]==-1:
                sysdata.write("\tcpu: {}; Maximum number of CPU is {}\n".format(multiprocessing.cpu_count(), multiprocessing.cpu_count()))
            else:
                sysdata.write("\tcpu: {}; Maximum number of CPU is {}\n".format(excecution_parameters["cpu"], multiprocessing.cpu_count()))
            sysdata.write("\tgpu: {}\n".format(excecution_parameters["gpu"]))
            
            sysdata.write("Solution Estimation\n")
            sysdata.write("\tstability: {}\n".format(excecution_parameters["stability"]))
            sysdata.write("\tmin_stability: {}\n".format(excecution_parameters["min_stability"]))
            sysdata.write("\tcombined_stability: {}\n".format(excecution_parameters["combined_stability"]))
            
            sysdata.write("COSMIC MATCH\n")
            sysdata.write("\topportunity_genome: {}\n".format(excecution_parameters["opportunity_genome"]))
            sysdata.write("\cosmic_version: {}\n".format(excecution_parameters["cosmic_version"]))
            sysdata.write("\tnnls_add_penalty: {}\n".format(excecution_parameters["nnls_add_penalty"]))
            sysdata.write("\tnnls_remove_penalty: {}\n".format(excecution_parameters["nnls_remove_penalty"]))
            sysdata.write("\tinitial_remove_penalty: {}\n".format(excecution_parameters["initial_remove_penalty"]))
            sysdata.write("\tde_novo_fit_penalty: {}\n".format(excecution_parameters["de_novo_fit_penalty"])) 
            sysdata.write("\trefit_denovo_signatures: {}\n".format(excecution_parameters["refit_denovo_signatures"]))
            sysdata.write("\tcollapse_to_SBS96: {}\n".format(excecution_parameters["collapse_to_SBS96"]))
            
            sysdata.write("\n-------Analysis Progress------- \n")
            sysdata.write("[{}] Analysis started: \n".format(str(start_time).split(".")[0]))
            

def sigProfilerExtractor(input_type, 
                         output, 
                         input_data, 
                         reference_genome="GRCh37", 
                         opportunity_genome = "GRCh37", 
                         cosmic_version=3.1,
                         context_type = "default", 
                         exome = False, 
                         minimum_signatures=1,
                         maximum_signatures=25,  
                         nmf_replicates=500, 
                         resample = True, 
                         batch_size=1, 
                         cpu=-1, 
                         gpu=False, 
                         nmf_init="random", 
                         precision= "single", 
                         matrix_normalization= "gmm", 
                         seeds= "random", 
                         min_nmf_iterations= 10000, 
                         max_nmf_iterations=1000000, 
                         nmf_test_conv= 10000, 
                         nmf_tolerance= 1e-15, 
                         nnls_add_penalty=0.05, 
                         nnls_remove_penalty=0.01,
                         de_novo_fit_penalty=0.02,
                         initial_remove_penalty=0.05,
                         refit_denovo_signatures=True,
                         collapse_to_SBS96=True,
                         clustering_distance="cosine",
                         export_probabilities=True,
                         make_decomposition_plots=True,
                         stability=0.8, 
                         min_stability=0.2, 
                         combined_stability=1.0,
                         get_all_signature_matrices= False): 
    memory_usage()
    """
    Extracts mutational signatures from an array of samples.
    
    
    Parameters
    ----------
    
    INPUT DATA:-
    
    input_type: A string. Type of input. The type of input should be one of the following:
            - "vcf": used for vcf format inputs.
            - "matrix": used for table format inputs using a tab seperated file.
             
        
    output: A string. The name of the output folder. The output folder will be generated in the current working directory. 
            
    input_data: A string. Name of the input folder (in case of "vcf" type input) or the input file (in case of "table"  type input). The project file or folder should be inside the current working directory. For the "vcf" type input,the project has to be a folder which will contain the vcf files in vcf format or text formats. The "text"type projects have to be a file.   
            
    reference_genome: A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".
       
    opportunity_genome: The build or version of the reference signatures for the reference genome. The default opportunity genome is GRCh37. If the input_type is "vcf", the genome_build automatically matches the input reference genome value.    
     
    context_type: A list of strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "SBS96,DBS78,ID83". 
    
    exome: Boolean, optional. Defines if the exomes will be extracted. The default value is "False".
    
    
    NMF RUNS:-
    
    minimum_signature: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1 
    
    maximum_signatures: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 10
    
    nmf_replicates: A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 100
    
    resample: Boolean, optional. Default is True. If True, add poisson noise to samples by resampling.  
    
    seeds: Boolean. Default is "random". If random, then the seeds for resampling will be random for different analysis.
                  If not random, then seeds will be obtained from a given path of a .txt file that contains a list of seed. 
    
    NMF RUNS:-
    
    matrix_normalization: A string. Method of normalizing the genome matrix before it is analyzed by NMF. Default is "log2". Other options are "gmm", "100X" or "no_normalization".         
    
    nmf_init: A String. The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'
              Default is 'nndsvd_min'.
    
    precision: A string. Values should be single or double. Default is single.
    
    min_nmf_iterations: An integer. Value defines the minimum number of iterations to be completed before NMF converges. Default is 2000.
    
    max_nmf_iterations: An integer. Value defines the maximum number of iterations to be completed before NMF converges. Default is 200000
    
    nmf_test_conv: An integer. Value definer the number number of iterations to done between checking next convergence.
            
    nmf_tolerance: A float. Value defines the tolerance to achieve to converge. 
    
    
    EXECUTION:-
    
    cpu: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available        processors. 
    
    gpu:Boolean, optional. Defines if the GPU resource will used if available. Default is False. If True, the GPU resource 
        will be used in the computation.

    batch_size: An integer. Will be effective only if the GPU is used. Defines the number of NMF replicates to be performed
              by each CPU during the parallel processing. Default is 1.
              
    
    SOLUTION ESTIMATION THRESH-HOLDS:-

    stability: A float. Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered. 

    min_stability: A float. Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered. 

    combined_stability: A float. Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered.            
    
    
    DECOMPOSITION:-
    
    de_novo_fit_penalty: Float, optional. Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to be assigned denovo signatures to a sample. 
    
    nnls_add_penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the strong (add) thresh-hold cutoff to be assigned COSMIC signatures to a sample. 
    
    nnls_remove_penalty: Float, optional. Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample.
     
    initial_remove_penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to be COSMIC assigned signatures to a sample.
    
    refit_denovo_signatures: Boolean, optional. Default is False. If True, then refit the denovo signatures with nnls.
    
    make_decomposition_plots: Boolean, optional. Defualt is True. If True, Denovo to Cosmic sigantures decompostion plots will be created as a part the results.

    
    OTHERS:-
    
    get_all_signature_matrices: A Boolean. If true, the Ws and Hs from all the NMF iterations are generated in the output.
    
    export_probabilities: A Boolean. Defualt is True. If False, then doesn't create the probability matrix.
    

    
    Returns
    -------
    To learn about the output, please visit https://osf.io/t6j7u/wiki/home/
    
    
    Examples
    --------
    
    Examples
    --------

    >>> from SigProfilerExtractor import sigpro as sig
    
    # to get input from vcf files
    >>> path_to_example_folder_containing_vcf_files = sig.importdata("vcf")
    >>> data = path_to_example_folder_containing_vcf_files # you can put the path to your folder containing the vcf samples
    >>> sig.sigProfilerExtractor("vcf", "example_output", data, minimum_signatures=1, maximum_signatures=3)
    
    Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
    Check the current working directory for the "example_output" folder.
    
    # to get input from table format (mutation catalog matrix)
    >>> path_to_example_table = sig.importdata("matrix")
    >>> data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog         matrix/table
    >>> sig.sigProfilerExtractor("matrix", "example_output", data, opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=3)
    
    Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
    Check the results in the "example_output" folder.
    """
    #record the start time
    start_time = datetime.datetime.now()
    
    #set the output variable
    out_put = output; 
    
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
    
    
    
    
    #format the project_name first:
    project = input_data  #will use this variable as the parameter for project argument in SigprofilerMatrixGenerator
    try:
        if project[-1] != "/":
            project_name = project.split("/")[-1]   #will use this variable as the parameter for project_name argument in SigprofilerMatrixGenerator
        else:
            project_name = project.split("/")[-2]
    except:
        project_name = "Input from DataFrame"
    
   
        
        
    excecution_parameters= {"input_type":input_type, 
                        "output":output, 
                        "input_data":input_data, 
                        "reference_genome":reference_genome, 
                        "opportunity_genome":opportunity_genome, 
                        "cosmic_version":cosmic_version,
                        "context_type":context_type,
                        "exome":exome,
                        "minimum_signatures":minimum_signatures, 
                        "maximum_signatures":maximum_signatures, 
                        "NMF_replicates":nmf_replicates, 
                        "cpu":cpu, 
                        "gpu":gpu, 
                        "batch_size":batch_size, 
                        "NMF_init":nmf_init,
                        "precision":precision,
                        "matrix_normalization":matrix_normalization,
                        "resample":resample, 
                        "seeds":seeds,
                        "min_NMF_iterations":min_nmf_iterations,
                        "max_NMF_iterations":max_nmf_iterations,
                        "NMF_test_conv": nmf_test_conv,
                        "NMF_tolerance": nmf_tolerance,
                        "nnls_add_penalty":nnls_add_penalty,
                        "nnls_remove_penalty":nnls_remove_penalty,
                        "initial_remove_penalty":initial_remove_penalty,
                        "de_novo_fit_penalty":de_novo_fit_penalty,
                        "refit_denovo_signatures":refit_denovo_signatures,
                        "collapse_to_SBS96":collapse_to_SBS96,
                        "dist":clustering_distance,
                        "export_probabilities":export_probabilities,
                        "make_decompostion_plots":make_decomposition_plots,
                        "stability":stability, 
                        "min_stability":min_stability, 
                        "combined_stability":combined_stability,
                        "get_all_signature_matrices":get_all_signature_matrices}
    
    
    
    ################################ take the inputs from the mandatory arguments ####################################
    input_type = input_type;
     
    #project = input_data   #the variable was already set above
        
    
    ################################ take the inputs from the general optional arguments ####################################
    startProcess=minimum_signatures ; 
    endProcess=maximum_signatures;
    
    #totalIterations=nmf_replicates
    cpu = cpu
    hierarchy = False #No use
    mtype=context_type
    #init=nmf_init
    wall=get_all_signature_matrices
    add_penalty=nnls_add_penalty
    remove_penalty=nnls_remove_penalty
    genome_build=opportunity_genome
    refgen=reference_genome
    refit_denovo_signatures
    #set the squence type ("genome" or "exome") for the tmb plot inside the make_final_solution function
    if exome==False:
        sequence="genome"
    if exome==True:
        sequence="exome"
    
    #setting seeds
    if seeds=="random":
         excecution_parameters["seeds"]=seeds
         replicates=list(range(1,nmf_replicates+1))
         seed=np.random.randint(0, 10000000, size=nmf_replicates)
         seeds=pd.DataFrame(list(zip(replicates, seed)), columns=["Replicates","Seeds"])
         seeds=seeds.set_index("Replicates")
         seeds.to_csv(out_put+"/Seeds.txt", sep="\t")
    else:
        try:
            excecution_parameters["seeds"]=seeds
            seeds=pd.read_csv(seeds,sep="\t", index_col=0)
            seeds.to_csv(out_put+"/Seeds.txt", sep="\t")
            seed=np.array(seeds["Seeds"])
            
            
            
        except:
            "Please set valid seeds"
    
   
    
    
    
    if input_type=="text" or input_type =="table" or input_type=="matrix":
        
        ################################### For text input files ######################################################
        
        text_file = project
        title = "" # set the title for plotting 
            
    
        if type(text_file)!=str:
            data=text_file
            excecution_parameters["input_data"]="Matrix["+str(data.shape[0])+" rows X "+str(data.shape[1])+ " columns]"
        else:
            data = pd.read_csv(text_file, sep="\t").iloc[:,:]
        
        
        if data.shape[0]==48:
            paths = cosmic.__path__[0]
            feature_map=pd.read_csv(paths+"/data/CN_classes_dictionary.txt", sep="\t", header=None)
            feature_order=pd.read_csv(paths+"/data/CNV_features.tsv", sep="\t", header=None)
            if list(data.iloc[:,0])==list(feature_order[0]):
                pass
            else:
                orderlist1=list(feature_map[0])
                orderlist2=list(feature_order[0])
                #sort the mutation types first step
                data["Mutation Types"]= pd.Categorical(data["Mutation Types"], orderlist1)
                data = data.sort_values("Mutation Types")
                data=data.reset_index()
                data=data.drop(columns='index')
                #sort the mutation types second step
                data["Mutation Types"]=feature_map[1]
                data["Mutation Types"]= pd.Categorical(data["Mutation Types"], orderlist2)
                data = data.sort_values("Mutation Types")
                
        
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
        elif mtypes[0] == "48":
            mtypes = ["CNV48"]
        elif mtypes[0]=="32":
            mtypes = ["SV32"]
        elif mtypes[0]=="96" or "288" or "384" or "1536":
            mtypes = ["SBS"+mtypes[0]]
        else:
            mtypes = ["CH"+mtypes[0]]
        
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
            mtypes = ["SBS96", "DBS78", "ID83"] 
        elif mtype == "default":
            mtypes = ["SBS96", "DBS78", "ID83"] 
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
    
    #recording context types
    excecution_parameters["context_type"]=",".join(mtypes) 
  
    
    record_parameters(sysdata, excecution_parameters, start_time)
    sysdata.close()      
    
    
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
        elif m.startswith("CNV"):
            m="CNV"
        elif m.startswith("SV"):
            m="SV"
        
        # Determine the types of mutation which will be needed for exporting and copying the files
        if not (m=="DINUC" or m.startswith("DBS") or m.startswith("ID") or m.startswith("CNV") or m.startswith("SV")):
            
            if m.startswith("SBS"):
                mutation_type = m
            elif m in ["96","288","384","1536"]:
                mutation_type="SBS"+m
            elif m.startswith("78"): 
                mutation_type="DBS78"
            elif m.startswith("83"):
                mutation_type="ID83"
            elif m.startswith("48"):
                mutation_type="CNV48"
            elif m.startswith("32"):
                mutation_type="SV32"
            else:
                mutation_type = "CH"+m
            
        else:
            if m == "DINUC" or m.startswith("DBS"):
                mutation_type = "DBS78"
            elif m== "ID" or m.startswith("ID"):
                mutation_type = "ID83"
            elif m== "CNV" or m.startswith("CNV"):
                mutation_type = "CNV48"
            elif m== "SV" or m.startswith("SV"):
                mutation_type = "SV32"
                
       
            
        if input_type=="vcf":
            
            try: 
                
                genomes = pd.DataFrame(data[m])
            except KeyError: 
                sysdata = open(out_put+"/JOB_METADATA.txt", "a")
                sysdata.write("Context {} is not available in the current vcf files".format(m)+"\n")
                print("Context {} is not available in the current vcf files".format(m))
                sysdata.close()
                
                continue
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
        
        
        
        
        # get the cutoff for normatization to handle the hypermutators 
        
        normalization_cutoff = sub.get_normalization_cutoff(genomes, manual_cutoff=100*genomes.shape[0])
        #print("Normalization Cutoff is :", normalization_cutoff)
        excecution_parameters["normalization_cutoff"]= normalization_cutoff
        
        #pass the seed values to inner funtions:
        excecution_parameters["seeds"]= seed
            
        
        
        
        if genomes.shape[1]<endProcess:
            endProcess=genomes.shape[1]
        
        #report the notmatlization criteria
        sysdata = open(out_put+"/JOB_METADATA.txt", "a")
        context_start_time=datetime.datetime.now()
        sysdata.write("\n##################################\n")
        sysdata.write("\n[{}] Analysis started for {}. Matrix size [{} rows x {} columns]\n".format(str(context_start_time).split(".")[0],mutation_type,genomes.shape[0],genomes.shape[1])) 
        if excecution_parameters["matrix_normalization"]=="gmm":
                sysdata.write("\n[{}] Normalization GMM with cutoff value set at {}\n". \
                              format(str(datetime.datetime.now()).split(".")[0], normalization_cutoff)) 
        elif excecution_parameters["matrix_normalization"]=="100X":
                sysdata.write("\n[{}] Normalization 100X with cutoff value set at {}\n". \
                              format(str(datetime.datetime.now()).split(".")[0],(genomes.shape[0]*100)))
        elif excecution_parameters["matrix_normalization"]=="log2":
            sysdata.write("\n[{}] Normalization Log2\n". \
                              format(str(datetime.datetime.now()).split(".")[0]))
        elif excecution_parameters["matrix_normalization"]=="none":
            sysdata.write("\n[{}] Analysis is proceeding without normalization\n". \
                          format(str(datetime.datetime.now()).split(".")[0]))
        else:
            sysdata.write("\n[{}] Normalization Custom with cutoff value set at {}\n". \
                              format(str(datetime.datetime.now()).split(".")[0],excecution_parameters["matrix_normalization"]))
            
            
        sysdata.close()        
        
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
            processes = sub.decipher_signatures(excecution_parameters,
                                                genomes= genomes, 
                                                mut_context=m,
                                                i = i)
            
            
            
            
            
            
            
            
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
                #sysdata.write("\nAnalysis of context type {} is ended successfully\n".format(m)) 
            #report progress to the system file:
            
            
            #Get total mutationation for each signature in reverse order and order the signatures from high to low mutation barden
            signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
            
            sorted_idx = np.argsort(-signature_total_mutations)
            processAvg = np.take(processAvg, sorted_idx, axis=1)
            exposureAvg = np.take(exposureAvg, sorted_idx, axis=0)
            signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
            processStd=np.take(processStd, sorted_idx, axis=1)
            exposureStd=np.take(exposureStd, sorted_idx, axis=0)
            clusterSilhouetteCoefficients=np.take(clusterSilhouetteCoefficients, sorted_idx, axis=0)
            
            
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
            sub.export_information(loopResults, m, layer_directory, index, colnames, wall=wall, sequence=sequence)
            
          
            
            all_similirities_list.append(all_similarities)
                #
            #similarity_dataframe["Total Signatures "+str(processes)] = cosine_similarities
            
            current_time_end = datetime.datetime.now()
            
            sysdata = open(out_put+"/JOB_METADATA.txt", "a")
            sysdata.write("\n[{}] {} de novo extraction completed for a total of {} signatures! \nExecution time:{}\n". \
                          format(str(datetime.datetime.now()).split(".")[0],mutation_type,processes,str(current_time_end-current_time_start).split(".")[0], current_time_end))
            sysdata.close()
            
        
        ################################################################################################################
        ########################################## Plot Stabiltity vs Reconstruction Error #############################        
        ################################################################################################################    
        # Print the Stabiltity vs Reconstruction Error as get the solution as well
        solution, all_stats = sub.stabVsRError(layer_directory+"/All_solutions_stat.csv", layer_directory, title, all_similirities_list, mtype=mutation_type, stability=stability, min_stability=min_stability, combined_stability=combined_stability)
        all_stats.insert(1, 'Stability (Avg Silhouette)', minimum_stabilities) #!!!!!!!!!!!!!!!!1 here minimum stability is avg stability
        all_stats=all_stats.set_index(["Signatures"])
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
        exposureAvg = information[solution-startProcess][1]
        processSTE = information[solution-startProcess][2]
        signature_stabilities = information[solution-startProcess][4]
        signature_total_mutations = information[solution-startProcess][5]  
        signature_stats = information[solution-startProcess][6] 
        all_similarities = information[solution-startProcess][7]
        
       
        # create the folder for the final solution/ De Novo Solution
        layer_directory1 = output+"/Suggested_Solution/"+mutation_type+"_De-Novo_Solution"
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
                       signature_total_mutations = signature_total_mutations,denovo_exposureAvg  = exposureAvg, \
                       signature_stats = signature_stats, add_penalty=add_penalty, remove_penalty=remove_penalty, \
                       initial_remove_penalty=initial_remove_penalty, refit_denovo_signatures=refit_denovo_signatures, \
                       de_novo_fit_penalty=de_novo_fit_penalty, sequence=sequence)    
          
        #try:
        # create the folder for the final solution/ Decomposed Solution
        
        layer_directory2 = output+"/Suggested_Solution/COSMIC_"+mutation_type+"_Decomposed_Solution"
        try:
            if not os.path.exists(layer_directory2):
                os.makedirs(layer_directory2)
        except: 
            print ("The {} folder could not be created".format("output"))
            
        originalProcessAvg=pd.DataFrame(processAvg, index=index)
        
        if processAvg.shape[0]==1536 and collapse_to_SBS96==True: #collapse the 1596 context into 96 only for the deocmposition 
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[1:8]).sum()
            index = genomes.index
            processAvg = np.array(processAvg)
            genomes = np.array(genomes)
            
            
        if processAvg.shape[0]==288 and collapse_to_SBS96==True: #collapse the 288 context into 96 only for the deocmposition 
            processAvg = pd.DataFrame(processAvg, index=index)
            processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
            genomes = pd.DataFrame(genomes, index=index)
            genomes = genomes.groupby(genomes.index.str[2:9]).sum()
            index = genomes.index
            processAvg = np.array(processAvg)
            genomes = np.array(genomes)
        
        originalProcessAvg.columns = listOfSignatures    
        final_signatures = sub.signature_decomposition(processAvg, m, layer_directory2, genome_build=genome_build, cosmic_version=cosmic_version, add_penalty=0.05, remove_penalty=0.01, mutation_context=mutation_context, make_decomposition_plots=make_decomposition_plots, originalProcessAvg=originalProcessAvg)
        
        # extract the global signatures and new signatures from the final_signatures dictionary
        globalsigs = final_signatures["globalsigs"]
        globalsigs = np.array(globalsigs)
        newsigs = final_signatures["newsigs"]
        try:    
            processAvg = np.hstack([globalsigs, newsigs])  
            allsigids = final_signatures["globalsigids"]+final_signatures["newsigids"]
        except: 
            processAvg=newsigs
            allsigids=final_signatures["newsigids"]
            
        attribution = final_signatures["dictionary"]
        background_sigs= final_signatures["background_sigs"]
        genomes = pd.DataFrame(genomes)
        
        
        exposureAvg = sub.make_final_solution(processAvg, genomes, allsigids, layer_directory2, m, index, colnames, \
                                cosmic_sigs=True, attribution = attribution, denovo_exposureAvg  = exposureAvg , \
                                background_sigs=background_sigs, add_penalty=add_penalty, remove_penalty=remove_penalty, \
                                initial_remove_penalty=initial_remove_penalty, genome_build=genome_build, \
                                collapse_to_SBS96=collapse_to_SBS96,sequence=sequence,export_probabilities=export_probabilities)
        
            
           
           
    
    sysdata = open(out_put+"/JOB_METADATA.txt", "a")
    end_time = datetime.datetime.now()
    sysdata.write("\n[{}] Analysis ended: \n".format(str(end_time).split(".")[0]))
    
    sysdata.write("\n-------Job Status------- \n")
    sysdata.write("Analysis of mutational signatures completed successfully! \nTotal execution time: "+str(end_time-start_time).split(".")[0]+" \nResults can be found in: "+" "+out_put+ " " +" folder")
    sysdata.close()

    print("\n\n \nYour Job Is Successfully Completed! Thank You For Using SigProfiler Extractor.\n ")
             


     
