[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/t6j7u/wiki/home/) 
[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Build Status](https://travis-ci.com/AlexandrovLab/SigProfilerExtractor.svg?branch=master)](https://travis-ci.com/AlexandrovLab/SigProfilerExtractor)

# SigProfilerExtractor
SigProfilerExtractor allows de novo extraction of mutational signatures from data generated in a matrix format. 
The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability 
for each signature to cause a specific mutation type in a cancer sample. The tool makes use of SigProfilerMatrixGenerator 
and SigProfilerPlotting. 

## INSTALLATION
In the commandline, please type the following line:
```
$pip install SigProfilerExtractor
```
Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37')
```
This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish.

open a python interpreter and import the SigProfilerExtractor module. Please see the examples of the functions. 

## FUNCTIONS

### importdata 
    
    
    Imports the path of example data.
    
    importdata(datatype="matrix")

    Example: 
    -------
    >>> from SigProfilerExtractor import sigpro as sig
    >>> path_to_example_table = sig.importdata("matrix")
    >>> data = path_to_example_table 
    This "data" variable can be used as a parameter of the "project" argument of the sigProfilerExtractor function.
    
    To get help on the parameters and outputs of the "importdata" function, please write down the following line:
    
    >>> help(sig.importdata)
        

### sigProfilerExtractor 
    
    
Extracts mutational signatures from an array of samples.

sigProfilerExtractor(input_type, out_put, input_data, reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, 
                         minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, resample = True, batch_size=1, cpu=-1, gpu=False, 
                         nmf_init="alexandrov-lab-custom", precision= "single", matrix_normalization= "100X", seeds= "none", 
                         min_nmf_iterations= 10000, max_nmf_iterations=1000000, nmf_test_conv= 10000, nmf_tolerance= 1e-15, nnls_penalty=0.05, get_all_signature_matrices= False): 


INPUT DATA:-
    
    **input_type**: A string. Type of input. The type of input should be one of the following:
            - "vcf": used for vcf format inputs.
            - "matrix": used for table format inputs using a tab seperated file.
             
        
    **out_put**: A string. The name of the output folder. The output folder will be generated in the current working directory. 
            
    **input_data**: A string. Name of the input folder (in case of "vcf" type input) or the input file (in case of "table"  type input). The project file or folder should be inside the current working directory. For the "vcf" type input,the project has to be a folder which will contain the vcf files in vcf format or text formats. The "text"type projects have to be a file.   
            
    **reference_genome**: A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".
       
    **opportunity_genome**: The build or version of the reference signatures for the reference genome. The default opportunity genome is GRCh37. If the input_type is "vcf", the genome_build automatically matches the input reference genome value.    
     
    **exome**: Boolean, optional. Defines if the exomes will be extracted. The default value is "False".
    
    
    NMF REPLICATES:-
    
    **minimum_signatures**: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1 
    
    **maximum_signatures**: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 10
    
    **nmf_replicates**: A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 100
    
    **resample**: Boolean, optional. Default is True. If True, add poisson noise to samples by resampling.  
    
    **seeds**: A string. It can be used to get reproducible resamples for the NMF replicates. A path of a tab separated .txt file containing the replicated id and preset seeds in a two columns dataframe can be passed through this parameter. The Seeds.txt file in the results folder from a previous analysis can be used for the seeds parameter in a new analysis. The Default value for this parameter is "none". When "none", the seeds for resampling will be random for different analysis.
                  
    
    NMF ENGINES:-
    
    **matrix_normalization**: A string. Method of normalizing the genome matrix before it is analyzed by NMF. Default is value is "100X". Other options are "gmm", "log2", "custom" or "none".           
    
    **nmf_init**: A String. The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'alexandrov-lab-custom'
              Default is 'alexandrov-lab-custom'.
    
    **precision**: A string. Values should be single or double. Default is single.
    
    **min_nmf_iterations**: An integer. Value defines the minimum number of iterations to be completed before NMF converges. Default is 10000.
    
    **max_nmf_iterations**: An integer. Value defines the maximum number of iterations to be completed before NMF converges. Default is 1000000.
    
    **nmf_test_conv**: An integer. Value definer the number number of iterations to done between checking next convergence. Default is 10000.
            
    **nmf_tolerance**: A float. Value defines the tolerance to achieve to converge. Default is 1e-15.
    
    
    EXECUTION:-
    
    **cpu**: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors. 
    
    **gpu**:Boolean, optional. Defines if the GPU resource will used if available. Default is False. If True, the GPU resources will be used in the computation.

    **batch_size**: An integer. Will be effective only if the GPU is used. Defines the number of NMF replicates to be performed by each CPU during the parallel processing. Default is 1.
              
              
    COSMIC:-
    
    **nnls_penalty**: Float, optional. Takes any positive float. Default is 0.05. Defines the thresh-hold cutoff to be assigned signatures to a sample. 
    
      **context_type**: A list of strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is ["96", "DINUC" , "ID"], where "96" is the SBS96 context, "DINUC"
            is the DINULEOTIDE context and ID is INDEL context. 
    
    OTHERS:-
    
    **get_all_signature_matrices**: A Boolean. If true, the Ws and Hs from all the NMF iterations are generated in the output.
    
    Returns
    -------
    To learn about the output, please visit https://osf.io/t6j7u/wiki/home/
    
Returns
-------
To learn about the output, please visit https://osf.io/t6j7u/wiki/home/  
    
```    
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
    
    To get help on the parameters and outputs of the "sigProfilerExtractor" function, please write down the following line:
    
    >>> help(sig.sigProfilerExtractor)
```

### decompose
    Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples
    
    Parameters: 
        
        signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation                         types and colunms are signature IDs. 
        activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs                       and colunms are signature IDs.
        samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types                  and colunms are sample IDs.
        output: A string. Path to the output folder.
        genome_build = A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
        verbose = Boolean. Prints statements. Default value is False. 
        
    Example:
         >>>from SigProfilerExtractor import decomposition as decomp
         >>>signatures = "path/to/dDe_Novo_Solution_Signatures.txt"
         >>>activities="path/to/De_Novo_Solution_Activities.txt"
         >>>samples="path/to/Samples.txt"
         >>>output="name or path/to/output.txt"
         >>>decomp.decompose(signatures, activities, samples, output, genome_build="GRCh37", verbose=False)   
         
     Values:
        The files below will be generated in the output folder--
        
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
        
### GPU support
```
If CUDA out of memory exceptions occur, it will be necessary to reduce the number of CPU processes used (the `cpu` parameter).

## For more information, help and examples, please visit: https://osf.io/t6j7u/wiki/home/

## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam.ucsd.edu
