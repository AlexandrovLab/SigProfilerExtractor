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
     
    **context_type**: A string of mutaion context name/names separated by comma (",") , optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "96,DINUC,ID", where "96" is the SBS96 context, "DINUC"
                  is the DINULEOTIDE context and ID is INDEL context.
                  
    **exome**: Boolean, optional. Defines if the exomes will be extracted. The default value is "False".
    
 
    
    NMF REPLICATES:-
    
    **minimum_signatures**: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1 
    
    **maximum_signatures**: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 25
    
    **nmf_replicates**: A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 100
    
    **resample**: Boolean, optional. Default is True. If True, add poisson noise to samples by resampling.  
    
    **seeds**: A string. It can be used to get reproducible resamples for the NMF replicates. A path of a tab separated .txt file containing the replicated id and preset seeds in a two columns dataframe can be passed through this parameter. The Seeds.txt file in the results folder from a previous analysis can be used for the seeds parameter in a new analysis. The Default value for this parameter is "none". When "none", the seeds for resampling will be random for different analysis.
                  
    
    NMF ENGINES:-
    
    **matrix_normalization**: A string. Method of normalizing the genome matrix before it is analyzed by NMF. Default is value is "gmm". Other options are, "log2", "custom" or "none".           
    
    **nmf_init**: A String. The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'
              Default is 'nndsvd_min'.
    
    **precision**: A string. Values should be single or double. Default is single.
    
    **min_nmf_iterations**: An integer. Value defines the minimum number of iterations to be completed before NMF converges. Default is 10000.
    
    **max_nmf_iterations**: An integer. Value defines the maximum number of iterations to be completed before NMF converges. Default is 1000000.
    
    **nmf_test_conv**: An integer. Value definer the number number of iterations to done between checking next convergence. Default is 10000.
            
    **nmf_tolerance**: A float. Value defines the tolerance to achieve to converge. Default is 1e-15.
    
    
    EXECUTION:-
    
    **cpu**: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors. 
    
    **gpu**:Boolean, optional. Defines if the GPU resource will used if available. Default is False. If True, the GPU resources will be used in the computation.

    **batch_size**: An integer. Will be effective only if the GPU is used. Defines the number of NMF replicates to be performed by each CPU during the parallel processing. Default is 1.
              
    
    SOLUTION ESTIMATION THRESH-HOLDS:-

    **stability**: Float, optional. Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered. 

    **min_stability**: Float, optional. Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered. 

    **combined_stability**: Float, optional. Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered. 
              
    
    DECOMPOSITION:-
    
    **de_novo_fit_penalty**: Float, optional. Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to assign denovo signatures to a sample. 
    
    **nnls_add_penalty**: Float, optional. Takes any positive float. Default is 0.05. Defines the strong (add) thresh-hold cutoff to assign COSMIC signatures to a sample. 
    
    **nnls_remove_penalty**: Float, optional. Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to assign COSMIC signatures to a sample.
     
    **initial_remove_penalty**: Float, optional. Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to COSMIC assign signatures to a sample.
    
    **refit_denovo_signatures**: Boolean, optional. Default is True. If True, then refit the denovo signatures with nnls.
    
    **make_decomposition_plots**: Boolean, optional. Defualt is True. If True, Denovo to Cosmic sigantures decompostion plots will be created as a part the results.
     
    
    OTHERS:-
    
    **get_all_signature_matrices**: A Boolean. If true, the Ws and Hs from all the NMF iterations are generated in the output.
    
    **export_probabilities**: A Boolean. Defualt is True. If False, then doesn't create the probability matrix.
    
    
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
    >>>def main_function():
        # to get input from vcf files
        >>> path_to_example_folder_containing_vcf_files = sig.importdata("vcf")
        >>> data = path_to_example_folder_containing_vcf_files # you can put the path to your folder containing the vcf     samples
        >>> sig.sigProfilerExtractor("vcf", "example_output", data, minimum_signatures=1, maximum_signatures=3)
    >>>if __name__="__main__":
       >>>main_function()
    
    Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
    Check the current working directory for the "example_output" folder.
    
    >>>def main_function():    
       # to get input from table format (mutation catalog matrix)
       >>> path_to_example_table = sig.importdata("matrix")
       >>> data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
       >>> sig.sigProfilerExtractor("matrix", "example_output", data, opportunity_genome="GRCh38", minimum_signatures=1,            maximum_signatures=3)
    >>>if __name__="__main__":
       >>>main_function()
    
    To get help on the parameters and outputs of the "sigProfilerExtractor" function, please write down the following line:
    
    >>> help(sig.sigProfilerExtractor)
```

### Estimation of the optimum solution
Estimate the optimum solution (rank) among different number of solutions (ranks): 
    
    Parameters: 
        
        base_csvfile: A string. Defualt is "All_solutions_stat.csv". Path to a  csv file that contains the statistics of all solutions. 
        All_solution: A string. Default is "All_Solutions". Path to a folder that contains the results of all solutions.
        genomes: A string. Default is Samples.txt. Path to a tab delimilted file that contains the mutation counts for all genomes given to different mutation types.
        output: A string. Default is "results". Path to the output folder.
        title: A string, optional. Default is "Selection_Plot". This sets the title of the selection_plot.pdf
        stability: Float, optional. Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered. 
        min_stability: Float, optional. Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered. 
        combined_stability: Float, optional. Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered. 
   
        
    Example:
         >>>from SigProfilerExtractor import estimate_best_solution as ebs
         >>>ebs.estimate_solution(base_csvfile="All_solutions_stat.csv", 
                    All_solution="All_Solutions", 
                    genomes="Samples.txt", 
                    output="results", 
                    title="Selection_Plot",
                    stability=0.8, 
                    min_stability=0.2, 
                    combined_stability=1.25):
       
         
     Values:
        The files below will be generated in the output folder--
        
        "All_solutions_stat.csv": a  csv file that contains the statistics of all solutions.
        "selection_plot.pdf": a plot that depict the Stability and Mean Sample Cosine Distance for different solutions.



### decompose
    Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples
    
    Parameters: 
        
        signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs. 
        activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
        samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
        output: A string. Path to the output folder.
        de_novo_fit_penalty: Float, optional. Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to be assigned denovo signatures to a sample.
        nnls_add_penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the strong (add) thresh-hold cutoff to be assigned COSMIC signatures to a sample.  
        nnls_remove_penalty: Float, optional. Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample.
        initial_remove_penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to be COSMIC assigned signatures to a sample.
        genome_build = A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
        verbose = Boolean. Prints statements. Default value is False. 
        
    Example:
         >>>from SigProfilerExtractor import decomposition as decomp
         >>>signatures = "path/to/De_Novo_Solution_Signatures.txt"
         >>>activities="path/to/De_Novo_Solution_Activities.txt"
         >>>samples="path/to/Samples.txt"
         >>>output="name or path/to/output"
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
        
### Activity stacked bar plot
    Generates a stacked bar plot showing activities in individuals
    
    Parameters: 
            
        activity_file: The standard output activity file showing the number of, or percentage of mutations attributed to  
        each sample. The row names should be samples while the column names should be signatures.
        output_file: path and full name of the output pdf file, including ".pdf"
        bin_size(optional): Number of samples plotted per page, recommended: 50
        
    Example:
         $ python plotActivity.py 50 sig_attribution_sample.txt test_out.pdf

### GPU support
```
If CUDA out of memory exceptions occur, it will be necessary to reduce the number of CPU processes used (the `cpu` parameter).

## For more information, help and examples, please visit: https://osf.io/t6j7u/wiki/home/

## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam@ucsd.edu
