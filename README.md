# SigProfilerExtractor
SigProfilerExtractor allows de novo extraction of mutational signatures from data generated in a matrix format. 
The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability 
for each signature to cause a specific mutation type in a cancer sample. The tool makes use of SigProfilerMatrixGenerator 
and SigProfilerPlotting. 

## PREREQUISITES

This module requires all the prerequisites of the 
[SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator/blob/master/README.md) and 
[SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting/blob/master/README.md) plus the following
packages referred in the table below:

| Name of packages | Version     | How to install  |
| ------------- |:-------------:| :-----:|
| nimfa         | 1.3.4 or more | ``` conda install -c cdeepakroy nimfa   ```|
| sigProfilerPlotting | 1.1 or more      |  ``` conda install sigProfilerPlotting ```|
| scipy          | 1.1 or more      |  ``` conda install scipy ``` |
| scikit-learn   | 0.19.1 or more | ``` conda install scikit-learn ```|
| SigProfilerMatrixGenerator |       |  ``` conda install SigProfilerMatrixGenerator ```|

"Conda" is a cross-platform, open source package manager that can be used to install different versions of software packages and libraries. To get help to install and use conda, please visit [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator/blob/master/README.md).SigProfilerExtractor supports all operating system. 

## INSTALLATION
In the commandline, please type the following line:
```
$pip install sigproextractor
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
    
    importdata(datatype="matobj")
    
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
        
    



### sigProfilerExtractor 
    
    
    Extracts mutational signatures from an array of samples.
    
    sigProfilerExtractor(input_type, out_put, project, refgen="GRCh37", startProcess=1, endProcess=10, totalIterations=8, 
    cpu=-1, hierarchy = False, mtype = ["default"],exome = False, indel_extended = False)
    
    Parameters
    ----------
    
    input_type: A string. Type of input. The type of input should be one of the following:
            - "vcf": used for vcf format inputs.
            - "text": used for text format inputs.
            - "matobj": used for matlab object format of inputs. 
        
    out_put: A string. The name of the output folder. The output folder will be generated in the current working
            directory. 
            
    project: A string. Name of the input folder (in case of "vcf" type input) or the input file 
            (in case of "text" or "matobj" type input). The project file or folder should be inside the current working directory. 
            For the "vcf" type input,the project has to be a folder which will contain the vcf files in vcf format or text formats.
            The "text" or "matobj" type projects have to be a file. "matobj" projects should have .mat extension.  
            
    refgen: A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only 
            if the input_type is "vcf".
            
    startProcess: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1 
    
    endProcess: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 10
    
    totalIterations: A positive integer, optional. The number of iteration to be performed to extract each number signature. 
            The default value is 8
            
    cpu: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will
            use all available processors. 
    
    hierarchy: boolean, optional. Defines if the signature will be extracted in a hierarchical fashion. The default value is "False".
    
    mtype: A list of strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. 
            The default value is ["96", "DINUC" , "INDEL"].
            
    exome: boolean, optional. Defines if the exomes will be extracted. The default value is "False".
    
    indel_extended: boolean, optional. Defines if the extended indels will be extracted. The default value is "False".
    
    
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
    "Analysis" which will further contain the solutions  in the layer (L) subdirectories. Everything else will be similar to
    the previously deccribed directory structures. The structure of the result folder is synopsized below:
        
        If Hierarchy is False:
            
        -Mutational Context folder
            -All solution folder
                -Signature folder
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
            -Final solution folder
                -De Novo Solution folder
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
                    -dendrogram plot file
                -Decomposed Solution folder
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
            -Analysis folder
                -Layer folder (L)
                    -All solution folder
                        -Signature folder
                            -exposure.txt file
                            -signature.txt file
                            -probabilities.txt file
                            -signature plot pdf file
                    -Selected solution folder
                        -exposure.txt file
                        -signature.txt file
                        -probabilities.txt file
                        -signature plot pdf file
                    -results_stat.csv file
                    -stability plot pdf
            -Final solution folder
                -De Novo Solution folder
                    -exposure.txt file
                    -signature.txt file
                    -probabilities.txt file
                    -signature plot pdf file
                    -dendrogram plot file
                -Decomposed Solution folder
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
    
    
## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam.ucsd.edu
