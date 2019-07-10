# SigProfilerExtractor
SigProfilerExtractor allows de novo extraction of mutational signatures from data generated in a matrix format. 
The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability 
for each signature to cause a specific mutation type in a cancer sample. The tool makes use of SigProfilerMatrixGenerator 
and SigProfilerPlotting. 

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

    Example: 
    -------
    >>> from sigproextractor import sigpro as sig
    >>> data = sig.importdata("text")
    This "data" variable can be used as a parameter of the "project" argument of the sigProfilerExtractor function.
    
    To get help on the parameters and outputs of the "importdata" function, please write down the following line:
    
    >>> help(sig.importdata)
        

### sigProfilerExtractor 
    
    
    Extracts mutational signatures from an array of samples.
    
    sigProfilerExtractor(input_type, out_put, project, refgen="GRCh37", startProcess=1, endProcess=10, totalIterations=8, 
    cpu=-1, hierarchy = False, mtype = ["default"],exome = False, indel_extended = False)
    
    
    
    Examples
    --------
    
    >>> from sigproextractor import sigpro as sig
    >>> data = sig.importdata("vcf")
    >>> sig.sigProfilerExtractor("vcf", "example_output", data, startProcess=1, endProcess=3)
    
    Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
    Check the current working directory for the "example_output" folder.
    
    To get help on the parameters and outputs of the "sigProfilerExtractor" function, please write down the following line:
    
    >>> help(sig.sigProfilerExtractor)
    
    
## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam.ucsd.edu
