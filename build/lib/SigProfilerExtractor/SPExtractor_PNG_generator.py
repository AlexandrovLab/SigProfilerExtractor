from SigProfilerExtractor import sigproEdited as sig
import os
import sys
import shutil

param_msg = """py_graph.py requires 5 inputs parameters:
(1) String - Path to directory containing SBS-96 SigProfilerMatrixGenerator formatted file Samples.txt
(2) String - The name of the directory for storing the output
(3) Integer - The start process number
(4) Integer - The end process number
(5) Integer - The number of iterations"""

#if len(sys.argv) != 6:
#    sys.exit(param_msg)

#input_dir = sys.argv[1] # Directory containing Samples.txt
#output_dir = sys.argv[2] # Directory to save the output to (inside of input_dir)
#start_proc = int(sys.argv[3])
#end_proc = int(sys.argv[4])
#total_iter = int(sys.argv[5])
# context = {"SBS", "DBS", "ID"}

png_outputdir_name = "cosmic_and_denovo_pngs"

SBS96_cosmic_output = "cosmic_and_denovo_pngs/SPExtractor_SBS96_Cosmic_png/"
SBS96_denovo_output = "cosmic_and_denovo_pngs/SPExtractor_SBS96_Signature_png/"
SBS96_decomp_path = "/SBS96/Suggested_Solution/Decomposed_Solution/"

SBS1536_cosmic_output = "cosmic_and_denovo_pngs/SPExtractor_SBS1536_Cosmic_png/"
SBS1536_denovo_output = "cosmic_and_denovo_pngs/SPExtractor_SBS1536_Signature_png/"
SBS1536_decomp_path = "/SBS1536/Suggested_Solution/Decomposed_Solution/"

DBS78_cosmic_output = "cosmic_and_denovo_pngs/SPExtractor_DBS78_Cosmic_png/"
DBS78_denovo_output = "cosmic_and_denovo_pngs/SPExtractor_DBS78_Signature_png/"
DBS78_decomp_path = "/DBS78/Suggested_Solution/Decomposed_Solution/"

ID83_cosmic_output = "cosmic_and_denovo_pngs/SPExtractor_ID83_Cosmic_png/"
ID83_denovo_output = "cosmic_and_denovo_pngs/SPExtractor_ID83_Signature_png/"
ID83_decomp_path = "/ID83/Suggested_Solution/Decomposed_Solution/"

def SBS96_gen_dirs():
    if not os.path.exists(SBS96_cosmic_output):
        os.mkdir(SBS96_cosmic_output)
    if not os.path.exists(SBS96_denovo_output):
        os.mkdir(SBS96_denovo_output)

def SBS1536_gen_dirs():
    if not os.path.exists(SBS96_cosmic_output):
        os.mkdir(SBS96_cosmic_output)
    if not os.path.exists(SBS1536_denovo_output):
        os.mkdir(SBS1536_denovo_output)

def DBS78_gen_dirs():
    if not os.path.exists(DBS78_cosmic_output):
        os.mkdir(DBS78_cosmic_output)
    if not os.path.exists(DBS78_denovo_output):
        os.mkdir(DBS78_denovo_output)

def ID83_gen_dirs():
    if not os.path.exists(ID83_cosmic_output):
        os.mkdir(ID83_cosmic_output)
    if not os.path.exists(ID83_denovo_output):
        os.mkdir(ID83_denovo_output)


# def SBS96_mv_dirs(output_dir):
#
#     try:
#         shutil.move(png_outputdir_name, output_dir + SBS96_decomp_path)
#     except:
#         shutil.rmtree(output_dir + SBS96_decomp_path + png_outputdir_name + "/" )
#         shutil.move(SBS96_cosmic_output, output_dir + SBS96_decomp_path)

def gen_png_plots(input_file_path, output_dir, start_proc, end_proc, total_iter, context):
    last_dir_loc = input_file_path.rindex("/")
    input_dir = input_file_path[0:last_dir_loc+1]
    input_file = input_file_path[last_dir_loc+1:]

    cwd = os.getcwd()

    # check if the directory already exists
    os.chdir(input_dir)
    if not os.path.exists(png_outputdir_name):
        os.mkdir(png_outputdir_name)

    if "SBS96" in context:
        SBS96_gen_dirs()
    elif "SBS1536" in context:
        SBS1536_gen_dirs()
    elif "DBS78" in context:
        DBS78_gen_dirs()
    elif "ID83" in context:
        ID83_gen_dirs()
    else:
        print("Program only supports SBS96, SBS1536, DBS78, and ID83 contexts.")
        os.chdir(cwd)
        return



    sig.sigProfilerExtractor("table", output_dir, input_file, genome_build="GRCh38", \
     startProcess=start_proc, endProcess=end_proc, totalIterations=total_iter)

    #SBS96_mv_dirs(output_dir)

    os.chdir(cwd)
