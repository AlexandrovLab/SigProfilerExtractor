import os
import pytest
import shutil

from SigProfilerExtractor import sigpro as sig


SPE_PATH = os.path.dirname(os.path.abspath(__file__))
SPE_PATH_SBS_INPUT = SPE_PATH + "/TextUnorderedInput/Samples_SBS.txt"
SPE_PATH_SBS_OUTPUT = SPE_PATH + "/TextOutput/SBS_output"

def test_output_text():
    # Specify the paths for actual and expected output files
    actual_reference_file = SPE_PATH + '/Reference_Signature_output/SBS96_De-Novo_Signatures.txt'
    expected_output_file = SPE_PATH_SBS_OUTPUT + '/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'

    # Run the function that generates the output
    sig.sigProfilerExtractor("matrix", SPE_PATH_SBS_OUTPUT, SPE_PATH_SBS_INPUT, seeds = SPE_PATH + "/seeds/SBS_seeds/Seeds.txt", minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)

    # Check if the actual output file matches the expected output file
    with open(actual_reference_file, 'r') as actual_file, open(expected_output_file, 'r') as expected_file:
        actual_content = actual_file.read()
        expected_content = expected_file.read()
    
    
    # Remove unwanted files in the SBS_output directory
    shutil.rmtree(SPE_PATH + "/TextOutput/SBS_output/SBS96/All_Solutions")
    os.remove(SPE_PATH + "/TextOutput/SBS_output/SBS96/All_solutions_stat.csv")
    os.remove(SPE_PATH + "/TextOutput/SBS_output/SBS96/Samples.txt")
    os.remove(SPE_PATH + "/TextOutput/SBS_output/SBS96/SBS96_selection_plot.pdf")  
    shutil.rmtree(SPE_PATH + "/TextOutput/SBS_output/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution")
    shutil.rmtree(SPE_PATH + "/TextOutput/SBS_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/")
    shutil.rmtree(SPE_PATH + "/TextOutput/SBS_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Solution_Stats/")
    os.remove(SPE_PATH + "/TextOutput/SBS_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS_96_plots_SBS96_De-Novo.pdf")
    os.remove(SPE_PATH + "/TextOutput/SBS_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures_SEM_Error.txt")

    assert actual_content == expected_content, "Output does not match expected content."


