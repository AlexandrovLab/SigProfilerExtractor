import os
import pytest
import shutil
import pandas as pd

from SigProfilerExtractor import sigpro as sig


SPE_PATH = os.path.dirname(os.path.abspath(__file__))
SPE_PATH_SBS_INPUT = SPE_PATH + "/TextUnorderedInput/Samples_SBS.txt"
SPE_PATH_SBS_OUTPUT = SPE_PATH + "/TextOutput/SBS_output/"

##########test ##
# def round_float(number):
#     return round(number, 4)  # Round to 5 decimal places after the fraction 

def test_output_text():
    # Specify the paths for actual and expected output files
    actual_reference_file = SPE_PATH + '/Reference_Signature_output/SBS96_De-Novo_Signatures.txt'
    expected_output_file = SPE_PATH_SBS_OUTPUT + '/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'
    # Read the seeds file and print its contents
    with open(SPE_PATH + "/seeds/SBS_seeds/Seeds.txt", 'r') as seeds_file:
        reference_seeds_content = seeds_file.read()
        print("reference Seeds content:")
        print(reference_seeds_content)
    # Run the function that generates the output
    sig.sigProfilerExtractor("matrix", SPE_PATH_SBS_OUTPUT, SPE_PATH_SBS_INPUT, seeds = SPE_PATH + "/seeds/SBS_seeds/Seeds.txt", minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)
    with open(SPE_PATH + "/TextOutput/SBS_output/Seeds.txt", 'r') as seeds_file:
        actual_seeds_content = seeds_file.read()
        print("actual Seeds content:")
        print(actual_seeds_content)
    # with open(actual_reference_file, 'r') as actual_file, open(expected_output_file, 'r') as expected_file:
    #     actual_lines = actual_file.readlines()
    #     expected_lines = expected_file.readlines()

    # assert len(actual_lines) == len(expected_lines), "Number of lines in actual and expected files are different."

    # for actual_line, expected_line in zip(actual_lines, expected_lines):
    #     # Skip header lines or any non-numeric values
    #     if actual_line.startswith("MutationType") or expected_line.startswith("MutationType"):
    #         continue
        
    #     actual_values = ["{:.4f}".format(float(val)) for val in actual_line.strip().split() if val.replace('.', '', 1).isdigit()]
    #     expected_values = ["{:.4f}".format(float(val)) for val in expected_line.strip().split() if val.replace('.', '', 1).isdigit()]

    #    assert actual_values == pytest.approx(expected_values), f"Output does not match expected content. Actual: {actual_values}, Expected: {expected_values}"
    
    # Read the actual and expected output files as dataframes
    actual_content = pd.read_csv(actual_reference_file, sep='\t')
    expected_content = pd.read_csv(expected_output_file, sep='\t')

    # Check if the actual dataframe matches the expected dataframe
    #assert actual_content.equals(expected_content), "Dataframes are not equal."
    assert actual_content.columns.tolist() == expected_content.columns.tolist(), "Column names are not the same."
