import os
import pytest
import shutil

from SigProfilerExtractor import sigpro as sig


SPE_PATH = os.path.dirname(os.path.abspath(__file__))
SPE_PATH_SBS_INPUT = SPE_PATH + "/TextUnorderedInput/Samples_SBS.txt"
# Define the path of the TextOutput folder within SPE_PATH
text_output_folder = os.path.join(SPE_PATH, "TextOutput")
# Define the path of the SBS_output folder within the TextOutput folder
SPE_PATH_SBS_OUTPUT = os.path.join(text_output_folder, "SBS_output")

def round_float(number):
    return round(number, 5)  # Round to 5 decimal places after the fraction 

def test_output_text():
    # Specify the paths for actual and expected output files
    actual_reference_file = os.path.join(SPE_PATH, "Reference_Signature_output", "SBS96_De-Novo_Signatures.txt")
    expected_output_file = os.path.join(SPE_PATH_SBS_OUTPUT, "SBS96", "Suggested_Solution", "SBS96_De-Novo_Solution", "Signatures", "SBS96_De-Novo_Signatures.txt")

    # Run the function that generates the output
    sig.sigProfilerExtractor("matrix", SPE_PATH_SBS_OUTPUT, SPE_PATH_SBS_INPUT, seeds=os.path.join(SPE_PATH, "seeds", "SBS_seeds", "Seeds.txt"), minimum_signatures=3, maximum_signatures=3, nmf_replicates=5, min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)

    with open(actual_reference_file, 'r') as actual_file, open(expected_output_file, 'r') as expected_file:
        actual_lines = actual_file.readlines()
        expected_lines = expected_file.readlines()

    assert len(actual_lines) == len(expected_lines), "Number of lines in actual and expected files are different."

    for actual_line, expected_line in zip(actual_lines, expected_lines):
        # Skip header lines or any non-numeric values
        if actual_line.startswith("MutationType") or expected_line.startswith("MutationType"):
            continue
        
        actual_values = [round_float(float(val)) for val in actual_line.strip().split() if val.replace('.', '', 1).isdigit()]
        expected_values = [round_float(float(val)) for val in expected_line.strip().split() if val.replace('.', '', 1).isdigit()]

        assert actual_values == expected_values, f"Output does not match expected content. Actual: {actual_values}, Expected: {expected_values}"

