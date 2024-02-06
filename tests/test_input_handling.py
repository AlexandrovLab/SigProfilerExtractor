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

    assert actual_content == pytest.approx(expected_content), "Output does not match expected content."

