import os
import pytest
import shutil
import csv
from SigProfilerExtractor import sigpro as sig


SPE_PATH = os.path.dirname(os.path.abspath(__file__))
SPE_PATH_SBS_INPUT = SPE_PATH + "/TextUnorderedInput/Samples_SBS.txt"
SPE_PATH_SBS_OUTPUT = SPE_PATH + "/TextOutput/SBS_output/"

def test_output_text():
    # Specify the paths for actual and expected output files
    actual_reference_file = SPE_PATH + '/Reference_Signature_output/SBS96_De-Novo_Signatures.txt'
    expected_output_file = SPE_PATH_SBS_OUTPUT + '/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'
    # Print the seeds value before running the function
    with open(SPE_PATH + "/seeds/SBS_seeds/Seeds.txt", 'r') as seeds_file:
        seeds_content = seeds_file.read()
        print(seeds_content, "Seeds value before running sigProfilerExtractor:")

    # Run the function that generates the output
    sig.sigProfilerExtractor("matrix", SPE_PATH_SBS_OUTPUT, SPE_PATH_SBS_INPUT, seeds = SPE_PATH + "/seeds/SBS_seeds/Seeds.txt", minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)
    # Print the seeds value after running the function
    with open(SPE_PATH + "/seeds/SBS_seeds/Seeds.txt", 'r') as seeds_file:
        seeds_content = seeds_file.read()
        print(seeds_content, "Seeds value after running sigProfilerExtractor:")

    # Check if the actual output file matches the expected output file
    with open(actual_reference_file, 'r') as actual_file, open(expected_output_file, 'r') as expected_file:
        # actual_content = actual_file.read()
        # expected_content = expected_file.read()
    # Compare headers (columns)
        actual_csv_reader = csv.reader(actual_file, delimiter='\t')
        expected_csv_reader = csv.reader(expected_file, delimiter='\t')
        #reading the headers (the first row) from both the actual and expected CSV files
        actual_header = next(actual_csv_reader)
        expected_header = next(expected_csv_reader)
        assert actual_header == expected_header, "Headers (columns) are not the same."

        # Compare the number of rows
        actual_rows = sum(1 for row in actual_csv_reader)
        expected_rows = sum(1 for row in expected_csv_reader)
        assert actual_rows == expected_rows, "Number of rows is not the same."

    #     # Reset file pointers to compare content in detail
    #     actual_file.seek(0)
    #     expected_file.seek(0)

    #     # Compare content row by row
    #     for actual_row, expected_row in zip(actual_csv_reader, expected_csv_reader):
    #         assert actual_row == expected_row, "Content in rows is not the same."
    
    # print("Test passed: Output matches expected content.")
    
    #assert actual_content == expected_content, "Output does not match expected content."


