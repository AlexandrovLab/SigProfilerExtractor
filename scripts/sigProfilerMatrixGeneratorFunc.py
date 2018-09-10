#!/usr/bin/env python3
import sigProfilerMatrixGenerator_withBED as matGen
import os
import re

def sigProfilerMatrixGeneratorFunc (project, genome, exome=False, indel=False, indel_extended=False, bed_file=None):
	'''
	Allows for the import of the sigProfilerMatrixGenerator.py function. Returns a dictionary
	with each context serving as the first level of keys. 

	Parameters:
			   project  -> unique name given to the current samples
				genome  -> reference genome 
				 exome  -> flag to use only the exome or not
				 indel  -> flag to create the matrix for the limited indel matrix
		indel_extended  -> flag to create the regular, extended INDEL matrix
			  bed_file  -> BED file that contains a list of ranges to be used in generating the matrices

	Returns:
			  matrices  -> dictionary (nested) of the matrices for each context

		example:
			matrices = {'96': {'PD1001a':{'A[A>C]A':23,
										 'A[A>G]A':10,...},
							  'PD1202a':{'A[A>C]A':23,
										 'A[A>G]A':10,...},...},
						'192':{'PD1001a':{'T:A[A>C]A':23,
										 'T:A[A>G]A':10,...},
							  'PD1202a':{'T:A[A>C]A':23,
										 'T:A[A>G]A':10,...},...},...}
	'''

	# Instantiates all of the required variables
	functionFlag = True
	bed = False
	bed_ranges = None
	
	matrices = {'96':None, '1536':None, '192':None, '3072':None, 'DINUC':None, '6':None, '12':None}

	ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
				  'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
				  'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
				  'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
				  'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
				  'NC_000087.7':'Y'}

	contexts = ['96', '1536', '192', '3072', 'DINUC']

	# Organizes all of the reference directories for later reference:
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	chrom_path = ref_dir + '/references/chromosomes/chrom_string/' + genome + "/"
	chromosome_TSB_path =  ref_dir + '/references/chromosomes/tsb/' + genome + "/"

	# Organizes all of the input and output directories:
	output_matrix = ref_dir + "/references/matrix/"
	vcf_path = ref_dir + '/references/vcf_files/' + project + "/"

	# Gathers all of the vcf files:
	vcf_files2 = os.listdir(vcf_path)
	vcf_files = []

	for file in vcf_files2:
		# Skips hidden files
		if file[0:3] == '.DS':
			pass
		else:
			vcf_files.append(file)

	file_name = vcf_files[0].split(".")
	file_extension = file_name[-1]

	output_path = ref_dir + "/references/vcf_files/single/"
	if not os.path.exists(output_path):
		os.makedirs(output_path)

	if file_extension == 'genome':
		os.system("bash convert_txt_files_to_simple_files.sh " + project)
	else:
		os.system("bash convert_" + file_extension + "_files_to_simple_files.sh " + project)
	vcf_files = os.listdir(ref_dir + '/references/vcf_files/single/')
	vcf_path = ref_dir + '/references/vcf_files/single/'

	# Include some kind of a flag for the INDEL option 
	sort_command_initial = "sort -t $'\t' -k 6,6n -k 6,6 -k 2,2 -k 7,7n "
	sort_command_initial_2 = " -o "
	sort_file = vcf_files[0]
	os.system(sort_command_initial + vcf_path + sort_file + sort_command_initial_2 + vcf_path + sort_file)

	print("Sorting complete...\nDetermining mutation type for each variant, one chromosome at a time.")

	# Adjusts the variables if the context is for INDELs
	if indel:
		contexts = ['INDEL']
		matrices = {'INDEL':None}

	print("Starting catalogue generation...")

	if bed_file != None:
		bed = True
		bed_file_path = ref_dir + "/references/vcf_files/BED/" + project + "/" + bed_file
		bed_ranges = matGen.BED_filtering(bed_file_path)

	# Creates the matrix for each context
	for context in contexts:
		if context != 'DINUC' and context != 'INDEL':
			matrix = matGen.catalogue_generator_single (vcf_path, vcf_files, chrom_path, chromosome_TSB_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges)
			matrices[context] = matrix

		elif context == 'DINUC':
			matrix = matGen.catalogue_generator_DINUC_single (vcf_path, vcf_files, chrom_path, chromosome_TSB_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges)
			matrices[context] = matrix

		elif context == 'INDEL':
			matrix = matGen.catalogue_generator_INDEL_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, indel_extended,functionFlag, bed, bed_ranges)
			matrices[context] = matrix

	# Deletes the temporary files and returns the final matrix
	print("Catalogue for " + context + " context is complete.")
	os.system("rm -r " + vcf_path)
	return(matrices)

	