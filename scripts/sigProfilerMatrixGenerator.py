#!/usr/bin/env python3

#This file is part of Mutational Signatures Project.

#Mutational Signatures Project: need info on project

#Copyright (C) 2018 Erik Bergstrom

#

#Mutational Signatures is free software: need to include distribtution

#rights and so forth

#

#Mutational Signatures is distributed in the hope that it will be useful,

#but WITHOUT ANY WARRANTY; without even the implied warranty of

#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

#GNU General Public License for more details [change to whatever group we should include.
 

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import os
import sys
import re
import argparse
import itertools
import pandas as pd
from itertools import chain
import time

################# Functions and references ###############################################
def perm(n, seq):
    '''
    Generates a list of all available permutations of n-mers
    '''
    permus = []
    for p in itertools.product(seq, repeat=n):
        permus.append("".join(p))
    return(permus)

def BED_filtering (bed_file_path):
    ranges = {}
    ranges_final = {}
    with open(bed_file_path) as f:
        next(f)
        for lines in f:
            line = lines.strip().split()
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            if chrom not in ranges.keys():
                ranges[chrom] = []
            ranges[chrom].append((start, end))

    for chroms in ranges.keys():
        ranges_final[chroms] = set(chain(*(range(start, end+1) for start, end in ranges[chroms])))

    return(ranges_final)

def catalogue_generator_single (vcf_path, vcf_files, chrom_path, chromosome_TSB_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges):
    '''
    Generates the mutational matrix for 96, 1536, 192, and 3072 context using a single vcf file with all samples of interest.

    Parameters:
                   vcf_path  -> path to vcf file of interest
                  vcf_files  -> actual vcf file
                 chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
                                file name: '1.txt', '2.txt', etc.
        chromosome_TSB_path  -> path to transcriptional strand information reference for each chromosome. Only necessary if
                                the desired context is 192 or 3072. Use "saveTSB_192.py" script to generate these files.
                    project  -> unique name given to the set of samples (ex. 'BRCA') 
              output_matrix  -> path where the final mutational matrix is stored
                    context  -> desired context (ex. 96, 1536, 192, 3072)
                 ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
                                for the mm10 assembly.
                        bed  -> parameter used to filter the mutations on a user-provided BED file
                 bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file

    Returns:
        None

    Outputs:
        Outputs a mutational matrix for the desired context within /references/matrix/[project]/

    '''

    # Small functions to provide reverse complements of TSB and sequence info:
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
    revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1'}[B] for B in x][::-1])
    
    # Provides the sorting order for the TSB matrices
    bias_sort = {'T':0,'U':1,'N':3,'B':2}
    tsb = ['T','U','N','B']
    bases = ['A','C','G','T']
    
    # Creates all possible mut_types
    size = 3
    if context == '1536' or context == '3072':
        size = 5
    mut_types_initial = perm(size, "ACGT")
    mut_types = []
    if context == '192' or context == '3072':
        for tsbs in tsb:
            for mut in mut_types_initial:
                current_base = mut[int(size/2)]
                if current_base == 'C' or current_base == 'T':
                    for base in bases:
                        if base != current_base:
                            mut_types.append(tsbs+":"+mut[0:int(size/2)] + "[" + current_base+">"+ base+"]"+mut[int(size/2)+1:])
    else:
        for mut in mut_types_initial:
            current_base = mut[int(size/2)]
            if current_base == 'C' or current_base == 'T':
                for base in bases:
                    if base != current_base:
                        mut_types.append(mut[0:int(size/2)] + "[" + mut[int(size/2)]+">"+ base+"]"+mut[int(size/2)+1:])


    types = []
    samples = []
    mutation_dict = {}
    flag = True
    i = 0
    file = vcf_files[0]
    sample_start = None

    if exome:
        exome_temp_file = "exome_temp.txt"
        exome_file = open(exome_temp_file, 'w')

    with open (vcf_path + file) as f:                
            # Skips any header lines
            for lines in f:
                if lines[0] == '#':
                    next(f)

                # Saves all relevant data from the file
                else:
                    line = lines.strip().split('\t')
                    sample = line[1]
                    chrom = line[5]
                    if chrom in ncbi_chrom.keys():
                        chrom = ncbi_chrom[chrom]
                    if len(chrom) > 1:
                        if chrom[0:3].upper() == 'CHR':
                            chrom = chrom[-1]

                    start = int(line[6])
                    ref = line[8]
                    mut = line[9]

                    if bed:
                        if start not in bed_ranges[chrom]:
                            print(chrom, start)
                            continue
 
                    if flag:
                        sample_start = sample
                        chrom_start = chrom

                        # Opens each chromosome string path and bias string path
                        try:
                            with open(chrom_path + chrom_start + ".txt") as f:
                                chrom_string = f.readline().strip()
                            if context == '192' or context == '3072':
                                with open (chromosome_TSB_path + chrom_start + "_192.txt", 'rb') as f:
                                    chrom_bias = f.read()
                        except:
                            print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
                            continue
                        flag = False

                    # Saves the sample name for later reference
                    if sample not in samples:
                        samples.append(sample)
                        mutation_dict[sample] = {}

                    # Opens the next chromosome data once a new chromosome is reached
                    if chrom != chrom_start:
                        print ("Chromosome " + chrom_start + " done")
                        chrom_start = chrom
                        try:
                            with open(chrom_path + chrom_start + ".txt") as f:
                                chrom_string = f.readline().strip()

                            if context == '192' or context == '3072':
                                with open(chromosome_TSB_path + chrom_start + "_192.txt", 'rb') as f:
                                    chrom_bias = f.read()
                            i += 1
                        except:
                            print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
                            continue

                        

                    # Pulls out the relevant sequence depending on the context
                    if context == '96' or context == '192':
                        sequence = chrom_string[start-2:start+1]

                    else:
                        sequence = chrom_string[start-3:start+2]

                    if context == '192' or context == '3072':
                        bias = chrom_bias[start]
                    char = sequence[int(len(sequence)/2)]

 
                    # Prints the sequence and position if the pulled sequence doesn't match
                    # the variant from the file
                    if char != ref and revcompl(char) != ref:
                        print (sequence + "\t" + chrom + "\t" + char + "\t" + ref + "\t" + str(start))
                
                    
                    # Saves the sequence/mutation type if it matched the reference  
                    else:
                        if char == ref:
                            if ref == 'A' or ref == 'G':
                                ref = revcompl(ref)
                                mut = revcompl(mut)
                                sequence = revcompl(sequence)

                        else:
                            ref = revcompl(ref)
                            mut = revcompl(mut) 
                            sequence = revcompl(sequence)
                            if context == '192' or context == '3072':
                                bias = revbias(str(bias))

                        mut_key = sequence[0:int(len(sequence)/2)] + '[' + ref + '>' + mut + ']' + sequence[int(len(sequence)/2+1):]



                        if context == '192' or context == '3072':
                            if bias == 0:
                                bias = 'N'
                            elif bias == 1:
                                bias = 'T'
                            elif bias == 2:
                                bias = 'U'
                            else:
                                bias = 'B'
                            mut_key = bias + ':' + mut_key

                        if mut_key not in types:
                            types.append(mut_key)
                        if mut_key not in mutation_dict[sample].keys():
                            mutation_dict[sample][mut_key] = 1
                        else:
                            mutation_dict[sample][mut_key] += 1

                        if exome:
                            exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + mut_key + "\n")
    print ("Chromosome " + chrom_start + " done")

    if exome:
        exome_file.close()
        os.system("sort -t $'\t' -k 2,2n -k 2,2 -k 3,3n " + exome_temp_file + " -o " + exome_temp_file)
        print("Beginning exome filtering. This may take a few moments...")
        mutation_dict, samples = exome_check (genome, exome_temp_file)
        os.system("rm " + exome_temp_file)

    # Calls the function to generate the final matrix
    if functionFlag:
        return(pd.DataFrame.from_dict(mutation_dict))
    else:
        matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed)


def catalogue_generator_DINUC_single (vcf_path, vcf_files, chrom_path, chromosome_TSB_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges):
    '''
    Generates the mutational matrix for the dinucleotide context.

    Parameters:
                   vcf_path  -> path to the vcf file of interest
                  vcf_files  -> vcf file of interest
                 chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
                                file name: '1.txt', '2.txt', etc.
        chromosome_TSB_path  -> path to transcriptional strand information reference for each chromosome. Only necessary if
                                the desired context is 192 or 3072. Use "saveTSB_192.py" script to generate these files.
                    project  -> unique name given to the set of samples (ex. 'BRCA') 
              output_matrix  -> path where the final mutational matrix is stored
                 ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
                                for the mm10 assembly.
                        bed  -> parameter used to filter the mutations on a user-provided BED file
                 bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file


    Returns:
        None

    Output:
        Saves the mutational matrix for the dinucleotide context

    '''

    # Small functions to provide reverse complements of TSB and sequence info:
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
    revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1'}[B] for B in x][::-1])
    

    # Instantiates the necessary variables/data structures
    dinucs = {}
    samples = []
    mutation_types = ['AA>CC','AA>CG','AA>CT','AA>GC','AA>GG','AA>GT','AA>TC','AA>TG','AA>TT',
                      'AC>CA','AC>CG','AC>CT','AC>GA','AC>GG','AC>GT','AC>TA','AC>TG','AC>TT',
                      'AG>CA','AG>CC','AG>CT','AG>GA','AG>GC','AG>GT','AG>TA','AG>TC','AG>TT',
                      'AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','CA>AC','CA>AG','CA>AT',
                      'CA>GC','CA>GG','CA>GT','CA>TC','CA>TG','CA>TT','CC>AA','CC>AG','CC>AT',
                      'CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AA','CG>AC','CG>AT',
                      'CG>GA','CG>GC','CG>TA','GA>AC','GA>AG','GA>AT','GA>CC','GA>CG','GA>CT',
                      'GA>TC','GA>TG','GA>TT','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA',
                      'TA>AC','TA>AG','TA>AT','TA>CC','TA>CG','TA>GC']


    types = []
    samples = []
    mutation_dict = {}
    flag = True
    i = 0
    file = vcf_files[0]
    sample_start = None

    if exome:
        exome_temp_file = "exome_temp.txt"
        exome_file = open(exome_temp_file, 'w')

    with open (vcf_path + file) as data:
            initial_line = data.readline()
            initial_line_data = initial_line.strip().split('\t')
            previous_sample = initial_line_data[1]
            previous_chrom = initial_line_data[5] 
            previous_start = int(initial_line_data[6])
            previous_ref = initial_line_data[8]
            previous_mut = initial_line_data[9]
                    
            
            for lines in data:
                # Skips any header lines
                if lines[0] == '#':
                    next(data)

                # Saves the relevant data from each line
                else:
                    line = lines.strip().split('\t')
                    sample = line[1]
                    chrom = line[5]
                    if chrom in ncbi_chrom.keys():
                        chrom = ncbi_chrom[chrom]
                    if len(chrom) > 1:
                        if chrom[0:3].upper() == 'CHR':
                            chrom = chrom[-1]

                    start = int(line[6])
                    ref = line[8]
                    mut = line[9]

                    if bed:
                        if start not in bed_ranges[chrom]:
                            continue

                    # Saves the sample name for later reference
                    if sample not in samples:
                        samples.append(sample)

                    if flag:
                        sample_start = sample
                        chrom_start = chrom

                        # Open each chromosome string path and bias string path
                        try:
                            with open(chrom_path + chrom_start + ".txt") as f:
                                chrom_string = f.readline().strip()
                            flag = False
                        except:
                            print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
                            continue

                    # Breaks the loop once a new chromosome is reached
                    if chrom != chrom_start:
                        print ("Chromosome " + chrom_start + " done")
                        chrom_start = chrom
                        with open(chrom_path + chrom_start + ".txt") as f:
                            chrom_string = f.readline().strip()
                        

                    # Grabs the slice of interest
                    else:
                        if start == previous_start + 1:
                            dinuc = previous_ref + ref + ">" + previous_mut + mut

                            if sample not in dinucs.keys():
                                dinucs[sample] = {}
                                for dinucl in mutation_types:
                                    dinucs[sample][dinucl]=0
                            if dinuc in mutation_types:
                                dinucs[sample][dinuc] += 1 
                            else:
                                dinuc = revcompl(previous_ref + ref) + '>' + revcompl(previous_mut + mut)
                                dinucs[sample][dinuc] += 1

                            if exome:
                                exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + dinuc + "\n")

                    previous_sample = sample
                    previous_chrom = chrom 
                    previous_start = start
                    previous_ref = ref
                    previous_mut = mut

    print ("Chromosome " + chrom_start + " done")

    if exome:
        exome_file.close()
        os.system("sort -t $'\t' -k 2,2n -k 2,2 -k 3,3n " + exome_temp_file + " -o " + exome_temp_file)
        print("Beginning exome filtering. This may take a few moments...")
        mutation_dict, samples = exome_check (genome, exome_temp_file)
        os.system("rm " + exome_temp_file)

    # Calls the function to generate the final mutational matrix
    if functionFlag:
        return(pd.DataFrame.from_dict(dinucs))
    else:
        matrix_generator_DINUC (output_matrix, samples, mutation_types, dinucs, project, exome, bed)
        
    

def catalogue_generator_INDEL_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, limited_indel, functionFlag, bed, bed_ranges):
    '''
    Generates the mutational matrix for the INDEL context.

    Parameters:
                   vcf_path  -> path to the vcf file of interest
                  vcf_files  -> vcf file of interest
                 chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
                                file name: '1.txt', '2.txt', etc.
                    project  -> unique name given to the set of samples (ex. 'BRCA') 
              output_matrix  -> path where the final mutational matrix is stored
                 ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
                                for the mm10 assembly.
                        bed  -> parameter used to filter the mutations on a user-provided BED file
                 bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file



    Returns:
        None

    Output:
        Saves the mutational matrix for the INDEL context

    '''

    # Instantiates the necessary variables/data structures
    indel_dict = {}
    samples = []
                        # Single point mutations
    indel_types = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
                   '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
                   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
                   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
                        # >1bp INDELS
                   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
                   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
                   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
                   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
                   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
                   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
                   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
                   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
                        #MicroHomology INDELS
                   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
                   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5', '2:Ins:M:1', 
                   '3:Ins:M:1', '3:Ins:M:2', '4:Ins:M:1', '4:Ins:M:2', '4:Ins:M:3', '5:Ins:M:1', 
                   '5:Ins:M:2', '5:Ins:M:3', '5:Ins:M:4', '5:Ins:M:5', 'complex', 'non_matching']

    if limited_indel:
        indel_types = indel_types[:-13]

    # Instantiates the remaining varibales/data structures
    i = 0
    chrom_string = None
    count = 0
    non_matching = 0
    complex_muts = 0

    if exome:
        exome_temp_file = "exome_temp.txt"
        exome_file = open(exome_temp_file, 'w')

    with open (vcf_path + vcf_files[0]) as data:
        first_flag = True

        # Saves the relevant data from each line
        for lines in data:
            line = lines.strip().split('\t')
            sample = line[1]
            chrom = line[5]
            if chrom in ncbi_chrom.keys():
                chrom = ncbi_chrom[chrom]
            if len(chrom) > 1:
                if chrom[0:3].upper() == 'CHR':
                    chrom = chrom[-1]

            start = int(line[6])
            ref = line[8]
            mut = line[9]
            sub_type_length = 0

            mut_type = None

            if bed:
                if start not in bed_ranges[chrom]:
                    continue

            if first_flag:
                initial_chrom = chrom
                try:
                    with open (chrom_path + initial_chrom + '.txt') as f:
                        chrom_string = f.readline().strip()
                    first_flag = False
                except:
                    print(initial_chrom + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
                    continue

            # Ensures that the variant matches the reference chromosome
            if ref[0] == chrom_string[start-1]:
                # Saves the mutation type for the given variant
                if len(ref) - len(mut) == len(ref)-1:
                    mut_type = 'Del'
                elif len(mut) - len(ref) == len(mut)-1:
                    mut_type = 'Ins'
                else:
                    if 'complex' not in indel_dict[sample].keys():
                        indel_dict[sample]['complex'] = 0
                    else:
                        indel_dict[sample]['complex'] += 1
                    continue
                
                # Opens the next chromosome when a new chromosome is reached in the file
                if chrom != initial_chrom:
                    initial_chrom = chrom
                    try:
                        with open (chrom_path + initial_chrom + '.txt') as f:
                            chrom_string = f.readline().strip()
                        i += 1
                        print(initial_chrom)
                    except:
                        print(chrom_start + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...")
                        continue
 


                type_sequence = ''
                
                # Pulls out the mutation subtype for deletions
                if mut_type == 'Del': 
                    type_sequence = ref[1:]   
                    
                    type_length = len(type_sequence)
                    sequence = type_sequence
                    pos = start + type_length 
                    pos_rev = start 
                    while pos_rev - type_length > 0 and chrom_string[pos_rev-type_length:pos_rev] == type_sequence:
                        sequence = chrom_string[pos_rev-type_length:pos_rev] + sequence
                        pos_rev -= type_length
                    while pos + type_length < len(chrom_string) and chrom_string[pos:pos + type_length] == type_sequence:
                        sequence += chrom_string[pos:pos+type_length]
                        pos += type_length
                    
                    # Pulls out possible microhomology deletions
                    if type_length > 1 and len(sequence) == type_length:
                        forward_homology = ref[1:-1]
                        reverse_homology = ref[2:]

                        
                        for_hom = False
                        pos = start + type_length
                        for i in range (len(forward_homology), 0, -1):
                            if chrom_string[pos:pos+i] == forward_homology[:i]:
                                sequence += forward_homology[:i]
                                mut_type += '_Micro_for'
                                for_hom = True
                                break

                        if for_hom != True:
                            pos = start
                            for i in range (len(reverse_homology), 0, -1):
                                if chrom_string[pos-i:pos] == reverse_homology[-i:]:
                                    sequence = reverse_homology[-i:] + sequence
                                    mut_type += '_Micro_rev'
                                    break
                
                # Pulls out the mutation subtype for insertions
                elif mut_type == 'Ins':        
                    type_sequence = mut[1:]
                    type_length = len(type_sequence)
                    sequence = type_sequence

                    pos = start
                    pos_rev = start
                    while pos_rev - type_length > 0 and chrom_string[pos_rev-type_length:pos_rev] == type_sequence:
                        sequence = chrom_string[pos_rev-type_length:pos_rev] + sequence
                        pos_rev -= type_length
                    while pos + type_length < len(chrom_string) and chrom_string[pos:pos + type_length] == type_sequence:
                        sequence += chrom_string[pos:pos+type_length]
                        pos += type_length
                    
                    # Pulls possible microhomology for insertions
                    if type_length > 1 and len(sequence) == type_length:
                        forward_homology = mut[1:-1]
                        reverse_homology = mut[2:]
                        
                        for_hom = False
                        pos = start
                        for i in range (len(forward_homology), 0, -1):
                            if chrom_string[pos:pos+i] == forward_homology[:i]:
                                sequence += forward_homology[:i]
                                mut_type += '_Micro_for'
                                for_hom = True
                                break

                        if for_hom != True:
                            pos = start
                            for i in range (len(reverse_homology), 0, -1):
                                if chrom_string[pos-i:pos] == reverse_homology[-i:]:
                                    sequence = reverse_homology[-i:] + sequence
                                    mut_type += '_Micro_rev'
                                    break

                # Saves the sample name for later reference
                if sample not in indel_dict.keys():
                    indel_dict[sample] = {}
                if sample not in samples:
                    samples.append(sample)

                # Instantiates variables used to create the unique INDEL keys
                indel_key_1 = None
                indel_key_2 = None
                indel_key_3 = None
                indel_key_4 = None
                indel_key = 'blah'

                output_sequence = None

                # Creates the INDEL key for all deletions
                if mut_type[0:3] == 'Del': 
                    indel_key_2 = 'Del'

                    # Includes deletions of >1 bp
                    if len(ref)-1 > 1: 
                        key_1 = len(ref)-1
                        if key_1 < 5:
                            indel_key_1 = key_1
                        else:
                            indel_key_1 = 5

                        # Only regular deleletions
                        if mut_type == 'Del': 
                            indel_key_3 = 'R'
                            key_4 = int(len(sequence)/key_1 - 1)
                            if key_4 < 5:
                                indel_key_4 = key_4
                            else: 
                                indel_key_4 = 5

                        # Only for microhomologies
                        else:
                            indel_key_3 = 'M'
                            key_4 = len(sequence) - (len(ref)-1) 
                            if key_4 > 5:
                                indel_key_4 = 5
                            elif key_4 < 0:
                                print(lines)
                            else:
                                indel_key_4 = key_4
                
                    # For deletions of 1bp
                    else:
                        indel_key_1 = 1
                        key_4 = len(sequence) -1 
                        if key_4 > 5:
                            indel_key_4 = 5
                        else:
                            indel_key_4 = key_4
                        
                        if ref[1] == 'C' or ref[1] == 'G':
                            indel_key_3 = 'C'
                        
                        else:
                            indel_key_3 = 'T'

                            
                # Creates the INDEL key for all insertions
                elif mut_type[0:3] == 'Ins':
                    indel_key_2 = 'Ins'

                    #Includes insertions of >1bp
                    if len(mut)-1 > 1:
                        key_1 = len(mut)-1
                        if key_1<5:
                            indel_key_1 = key_1
                        else:
                            indel_key_1 = 5
                            
                        # Only regular insertions
                        if mut_type == 'Ins':
                            indel_key_3 = 'R'
                            key_4 = int(len(sequence)/key_1 - 1)
                            if key_4 < 5:
                                indel_key_4 = key_4
                            else:
                                indel_key_4 = 5
                        # Only for microhomologies
                        else:
                            indel_key_3 = 'M'
                            key_4 = len(sequence) - (len(mut)-1) 
                            if key_4 >= 5:
                                indel_key_4 = 5
                            elif key_4 < 0:
                                print(lines)
                            else:
                                indel_key_4 = key_4
                            
                    # Includes insertions of 1bp
                    else:
                        indel_key_1 = 1
                        key_4 = len(sequence)-1
                        if key_4 >= 5:
                            indel_key_4 = 5
                        else:
                            indel_key_4 = key_4
                            
                        if mut[1] == 'C' or mut[1] == 'G':
                            indel_key_3 = 'C'
                        else:
                            indel_key_3 = 'T'

                # Counts the number of "complex" mutations
                else:
                    non_matching += 1

                # Creates the final INDEl key and saves it into the data structure
                if limited_indel and indel_key_2 == 'Ins' and indel_key_3 == 'M':
                        indel_key = str(indel_key_1) + ':' + indel_key_2 + ':' + 'R' + ':' + '0'

                else:        
                    indel_key = str(indel_key_1) +':'+indel_key_2+':'+indel_key_3+':'+str(indel_key_4)

                if indel_key not in indel_dict[sample].keys():
                    indel_dict[sample][indel_key] = 1
                else:
                    indel_dict[sample][indel_key] += 1

                if exome:
                    exome_file.write(sample + '\t' + chrom + '\t' + str(start) + '\t' + mut_key + "\n")
            else:
                if 'non_matching' not in indel_dict[sample].keys():
                    indel_dict[sample]['non_matching'] = 1
                else:
                    indel_dict[sample]['non_matching'] += 1

    # Prints the total number of complex mutations
    print(non_matching)

    if exome:
        exome_file.close()
        os.system("sort -t $'\t' -k 2,2n -k 2,2 -k 3,3n " + exome_temp_file + " -o " + exome_temp_file)
        print("Beginning exome filtering. This may take a few moments...")
        mutation_dict, samples = exome_check (genome, exome_temp_file)
        os.system("rm " + exome_temp_file)


    # Calls the function to generate the final mutational matrix
    if functionFlag:
        return(pd.DataFrame.from_dict(indel_dict))
    else:
        matrix_generator_INDEL(output_matrix, samples, indel_types, indel_dict, project, exome, limited_indel, bed)

def exome_check (genome, exome_temp_file):
    '''
    Filters the variants for those present within the exome. 

    Parameters:
                 genome  -> name of the genome of interest (ex: GRCh37)
        exome_temp_file  -> The temporary file that contains all of the variants used for filtering

    Returns:
          mutation_dict  -> updated mutation dictionary for each sample for each mutation type post filtering
                samples  -> updated list of samples that still contain mutations post filtering

    ''' 

    # Instantiates the relevant variables/data structures
    base_cushion = 200
    mutation_dict = {}
    samples = []

    current_dir = os.getcwd()
    ref_dir = re.sub('\/scripts$', '', current_dir)


    initial = True
    udpate_chrom = False
    exome_file = ref_dir + "/references/chromosomes/exome/" + genome + "_exome.interval_list"


    with open(exome_temp_file) as f, open(exome_file) as exome:
        previous_chrom_ref = None
        previous_chrom_start = None
        previous_chrom_end = None

        chrom_ref = None
        start_ref = None
        end_ref = None

        read = True

        for lines in f:
            # Saves the relevant data for the current variant for later reference
            line = lines.strip().split()
            sample = line[0]
            chrom = line[1]
            start = int(line[2])
            mut_type = line[3]

            # Saves a value for the x and y chromosomes as a numeric reference
            if chrom == 'X':
                chrom_value = -1
            elif chrom == 'Y':
                chrom_value = 0
            else:
                chrom_value = int(chrom)


            if initial:
                chrom_start = chrom
                initial = False

            stop = False
            while not stop:
                if chrom == previous_chrom_ref:
                    if start >= previous_chrom_start - base_cushion and start <= previous_chrom_end + base_cushion:
                        if sample not in mutation_dict.keys():
                            samples.append(sample)
                            mutation_dict[sample] = {}
                        if mut_type not in mutation_dict[sample].keys():
                            mutation_dict[sample][mut_type] = 1
                        else:
                            mutation_dict[sample][mut_type] += 1
                        read = True
                        break

                if read:
                    lines2 = exome.readline()
                try:
                    if lines2[0] == "@":
                        continue
                except:
                    break
                else:
                    line2 = lines2.strip().split('\t')
                    chrom_ref = line2[0]
                    start_ref = int(line2[1])
                    end_ref = int(line2[2])

                    if chrom_ref == 'X':
                        ref_chrom_value = -1
                    elif chrom_ref == 'Y':
                        ref_chrom_value = 0
                    else:
                        ref_chrom_value = int(chrom_ref)

                    if chrom == chrom_ref:

                        if start > (start_ref - base_cushion and end_ref + base_cushion):
                            read = True
                            continue
                        elif start >= start_ref - base_cushion and start <= end_ref + base_cushion: 
                            if sample not in mutation_dict.keys():
                                samples.append(sample)
                                mutation_dict[sample] = {}
                            if mut_type not in mutation_dict[sample].keys():
                                mutation_dict[sample][mut_type] = 1
                            else:
                                mutation_dict[sample][mut_type] += 1
                            read = True
                            break
                        elif start < (start_ref - base_cushion):
                            read = False
                            break


                    else:
                        if chrom_value < ref_chrom_value:
                            read = False
                            break
                        elif chrom_value > ref_chrom_value:
                            read = True
                            continue



            chrom_start = chrom
            previous_chrom_ref = chrom_ref
            previous_chrom_start = start_ref
            previous_chrom_end = end_ref

    
    print("Exome filtering is complete. Proceeding with the final catalogue generation...")
    return(mutation_dict, samples)




def matrix_generator (context, output_matrix, project, samples, bias_sort, mutation_dict, types, exome, mut_types, bed):
    '''
    Writes the final mutational matrix given a dictionary of samples, mutation types, and counts

    Parameters:
                    context  -> desired context (ex. 96, 1536, 192, 3072)
              output_matrix  -> path where the final mutational matrix is stored
                    project  -> unique name given to the set of samples (ex. 'BRCA') 
                     samples -> a list of all sample names
                   bias_sort -> dictionary that provides the sorting order for the TSB matrices
             muatation_dict  -> dictionary with the counts for each mutation type for each sample
                      types  -> list of the mutation types for the given context 
                      exome  -> Boolean for whether the catalogue should be generated across the whole
                                genome or just the exome
                  mut_types  -> list with all possible mutation types for the given context
                        bed  -> parameter used to filter the mutations on a user-provided BED file



    Returns:
        None

    Output:
        Write the final mutational matrix for 96, 192, 1536, 3072 contexts

    '''

    mut_types_six = ['C>A','C>G','C>T','T>A','T>C','T>G']
    mut_types_twelve = ['T:C>A','T:C>G','T:C>T','T:T>A','T:T>C','T:T>G',
                        'U:C>A','U:C>G','U:C>T','U:T>A','U:T>C','U:T>G',
                        'B:C>A','B:C>G','B:C>T','B:T>A','B:T>C','B:T>G',
                        'N:C>A','N:C>G','N:C>T','N:T>A','N:T>C','N:T>G']

    mut_count_six = {'C>A':{},'C>G':{},'C>T':{},'T>A':{},'T>C':{},'T>G':{}}
    mut_count_twelve = {'T:C>A':{},'T:C>G':{},'T:C>T':{},'T:T>A':{},'T:T>C':{},'T:T>G':{},
                        'U:C>A':{},'U:C>G':{},'U:C>T':{},'U:T>A':{},'U:T>C':{},'U:T>G':{},
                        'B:C>A':{},'B:C>G':{},'B:C>T':{},'B:T>A':{},'B:T>C':{},'B:T>G':{},
                        'N:C>A':{},'N:C>G':{},'N:C>T':{},'N:T>A':{},'N:T>C':{},'N:T>G':{}}

    if bed:    
        file_prefix = project + "_BED" 
    else:
        file_prefix = project

    if exome:
        output_file_matrix = output_matrix + file_prefix + ".exome.mut" + context 
        if context == '96':
            output_file_matrix_six = output_matrix + file_prefix + ".exome.mut6"
        elif context == '192':
            output_file_matrix_twelve = output_matrix + file_prefix + ".exome.mut12"
    else:
        output_file_matrix = output_matrix + file_prefix + ".genome.mut" + context
        if context == '96':
            output_file_matrix_six = output_matrix + file_prefix + ".genome.mut6"
        elif context == '192':
            output_file_matrix_twelve = output_matrix + file_prefix + ".genome.mut12"



    with open (output_file_matrix, 'w') as out:

        # Prints all of the sample names into the first line of the file
        print ('MutationType\t', end='', flush=False, file=out)  
        for sample in samples:
            print (sample + '\t', end='', flush=False, file=out)
        print(file=out)

        if context == '192' or context == '3072':
            try:
                types = sorted(mut_types, key=lambda val: (bias_sort[val[0]], val[2:]))
            except:
                print(mut_types)
        else:
            types = mut_types

        # Prints the mutation count for each mutation type across every sample
        for mut_type in types:
            print (mut_type + '\t', end='', flush =False, file=out)
            for sample in samples:
                if context == '96':
                    subType = mut_type[2:5]
                    if sample not in mut_count_six[subType].keys():
                        if mut_type in mutation_dict[sample].keys():
                            mut_count_six[subType][sample] = mutation_dict[sample][mut_type]
                        else:
                            mut_count_six[subType][sample] = 0
                    else:
                        if mut_type in mutation_dict[sample].keys():
                            mut_count_six[subType][sample] += mutation_dict[sample][mut_type]
                        else:
                            mut_count_six[subType][sample] = 0

                elif context == '192':
                    subType = mut_type[0:2] + mut_type[4:7]
                    if sample not in mut_count_twelve[subType].keys():
                        if mut_type in mutation_dict[sample]:
                            mut_count_twelve[subType][sample] = mutation_dict[sample][mut_type]
                        else:
                            mut_count_twelve[subType][sample] = 0
                    else:
                        if mut_type in mutation_dict[sample]:
                            mut_count_twelve[subType][sample] += mutation_dict[sample][mut_type]
                        else:
                            mut_count_twelve[subType][sample] += 0

                if mut_type in mutation_dict[sample].keys():
                    print (str(mutation_dict[sample][mut_type]) + '\t', end='', file=out)
                else:
                    print ('0\t', end='', file=out)
            print(file=out)


    # sorts the 96 and 1536 matrices by mutation type
    if context == '96' or context == '1536':
        command1 = "head -1 " + output_file_matrix + " > " + output_matrix + "a.tmp;"
        command2 = "tail -n+2 " + output_file_matrix + " | sort -n >> " + output_matrix + "a.tmp;"

        
        os.system(command1)
        os.system(command2)
        os.system("cat " + output_matrix + "a.tmp > " + output_file_matrix)
        os.system("rm " + output_matrix + "a.tmp")
    
    # Generates the matrix for the 6 SNV mutation types
    if context == '96':
        with open(output_file_matrix_six, 'w') as out:
            print ('MutationType\t', end='', flush=False, file=out)  
            for sample in samples:
                print (sample + '\t', end='', flush=False, file=out)
            print(file=out)
            for types in mut_types_six:
                print (types + '\t', end='', flush =False, file=out)
                for sample in samples:
                    print(str(mut_count_six[types][sample]) + "\t",end='',file=out)
                print(file=out)
        print("Catalogue for 6 context is complete.")

    # Generates the matrix for the 12 SNV mutation types
    elif context == '192':
        with open(output_file_matrix_twelve, 'w') as out:
            print ('MutationType\t', end='', flush=False, file=out)  
            for sample in samples:
                print (sample + '\t', end='', flush=False, file=out)
            print(file=out)
            for types in mut_types_twelve:
                print (types + '\t', end='', flush =False, file=out)
                for sample in samples:
                    print(str(mut_count_twelve[types][sample]) + "\t",end='',file=out)
                print(file=out)
        print("Catalogue for 12 context is complete.")

def matrix_generator_INDEL (output_matrix, samples, indel_types, indel_dict, project, exome, limited_indel, bed):
    '''
    Writes the final mutational matrix for INDELS given a dictionary of samples, INDEL types, and counts

    Parameters:
              output_matrix  -> path where the final mutational matrix is stored
                    samples  -> a list of all sample names
                indel_types  -> list of the INDEL types 
                 indel_dict  -> dictionary with the counts for each INDEL type for each sample
                    project  -> unique name given to the set of samples (ex. 'BRCA') 
                      exome  -> Boolean for whether the catalogue should be generated across the whole
                                genome or just the exome
                        bed  -> parameter used to filter the mutations on a user-provided BED file

    Returns:
        None

    Output:
        Write the final mutational matrix for INDELS

    '''
    if bed:
        file_prefix = project + "_BED"
    else:
        file_prefix = project

    if exome:
        if limited_indel:
            output_file_matrix = output_matrix + file_prefix + ".exome.mutLimitedINDEL"
        else: 
            output_file_matrix = output_matrix + file_prefix + ".exome.mutINDEL" 
    else:
        if limited_indel:
            output_file_matrix = output_matrix + file_prefix + ".genome.mutLimitedINDEL"
        else:
            output_file_matrix = output_matrix + file_prefix + ".genome.mutINDEL"

    with open (output_file_matrix, 'w') as out:
        # Prints all of the sample names into the first line of the file
        print ('MutationType\t', end='', flush=False, file=out)  
        samples.sort()
        for sample in samples:
            print (sample + '\t', end='', flush=False, file=out)
        print(file=out)
        
        # Prints the mutation count for each INDEL type across every sample
        for indel in indel_types:
            print (indel + '\t', end='', flush =False, file=out)
            for sample in samples:
                if indel in indel_dict[sample].keys():
                    print (str(indel_dict[sample][indel]) + '\t', end='', file=out)
                else:
                    print ('0\t', end='', file=out)
            print(file=out)


def matrix_generator_DINUC (output_matrix, samples, mutation_types, dinucs, project, exome, bed):
    '''
    Writes the final mutational matrix for INDELS given a dictionary of samples, INDEL types, and counts

    Parameters:
              output_matrix  -> path where the final mutational matrix is stored
                    samples  -> a list of all sample names
             mutation_types  -> list of the DINUC types 
                     dinucs  -> dictionary with the counts for each DINUC type for each sample
                    project  -> unique name given to the set of samples (ex. 'BRCA') 
                      exome  -> Boolean for whether the catalogue should be generated across the whole
                                genome or just the exome
                        bed  -> parameter used to filter the mutations on a user-provided BED file

    Returns:
        None

    Output:
        Write the final mutational matrix for DINUCs

    '''
    if bed:
        file_prefix = project + "_BED"
    else:
        file_prefix = project

    if exome:
        output_file_matrix = output_matrix + file_prefix + ".exome.mutDINUC" 
    else:
        output_file_matrix = output_matrix + file_prefix + ".genome.mutDINUC"

    with open (output_file_matrix, 'w') as out:
        # Prints all of the sample names into the first line of the file
        print ('MutationType\t', end='', flush=False, file=out)  
        for sample in samples:
            print (sample + '\t', end='', flush=False, file=out)
        print(file=out)
        
        # Prints the mutation count for each INDEL type across every sample
        for dinuc in mutation_types:
            print (dinuc + '\t', end='', flush =False, file=out)
            for sample in samples:
                try:
                    if dinuc in dinucs[sample].keys():
                        print (str(dinucs[sample][dinuc]) + '\t', end='', file=out)
                    else:
                        print ('0\t', end='', file=out)
                except:
                    print('0\t', end='', file=out)
            print(file=out)
    


def main():
    start = time.time()

    ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
                  'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
                  'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
                  'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
                  'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
                  'NC_000087.7':'Y'}

    contexts = ['96', '1536', '192', '3072', 'DINUC']

    exome = False
    indel = False
    limited_indel = False
    functionFlag = False
    bed = False
    bed_file = None
    bed_ranges = None

    parser = argparse.ArgumentParser(description="Provide the necessary arguments to create the desired catalogue.")
    parser.add_argument("--project", "-p",help="Provide a unique name for your samples. (ex: BRCA)")
    parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
    parser.add_argument("-e", "--exome", help="Optional parameter instructs script to create the catalogues using only the exome regions. Whole genome context by default", action='store_true')
    parser.add_argument("-i", "--indel", help="Optional parameter instructs script to create the catalogue for limited INDELs", action='store_true')
    parser.add_argument("-ie", "--extended_indel", help="Optional parameter instructs script to create the catalogue for extended INDELs", action='store_true')
    parser.add_argument("-b", "--bed", nargs='?', help="Optional parameter instructs script to simulate on a given set of ranges (ex: exome). Whole genome context by default")
    args=parser.parse_args()
    project = args.project
    genome = args.genome
    #context = args.context
    
    if args.exome:
        exome = True

    if args.extended_indel:
        indel = True

    if args.indel:
        indel = True
        limited_indel = True

    if args.bed:
        bed = True
        bed_file = args.bed


        

    # Organizes all of the reference directories for later reference:
    current_dir = os.getcwd()
    ref_dir = re.sub('\/scripts$', '', current_dir)
    chrom_path = ref_dir + '/references/chromosomes/chrom_string/' + genome + "/"
    chromosome_TSB_path =  ref_dir + '/references/chromosomes/tsb/' + genome + "/"



    # Organizes all of the input and output directories:
    output_matrix = ref_dir + "/references/matrix/"
    vcf_path = ref_dir + '/references/vcf_files/' + project + "/"
    bed_path = ref_dir + '/references/vcf_files/BED/' + project + "/"

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

    if indel:
        contexts = ['INDEL']

    print("Starting catalogue generation...")

    if bed:
        bed_file_path = ref_dir + "/references/vcf_files/BED/" + project + "/" + bed_file
        bed_ranges = BED_filtering(bed_file_path)

    for context in contexts:
        # Single file:
        if context != 'DINUC' and context != 'INDEL':
            catalogue_generator_single (vcf_path, vcf_files, chrom_path, chromosome_TSB_path, project, output_matrix, context, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges)

        elif context == 'DINUC':
            catalogue_generator_DINUC_single (vcf_path, vcf_files, chrom_path, chromosome_TSB_path, project, output_matrix, exome, genome, ncbi_chrom, functionFlag, bed, bed_ranges)

        elif context == 'INDEL':
            catalogue_generator_INDEL_single (vcf_path, vcf_files, chrom_path, project, output_matrix, exome, genome, ncbi_chrom, limited_indel, functionFlag, bed, bed_ranges)

        print("Catalogue for " + context + " context is complete.")
    os.system("rm -r " + vcf_path)
    end = time.time()
    print("Job took ",str(end-start), " seconds.")


if __name__ == '__main__':
    main()