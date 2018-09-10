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
import time
import re
import argparse

start_time = time.time()

def save_tsb (chromosome_string_path, transcript_path, output_path):
    '''
    Creates binary files that contain the transcriptional information at a given 
    base. The transcriptional information entails Transcribed, Untranscribed, 
    Bi-directionally transcribed, and Non-transcribed. These files are required to
    simulate transcriptional strand bias.

    Input:
        chromosome_string_path -> path to chromosomes saved in individual text files
        transcript_path -> path to transcript data. The following should be used as the format
                           (additional information may be saved, however it will not
                           be used when creating the files):

        Gene stable ID  Transcript stable ID    Chromosome/scaffold name    Strand  Transcript start (bp)   Transcript end (bp)

        output_path -> path where the binary transcript files are saved.

    Output:
        Binary files that contain the transciptional information at each base for
        each chromosome. The files are saved as "[chrom]_192.txt"

    '''

    # Retrieves the transcript file 
    transcript_files = os.listdir(transcript_path)

    if len(transcript_files) < 2:
        # Instantiates required flag variables
        initial_flag = True
        initial_chrom = None
        for file in transcript_files:
            file_name = file.split("_")
            chrom = file_name[0]

            # Skips hidden files
            if chrom == '.DS':
                pass
            else:

                # Sorts the transcript file by chromosome and transcript range
                os.system("sort -t $'\t' -k 3,3n -k 3,3 -k 5,5 -k 6,6 " + transcript_path+file + " -o " + transcript_path+file)
                with open (transcript_path + file) as f:
                    next(f)

                    # Saves all transcriptional information for a chromosome into a unique file
                    for lines in f:
                        line = lines.strip()
                        line_split = line.split()
                        chrom = line_split[2]

                        if initial_flag:
                            initial_chrom = chrom
                            initial_flag = False
                            out = open (transcript_path+chrom + "_transcripts.txt", 'w')

                        if chrom == initial_chrom:
                            print(line, file=out)

                        else:
                            out.close()
                            initial_chrom = chrom
                            out = open(transcript_path+chrom+"_transcripts.txt",'w')
                    out.close()

            # Removes the initial transcript file
            os.system("rm " + transcript_path + file)

    # Retrieves the chromosome based list of transcript files 
    transcript_files = os.listdir(transcript_path)

    print("Creating the transcriptional reference files now. This may take awhile...")


    for files in transcript_files:
        file_name = files.split("_")
        chrom = file_name[0]
        if chrom == '.DS':
            pass
        else:
            try:
                if os.path.exists(output_path + chrom + "_192.txt"):
                    continue

                else:
                    # Create a binary file for the output.
                    outFile = open (output_path + chrom + "_192.txt", "wb")


                    # Instantiates all of the required parameters.
                    transcripts = None
                    bi_start = 0
                    bi_end= None
                    I_chrom = None
                    I_T_start = None
                    I_T_end = None
                    I_strand = None
                    location = 0
                    chrom_length = None

                    # Save the length for the current chromosome.
                    with open (chromosome_string_path + chrom + ".txt") as f:
                        chrom_length = len(f.readline())

                    # Establishes a unique binary number for each transcript type.
                    Non_transcribed = '00000000'
                    Transcribed = '00000001'
                    Untranscribed = '00000010'
                    BiDirectional = '00000011'

                    # Organizes the binary/integer formats for the data structures.
                    nN = int(Non_transcribed, 2)
                    nT = int(Transcribed, 2)
                    nU = int(Untranscribed, 2)
                    nB = int(BiDirectional, 2)
                    N = bytes([nN])
                    T = bytes([nT])
                    U = bytes([nU])
                    B = bytes([nB])
                    stringTest = ''


                    l = 0
                    line_count = 1
                    pointer = 0
                    with open(transcript_path + files, 'r') as f:
                        all_lines2 = f.readlines()
                        all_lines = all_lines2[1:]
                        for lines in all_lines:
                            
                            # Handles the first line separately. 
                            if pointer == 0:
                                first_line = lines.split()
                                location = int(first_line[4])
                                I_chrom = first_line[2][3:]
                                I_T_start = int(first_line[4])
                                I_T_end = int(first_line[5])
                                I_strand = first_line[3]
                                
                                # Saves Non-transcribed data up until
                                # the first transcript.
                                for i in range(0, location, 1):
                                    outFile.write(N)
                                    l += 1
                                pointer = 1

                            # Handles all other lines after the first.
                            else:
                                first_line = lines.split()
                                I_chrom = first_line[2][3:]
                                I_T_start = int(first_line[4])
                                I_T_end = int(first_line[5])
                                I_strand = first_line[3]

                                # Saves Non-transcribed data up until the 
                                # next transcript.
                                if I_T_start > location:
                                    for i in range(location, I_T_start, 1):
                                        outFile.write(N)
                                        l += 1
                                    location = I_T_start
                                
                            # Reads the subsequent line to look for bi-directional
                            # transcription sites.
                            for line2 in all_lines[line_count:]:
                                next_line = line2.split()
                                c_chrom = next_line[2][3:]
                                c_start = int(next_line[4])
                                c_end= int(next_line[5])
                                c_strand = next_line[3]
                                
                                # Breaks to the next line if the two transcripts 
                                # don't overlap.
                                if c_start > I_T_end:
                                    break

                                # Checks if the transcripts are on opposite strands
                                # if they overlap.    
                                else:
                                    if c_strand != I_strand:
                                        bi_start = c_start
                                        if c_end > I_T_end:
                                            bi_end = I_T_end
                                        else:
                                            bi_end = c_end

                                    # Saves Un/Transcribed data up until the start of
                                    # the bi-directional transcription.
                                    if bi_start != 0 and bi_start > location:
                                        for i in range (location, bi_start, 1):
                                            if I_strand == '1':
                                                outFile.write(T)
                                                l += 1
                                            else:
                                                outFile.write(U)
                                                l += 1
                                        location = bi_start

                                        # Saves Bi-directional data for the length
                                        # of the opposing transcripts overlap
                                        for i in range (location, bi_end, 1):
                                            outFile.write(B) 
                                            l += 1
                                        location = bi_end
                                    bi_start = 0
                                    bi_end = 0

                            # Saves Un/Transcribed data up to the end of the current
                            # transcript if data has not already been saved for these bases.
                            if I_T_end > location:
                                for i in range(location, I_T_end, 1):
                                    if I_strand == '1':
                                        outFile.write(T)
                                        l += 1
                                    else:
                                        outFile.write(U)
                                        l +=1
                                location = I_T_end
                            line_count += 1

                    # Save Non-transcribed data for the remainder of the chromosome after
                    # the final transcript is analyzed         
                    if location < chrom_length:
                        for i in range(location, chrom_length,1):
                            outFile.write(N) 
                            l += 1        
                        
                    outFile.close()
                    print("The transcript reference file for chromosome " +chrom + " has been created.")

            except:
                print(chrom + " has been skipped. This chromosome is not supported.")

    end_time = time.time()
    print("Transcript files created.\n Job took: ", end_time-start_time, " seconds")

def main ():
    parser = argparse.ArgumentParser(description="Provide the necessary arguments to create the transcriptional strand bias reference strings.")
    parser.add_argument("--genome", "-g",help="Provide a reference genome. (ex: GRCh37, GRCh38, mm10)")
    args=parser.parse_args()
    genome = args.genome

    script_dir = os.getcwd()
    ref_dir = re.sub('\/scripts$', '', script_dir)

    chromosome_string_path = ref_dir + "/references/chromosomes/chrom_string/" + genome + "/"
    transcript_path = ref_dir + "/references/chromosomes/transcripts/" + genome + "/"
    output_path = ref_dir + "/references/chromosomes/tsb/" + genome + "/"
    if os.path.exists(output_path) == False:
        os.makedirs(output_path)

    save_tsb(chromosome_string_path, transcript_path, output_path)

if __name__ == '__main__':
    main()


