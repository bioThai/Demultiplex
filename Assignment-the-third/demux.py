#!/bin/python
import argparse
import gzip     #needed to open/read gzipped files and write to gzipped files

def get_args():
    '''Defines/sets possible command line arguments for script'''
    parser = argparse.ArgumentParser("A program to demultiplex FASTQ data")
    parser.add_argument("-r", nargs="+", help="Specifies input FASTQ read file name(s). Note: Only the first 2 files specified will be parsed.", type=str, required=True)
    parser.add_argument("-i", nargs="+", help="Specifies input FASTQ index file name(s). Note: Only the first 2 files specified will be parsed.", type=str, required=True)
    parser.add_argument("-t", help="Specifies text file containing the known reference indexes that the indexes in the FASTQ index files will be compared to.", type=str, required=True)
    parser.add_argument("-q", help="Specifies quality score mimimum value to use as cutoff for index file qscores (default=30).", type=int, default=30)
    return parser.parse_args()

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    phred_value = ord(letter) - 33
    return phred_value

def rev_comp(seq: str) -> str:
    '''Takes an input sequence string (written 5' to 3' from left to right) and returns the reverse complement as a string (written from 5' to 3' from left to right).'''
    # local variables
    comp_seq: str = ""      #holds the complementary sequence to the input sequence
    rev_comp_seq: str = ""  #holds the reverse complementary sequence to the input sequence

    for char in seq.upper():
        if char == "A":
            comp_seq += "T"
        elif char == "T":
            comp_seq += "A"
        elif char == "G":
            comp_seq += "C"
        elif char == "C":
            comp_seq += "G"
    rev_comp_seq = comp_seq[::-1]   #reverse the comp_seq string: slice the comp_seq string from beginning to end, but step the increments by -1 (go backwards through the string)
    return rev_comp_seq

def get_ref_indexes(ref_indexes_file: str) -> dict:
    '''Takes a text file of known reference indexes and puts the index names (A1, A2, B1, C1, etc) and index sequences in a dictionary. Returns the dictionary.'''
    # local variables
    index_name: str = ""
    index_seq: str = ""
    ref_indexes_dict: dict = {} #keys: index sequences, values: index names

    with open(ref_indexes_file, "r") as fh:
        fh.readline()   #read the first line in file but don't store it anywhere. Just want to move the file-reading pointer past the first line (First line is just a header line, not needed)
        for line in fh:
            line = line.strip()
            index_name = line.split()[3]
            index_seq = line.split()[4]
            if index_seq not in ref_indexes_dict:
                ref_indexes_dict[index_seq] = index_name
    return ref_indexes_dict

def get_output_files_dict(read1_file: str, read2_file: str, ref_indexes_dict: dict) -> dict:
    '''Create a dictionary to keep track of 52 output filenames: 26 FASTQ files for each of the two input biological read FASTQ files (total of 52 FASTQ files). For each input read file, there are 24 output FASTQ files, with each file containing all the correctly-indexed reads for a specific index-pair (correct index sequences are listed in indexes.txt file on Talapas). In addition, for each input read file, there is also a 25th output FASTQ file for all the reads with hopped indexes, and a 26th output FASTQ file for all the reads with indexes that don't match the indexes in the indexes.txt file and/or indexes with low-quality scores.'''
    # local variables
    output_files_dict: dict = {}    #keys: reference index sequence, values: list [read1 output file for index, read2 output file for index]
    read1_outfile: str = ""
    read2_outfile: str = ""

    #split the absolute filepaths stored in read1_file and read2_file based on "/", and extract only the name of the file itself (relative to current directory)
    read1_file = read1_file.split("/")[-1]
    read2_file = read2_file.split("/")[-1]

    #create output filenames for each index, for each input read file (48 filenames)
    for ref_index_seq, ref_index_name in list(ref_indexes_dict.items()):
        read1_outfile = ref_index_name + "_" + read1_file
        read2_outfile = ref_index_name + "_" + read2_file
        if ref_index_seq not in output_files_dict:
            output_files_dict[ref_index_seq] = []
        output_files_dict[ref_index_seq].extend([read1_outfile, read2_outfile])

    #create 2 filenames for swapped output files
    read1_outfile = "swapped_" + read1_file
    read2_outfile = "swapped_" + read2_file
    if "swapped" not in output_files_dict:
        output_files_dict["swapped"] = []
    output_files_dict["swapped"].extend([read1_outfile, read2_outfile])  

    #create 2 filenames for unknown/low-qscore output files
    read1_outfile = "unknown_lowQ_" + read1_file
    read2_outfile = "unknown_lowQ_" + read2_file
    if "unknown_lowQ" not in output_files_dict:
        output_files_dict["unknown_lowQ"] = []
    output_files_dict["unknown_lowQ"].extend([read1_outfile, read2_outfile])  
    
    return output_files_dict

def get_record(input_filehandlers_list: list) -> list:
    '''Read 4 lines (1 record) from each of the 4 input filehandlers in a given list, store lines in 2d list, and return the reference to the 2d list'''
    # local variables
    row: int = 0
    col: int = 0
    line: str = ""
    record_list: list = [["" for i in range(4)] for j in range(4)]     
        #2D list of strings (4 rows, 4 cols) to temporarily hold a record from each of the 4 input files  
        #record_list[0]: holds record for input biological read file 1
        #record_list[1]: holds record for input biological read file 2
        #record_list[2]: holds record for input index read file 1
        #record_list[3]: holds record for input index read file 2
        #for each record from an input file, 
            #record_list[i][0]: holds header line
            #record_list[i][1]: holds sequence line
            #record_list[i][2]: holds "+" line
            #record_list[i][3]: holds quality score line

    row = 0
    while row < 4:
        #print(input_filehandlers_list[row]) 
        col = 0
        while col < 4:
            line = input_filehandlers_list[row].readline()
            line = line.strip("\n")
            #print(line)
            record_list[row][col] = line
            col += 1
        row += 1
    #print(record_list)    
    return record_list

def parse_input_files(read1_file: str, read2_file: str, index1_file: str, index2_file: str, ref_indexes_dict: dict, output_files_dict: dict, qscore_cutoff: int):
    '''Takes 2 FASTQ read files and 2 FASTQ index files and demultiplexes them based on indexes in a given reference index file and based on a given quality score cutoff.'''
    # local variables
    is_index_qscore_below_cutoff: bool = False

    temp_records_list: list = []
    #2D list of strings (4 rows, 4 cols) to temporarily hold a record from each of the 4 input files  
        #record_list[0]: holds record for input biological read file 1
        #record_list[1]: holds record for input biological read file 2
        #record_list[2]: holds record for input index read file 1
        #record_list[3]: holds record for input index read file 2
        #for each record from an input file, 
            #record_list[i][0]: holds header line
            #record_list[i][1]: holds sequence line
            #record_list[i][2]: holds "+" line
            #record_list[i][3]: holds quality score line
    input_fh_list: list = []    #holds the filehandlers for each of the 4 input FASTQ files
    output_fh_dict: dict = {}   #keys: holds index sequences of each reference index from output_files_dict. For swapped and unknown/low-qscore files, key holds "swapped" or "unknown_lowQ".
                                #values: In a list w/2 elements (element 0 for read1 filehandlers, element 1 for read2 filehandlers), holds all the filehandlers of output files that get opened
    record_counters_dict: dict = {} #keys: holds index sequences of each reference index from output_files_dict. For swapped and unknown/low-qscore files, key holds "swapped" or "unknown_lowQ".
                                    #values: Holds an integer counter of number of read-pairs with each index sequence. 
    record_counts_filename: str = "demux_final_report.tsv"
    sum_of_reads: int = 0
    percent_of_reads: float = 0

    read_record_count: int = 0  #holds count of number of records read from input read files

    #initialize values of output_fh_dict
    for key in output_files_dict:
        if key not in output_fh_dict:
            output_fh_dict[key] = []  

    #initialize values of record_counters_dict
    for key in output_files_dict:
        if key not in record_counters_dict:
            record_counters_dict[key] = 0
    
    #open all gzipped FASTQ read and index files
    read1_fh = gzip.open(read1_file, "rt")
    read2_fh = gzip.open(read2_file, "rt")
    index1_fh = gzip.open(index1_file, "rt")
    index2_fh = gzip.open(index2_file, "rt")
    input_fh_list = [read1_fh, read2_fh, index1_fh, index2_fh]
    
    #open all output files for writing
    for ref_index_seq in output_files_dict:
        for output_file in output_files_dict[ref_index_seq]:
            #fh = gzip.open(output_file, 'wt')
            fh = open(output_file, "w") #much faster than outputting/writing into gzipped files
            output_fh_dict[ref_index_seq].append(fh)

    #demultiplex each biological read record in each input file
    while True:    
        temp_records_list = get_record(input_fh_list)
        is_index_qscore_below_cutoff = False
        
        #if EOF has been reached:
        if temp_records_list[0][0] == "":
            break
    
        read_record_count += 1
        if read_record_count % 1000000 == 0:
            print("Records read: ", read_record_count)
            print("Lines read: ", 4 * read_record_count)
            print()

        #modify header lines of input biological read/record
        temp_records_list[0][0] += " " + temp_records_list[2][1] + "-" + temp_records_list[3][1]
        temp_records_list[1][0] += " " + temp_records_list[2][1] + "-" + temp_records_list[3][1]
        
        #check if any qscore in any index records is below qscore threshold
        for char in temp_records_list[2][3]:
            if convert_phred(char) < qscore_cutoff:
                is_index_qscore_below_cutoff = True
                break
        for char in temp_records_list[3][3]:
            if convert_phred(char) < qscore_cutoff:
                is_index_qscore_below_cutoff = True
                break

        if (
            ("N" in temp_records_list[2][1] or "N" in temp_records_list[3][1])
            or (temp_records_list[2][1] not in ref_indexes_dict or rev_comp(temp_records_list[3][1]) not in ref_indexes_dict)
            or (is_index_qscore_below_cutoff == True)
            ):
            
            #write each biological read/record to corresponding unknown_lowQ output file
            for i in range(2):
                for j in range(4):
                    output_fh_dict["unknown_lowQ"][i].write(temp_records_list[i][j] + "\n")
                    #print(temp_records_list[i][j])

            #increment counter
            record_counters_dict["unknown_lowQ"] += 1

        else:
            if temp_records_list[2][1] == rev_comp(temp_records_list[3][1]):
                #Write each biological read/record to corresponding output bucket file
                for i in range(2):
                    for j in range(4):
                        output_fh_dict[temp_records_list[2][1]][i].write(temp_records_list[i][j] + "\n")
                        #print(temp_records_list[i][j])

                #Increment counter of records in for corresponding bucket (<bucket>_counter += 1)
                record_counters_dict[temp_records_list[2][1]] += 1
            else:
                # Write record to appropriate "swapped" output file
                for i in range(2):
                    for j in range(4):
                        output_fh_dict["swapped"][i].write(temp_records_list[i][j] + "\n")
                        #print(temp_records_list[i][j])
                        
                # Increment counter of records for swapped reads
                record_counters_dict["swapped"] += 1

    #close all FASTQ read and index files, and all output files
    read1_fh.close()
    read2_fh.close()
    index1_fh.close()
    index2_fh.close()
    for ref_index_seq in output_fh_dict:
        for output_fh in output_fh_dict[ref_index_seq]:
            output_fh.close()

    #write final demux stats to output TSV file
    with open(record_counts_filename, "w") as output_stats_fh:
        output_stats_fh.write("Bucket" + "\t" + "Number_of_Read_Pairs" + "\t" + "Percentage_of_Reads" + "\n")
        sum_of_reads = sum(list(record_counters_dict.values()))
        for index_seq in record_counters_dict:
            percent_of_reads = (record_counters_dict[index_seq] / sum_of_reads) * 100
            output_stats_fh.write(index_seq + "\t" + str(record_counters_dict[index_seq]) + "\t" + str(percent_of_reads) + "\n")

    #plot distribution/percentages of reads from each bucket


def main():
    '''Main function, drives the order of execution for script'''
    #local variables
    args = get_args()
    readfile_list: list = args.r
    indexfile_list: list = args.i
    ref_indexes_file: str = args.t
    qscore_cutoff: int = args.q
    ref_indexes_dict: dict = {}     #keys: reference index sequences, values: reference index names
    output_files_dict: dict = {}    #keys: name of output FASTQ files, values: number of FASTQ records in file

    ref_indexes_dict = get_ref_indexes(ref_indexes_file)
    #print(ref_indexes_dict)
    output_files_dict = get_output_files_dict(readfile_list[0], readfile_list[1], ref_indexes_dict)
    #print(output_files_dict)
    #print(len(list(output_files_dict.keys())))
    parse_input_files(readfile_list[0], readfile_list[1], indexfile_list[0], indexfile_list[1], ref_indexes_dict, output_files_dict, qscore_cutoff)


if __name__ == "__main__":
    main()
    