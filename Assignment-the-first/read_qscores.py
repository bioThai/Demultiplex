#!/bin/python
import argparse
import gzip     #needed to read g-zipped files
import matplotlib.pyplot as plt
import Bioinfo

def get_args():
    '''Defines/sets possible command line arguments for script'''
    parser = argparse.ArgumentParser("A program to process FASTQ data, get average quality score for each read position, and plot distribution")
    parser.add_argument("-f", nargs="+", help="specifies input FASTQ filename(s)", type=str, required=True)
    parser.add_argument("-o", help="specifies output file prefix", type=str, required=True)
    return parser.parse_args()

def parse_file(in_filename: str, out_file_prefix: str):
    '''Open gzipped FASTQ file, put all qscores in file into arrays to track qscores for each nucleotide position.'''
    # local variables
    line_num: int = 0   #keep track of what line number of input file you're on
    read_num: int = 0   #keep track of number of records/reads in input file
    nuc_pos: int = 0    #keep track of nucleotide position
    converted_phred: int = 0
    sum_qscores_list: list = []     #keep track of the running total of qscores for each nucleotide position
    mean_qscores_list: list = []    #holds the mean qscore for each nucleotide position
    output_fname: str = "means_" + out_file_prefix + ".tsv"
    plot_name: str = "plot_" + out_file_prefix + ".png"

    with gzip.open(in_filename, 'rt') as input_file, open(output_fname, 'w') as output_file:
        for line in input_file:
            line_num += 1
            line = line.strip()
            if line_num % 4 == 0:   #grab the qscore line in filename
                nuc_pos = 0
                for char in line:
                    converted_phred = Bioinfo.convert_phred(char)
                    #if the sum_qscores_list doesn't yet have an element for each read position in the file:
                    if len(sum_qscores_list) != len(line):
                        sum_qscores_list.append(converted_phred)
                    else:
                        sum_qscores_list[nuc_pos] += converted_phred
                    nuc_pos += 1
                read_num += 1   #once whole qscore line is read, increment the read number    
        mean_qscores_list = [sum_qscore / read_num for sum_qscore in sum_qscores_list]    #divide each element in sum_qscores_list by the number of reads in the file to get the avg qscore for each nucleotide position. Put each of these avg scores into a mean_qscore_list.    
    
        #write the mean qscore for each nucleotide position into an output file
        output_file.write("Position\tMean\n")
        for i in range(len(mean_qscores_list)):
            output_file.write(str(i) + "\t" + str(mean_qscores_list[i]) + "\n")
    
    #plot mean qscore distribution
    plt.bar(range(len(mean_qscores_list)), mean_qscores_list)
    plt.title("Mean Quality Score for all Nucleotide Positions")
    plt.xlabel("Nucleotide Position")
    plt.ylabel("Mean Quality Score")
    plt.savefig(plot_name)

def main():
    '''Main function, drives the order of execution for script'''
    # local variables
    args = get_args()
    output_file_prefix = args.o

    for input_filename in args.f:
        parse_file(input_filename, output_file_prefix)

if __name__ == "__main__":
    main()
