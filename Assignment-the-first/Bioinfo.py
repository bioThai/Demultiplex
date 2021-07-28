####### Handy constants: #######
DNA_BASES: set = set('ATGCatcg')
RNA_BASES: set = set('AUGCaucg')


####### Functions: #######
def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    phred_value = ord(letter) - 33
    return(phred_value)

def qual_score(phred_score: str) -> float:
    """Takes an unmodified phred score string of letters, iterates through each letter, converts each letter to a numerical phred value, and calculates the average quality score of the whole phred string."""
    index = 0
    converted_score_sum = 0
    for letter in phred_score:
        converted_score = convert_phred(phred_score[index])
        converted_score_sum += converted_score
        index += 1
    avg_score = converted_score_sum / len(phred_score)
    return avg_score

def validate_DNA_seq(seq: str) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("T") + seq.count("G")

def gc_content(DNA: str) -> float:
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_DNA_seq(DNA), "String contains invalid characters"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

#To use pythag function, must import some Python modules first:
import math
from typing import Union
def pythag(a: Union[int, float], b: Union[int, float]) -> float:
    '''Computes the length of the hypotenuse of a right triangle given the lengths
    of the other two sides.'''
    c = math.sqrt(a**2 + b**2)
    return c

#Protein FASTA filter (must import re first):
import re
def filter_longest_recs(filename: str):
    '''Reads through a fasta file of proteins, extracts the longest protein record for each gene, and outputs these records into another fasta file.'''
    #local variables
    output_fname: str = "filtered_" + filename
    prot_stable_ID: str = ""    #protein-stable ID
    gene_stable_ID: str = ""
    gene_name: str = ""
    curr_seq: str = ""      #holds current sequence line being read
    gene_dict1: dict = {}   #key:value -> geneID:new longest protein sequence
    gene_dict2: dict = {}   #key:value -> geneID:[protein ID, gene_name]
    output_header: str = ""
    output_seq: str = ""

    with open(filename, "r") as input_fh, open(output_fname, "w") as output_fh:
        for line in input_fh:
            line = line.strip()
            if line.startswith(">"):    #if header line is being read
                curr_seq = ""           #for each header line read, reset the curr_seq back to "" for the next seq read

                prot_stable_ID = re.search("(^>[a-zA-Z]+[Pp][0-9]+)([.][0-9]+)", line).group(1)
                gene_stable_ID = re.search("(gene:)([a-zA-Z]+[Gg][0-9]+)([.][0-9]+)", line).group(2)
                
                search_res = re.search("(gene_symbol:)([a-zA-Z0-9_.:-]+)", line)
                if search_res != None:
                    gene_name = search_res.group(2)
                else:
                    gene_name = ""
                
                if gene_stable_ID not in gene_dict1:
                    gene_dict1[gene_stable_ID] = ""
                    gene_dict2[gene_stable_ID] = ["",""]
            else:
                curr_seq += line
                #if length of protein seq that was just read > length of protein sequence previously read (that also belonged to same gene_ID from the most recently-read header line)
                #If protein seq spans multiple lines, this code still works. 
                # The curr_seq string just grows with each seq line read, and the protein-stable ID, gene name, and gene ID have not changed since the most recent header line was read, so it's ok to keep updating the dictionaries w/these values .
                if len(curr_seq) > len(gene_dict1[gene_stable_ID]):
                    gene_dict1[gene_stable_ID] = curr_seq
                    gene_dict2[gene_stable_ID] = [prot_stable_ID, gene_name]
        
        #write filtered fasta records to output file
        for genID in gene_dict2:
            output_header = gene_dict2[genID][0] + " " + genID + " " + gene_dict2[genID][1]
            output_seq = gene_dict1[genID] 
            output_fh.write(output_header + "\n" + output_seq + "\n")

# convert FASTA files with multiple sequence lines per record into FASTA files with each record's sequence on just 1 line
#### Note: This function has not yet been tested. It needs to be tested/verified before you use it for the first time!
def oneline_fasta(filename: str):
    '''Reads through a FASTA file and concatenates each multi-line sequence into a single sequence line for each record.'''
    #local variables
    output_fname: str = "oneline_" + filename
    curr_header: str = ""
    curr_seq: str = ""      #holds current sequence line being read
    record_num: int = 0
    header_lines_list: list = []
    seq_lines_list: list = []

    with open(filename, "r") as input_fh, open(output_fname, "w") as output_fh:
        for line in input_fh:
            line = line.strip()
            
            if line.startswith(">"):    #if header line is being read
                record_num += 1
                curr_seq = ""           #for each header line read, reset the curr_seq back to "" for the next seq read
                curr_header = line
                header_lines_list.append(curr_header)
            else:
                curr_seq += line

            if curr_seq != "" and len(seq_lines_list) == record_num - 1:
                seq_lines_list.append(curr_seq)
            elif curr_seq != "" and len(seq_lines_list) == record_num:
                seq_lines_list[record_num - 1] = curr_seq    

        #write FASTA headers and concatenated sequences to output file
        for i in range(len(header_lines_list)):
            output_fh.write(header_lines_list[i] + "\n" + seq_lines_list[i] + "\n")



if __name__ == "__main__":
    print("final_Bioinfo.py is being run directly.")
else:
    print("final_Bioinfo.py has been imported for module use.")