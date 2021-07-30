# Assignment the First

## Part 1
1. Be sure to upload your Python script. (See ```read_qscores.py```)

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
        
       Distribution for 1294_S1_L008_R1_001.fastq.gz:
       ![plot1](https://github.com/bioThai/Demultiplex/blob/55257c3ae4d0fc6eb938940225514aca83ed6d32/Assignment-the-first/plot_1294_S1_L008_R1_001.fastq.gz.png)
       
       Distribution for 1294_S1_L008_R2_001.fastq.gz:
       ![plot2](https://github.com/bioThai/Demultiplex/blob/55257c3ae4d0fc6eb938940225514aca83ed6d32/Assignment-the-first/plot_1294_S1_L008_R2_001.fastq.gz.png)
       
       Distribution for 1294_S1_L008_R3_001.fastq.gz:
       ![plot3](https://github.com/bioThai/Demultiplex/blob/55257c3ae4d0fc6eb938940225514aca83ed6d32/Assignment-the-first/plot_1294_S1_L008_R3_001.fastq.gz.png)
       
       Distribution for 1294_S1_L008_R4_001.fastq.gz:
       ![plot4](https://github.com/bioThai/Demultiplex/blob/55257c3ae4d0fc6eb938940225514aca83ed6d32/Assignment-the-first/plot_1294_S1_L008_R4_001.fastq.gz.png)
    
    
    2. What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.
    
       A good quality score cutoff for index reads might be around 35 since index reads need to be very certain and you don't want any undetermined (N) base calls that would prevent you from accurately matching indexes and identifying corresponding biological reads. From the two index read files above, the first one or two base calls for many of the index sequences are Ns, and these nucleotide positions still have average Qscores above 30. The rest of the nucleotide positions in the index reads have average Qscores of at least 35. So, to be safe, a Qscore minimum around 35 might be a good cutoff for index reads.
       
       A good quality score cutoff for biological read pairs might be a little less than for index reads, so perhaps a qscore of around 30 might work. If you set the Qscore for the biological read pairs too high, you might have less data for downstream analysis.
    
    
    3. How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)
   
       Command (on Talapas/Linux): 
       ```
       zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n "2~4p" | grep -c "N"
       ```
       Number of indexes with undetermined base calls: 7304664
       
       
    
## Part 2
1. Define the problem

   There are biological reads in the provided FASTQ files which have incorrect indexes due to index hopping/swapping. Other reads may have missing/unknown/low-quality-score indexes. We need to demultiplex all reads, figure out which reads have hopped indexes, which reads have unknown/low quality indexes, and separate these reads from each correctly-indexed read. In addition, for each biological read file, each correctly-indexed read needs to be separated into different bins/files based on their indexes. As each read is processed, the sequence of the index pairs corresponsiding to each read need to be appended to the end of the header line in the FASTQ file.

2. Describe output

   The output of the demultiplexing algorithm is 26 FASTQ files for each of the two input biological read FASTQ files (total of 52 FASTQ files). For each input read file, there are 24 output FASTQ files, with each file containing all the correctly-indexed reads for a specific index-pair (correct index sequences are listed in ```indexes.txt``` file on Talapas). In addition, for each input read file, there is also a 25th output FASTQ file for all the reads with hopped indexes, and a 26th output FASTQ file for all the reads with indexes that don't match the indexes in the indexes.txt file and/or indexes with low-quality scores.
   
   The algorithm should also output useful data such as the number of records/reads in each of the above output files.


3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

   See uploaded files in repo.


4. Pseudocode

   ```
   Parse indexes.txt file to get the correct index sequences you're supposed to see in the FASTQ read files.
        Put each index into a dictionary: 
            key = the "name" (A1, B1, C9, etc) of each index, value = actual index sequence
        For each input biological read file (R1, R4 files), 
            create/open files to hold the FASTQ reads for each index in the dictionary
            create/open a file to hold the FASTQ reads w/swapped indexes
            create/open a file to hold the FASTQ reads for 
           
   ```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
