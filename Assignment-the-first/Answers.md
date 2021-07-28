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

    There are reads in the FASTQ files with incorrect indexes due to index hopping/swapping. We need to figure out which reads have hopped indexes, and separate the reads which have hopped


2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
