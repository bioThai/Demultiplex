#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --output=demux_wrapper%j.out
#SBATCH --error=demux_wrapper%j.err

conda activate bgmp_py39

script_file="/projects/bgmp/tnguye14/bioinfo/Bi622/Demultiplex/Assignment-the-third/demux.py"
read1_file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
read2_file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
index1_file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
index2_file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
ref_indexes_file="/projects/bgmp/shared/2017_sequencing/indexes.txt"

/usr/bin/time -v python $script_file -r $read1_file $read2_file -i $index1_file $index2_file -t $ref_indexes_file

exit
