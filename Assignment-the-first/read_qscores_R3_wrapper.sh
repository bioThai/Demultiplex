#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --time=16:00:00
#SBATCH --output=r3_wrapper%j.out
#SBATCH --error=r3_wrapper%j.err

conda activate bgmp_py39
script_file="/projects/bgmp/tnguye14/bioinfo/Bi622/demultiplex/part1/read_qscores.py"
input_file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
output_file_prefix="1294_S1_L008_R3_001.fastq.gz"

/usr/bin/time -v python $script_file -f $input_file -o $output_file_prefix


exit