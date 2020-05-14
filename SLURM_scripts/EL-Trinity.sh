#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meschedl@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/


module load Trinity/2.8.4-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load Jellyfish/2.2.6-foss-2016b
module load Salmon/0.10.2-foss-2016b-Python-2.7.12


Trinity --seqType fq --SS_lib_type F --min_contig_length 300 --max_memory 500G \
--single EL_29_4_C1_S49_L005_R1_001.fastq.gz.trim.fq.gz,\
EL_33_4_C1_S25_L003_R1_001.fastq.gz.trim.fq.gz,\
EL_29_B1_S2_L001_R1_001.fastq.gz.trim.fq.gz,\
EL_33_B2_S79_L007_R1_001.fastq.gz.trim.fq.gz,\
EL_29_G3_S65_L006_R1_001.fastq.gz.trim.fq.gz,\
EL_33_G1_S26_L003_R1_001.fastq.gz.trim.fq.gz,\
EL_29_L2_S11_L001_R1_001.fastq.gz.trim.fq.gz,\
EL_33_L3_S86_L008_R1_001.fastq.gz.trim.fq.gz,\
EL_EGGS_2_S9_L001_R1_001.fastq.gz.trim.fq.gz \
--CPU 20 --full_cleanup
