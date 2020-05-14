#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/

module load TransDecoder/5.5.0-foss-2019a

TransDecoder.LongOrfs -t cdhit90El.Trinity.fasta
