#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=5G
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/

module load CD-HIT/4.8.1-foss-2018b

cd-hit-est -i EL.Trinity.fasta -o cdhit90El.Trinity.fasta -n 8 -c 0.90 -M 5000 -T 10 -d 0 -g 1
