#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=4GB
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/Orthofinder/

module load DIAMOND/0.9.25-foss-2018b
module load  MCL/14.137-foss-2018b
module load FastME/2.1.6.1-foss-2018b
module load BLAST+/2.8.1-foss-2018b
module load OrthoFinder/2.3.3-foss-2018b-Python-2.7.15

orthofinder.py -f ./Fastas/ -t 10
