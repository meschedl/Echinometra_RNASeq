#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=20G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trimmomatic/0.38-Java-1.8

for fastq in *.fastq.gz
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE $fastq $fastq.trim.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 HEADCROP:11 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:40

done
