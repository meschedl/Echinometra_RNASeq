#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=4GB
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/

module load BLAST+/2.9.0-iimpi-2019a


#Prepare the protein database for blast searches by:
makeblastdb -in uniprot_sprot.pep -dbtype prot



#Search Transdecoder-predicted proteins
blastp -query LORF_EL90.pep.fasta -db uniprot_sprot.pep -num_threads 10 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > ELblastp.outfmt6
