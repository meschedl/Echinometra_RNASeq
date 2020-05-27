### Trinotate

Need to run sequence analyses for the transcriptome before Trinotate will compile them into a database.

Steps:
1. Build protein databases
2. Run BLASTp
3. Run BLASTx
4. Run HMMR
5. Compile into annotation database

Create directory to work in:  
 `cd /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD`  
`mkdir Trinnotate`

1. Build protein databases

Had problems running Trinotate perl script, had to specify exact path to script in Bluewaves. It will not run correctly if you just download the script from GitHub. All Bluewaves software lives in these directories /opt/software/

In interactive mode:  

`/opt/software/Trinotate/3.2.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate`  
Created Pfam-A.hmm.gz  Trinotate.sqlite  uniprot_sprot.dat.gz  uniprot_sprot.pep

2. Capture BLAST homologies for protein

This is run on both the open-reading frame peptide file from Transdecoder.

Navagate to peptide file:  
`cd /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/cdhit90El.Trinity.fasta.transdecoder_dir`


Transfer it to Trinotate directory:  
`cp longest_orfs.pep ../../LORF_EL90.pep.fasta`  
`cd ../..`  
`mv LORF_EL90.pep.fasta Trinnotate/`  


Make job script that prepares protein database and runs BLASTp:  

Options for E value:  
The BLAST E-value is the number of expected hits of similar quality (score) that could be found just by chance.
E-value of 10 means that up to 10 hits can be expected to be found just by chance, given the same size of a random database. Default for Trinotate suggests 1e-3 (below .01 should be good).  
The also recommend using --max_target_seqs 1 for all BLAST, however BLAST recommends 5. I did 1, but this may need to change.


`nano BLAST_pep.sh`
 ```
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
```
`sbatch BLAST_pep.sh`
Submitted batch job 1131108

3. Capture BLAST homologies for DNA

Navigate to CD-HIT filtered assembly and copy to Trinotate directory:  
`cd /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD`  
`cp cdhit90El.Trinity.fasta Trinnotate/`
`cd Trinnotate/`

Make job script that runs BLASTx:   
`nano BLAST_X.sh`
```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=4GB
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/

module load BLAST+/2.9.0-iimpi-2019a

#Search Trinity transcripts
blastx -query cdhit90El.Trinity.fasta -db uniprot_sprot.pep -num_threads 10 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > ELblastx.outfmt6
```
`sbatch BLAST_X.sh`  
Submitted batch job 1132885

4. Running HMMER to identify protein domains

this uses the transdecoder pep file

nano HMMR.sh
 ```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=4GB
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/

module load HMMER/3.2.1-foss-2018b

# Uncompress and prepare the Pfam database for use with 'hmmscan' like so:

gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

# Running HMMER to identify protein domains
hmmscan --cpu 10 --domtblout TrinotatePFAM.out Pfam-A.hmm LORF_EL90.pep.fasta > ELpfam.log
```
sbatch HMMR.sh
Submitted batch job 1133164


Now  I have:   
 ELpfam.log an TrinotatePFAM.out from HMMR  
 ELblastx.outfmt6 from Blastx  
 and ELblastp.outfmt6 from Blastp  

 **Now to load them into the annotation report**  
 https://github.com/Trinotate/Trinotate.github.io/wiki/Loading-generated-results-into-a-Trinotate-SQLite-Database-and-Looking-the-Output-Annotation-Report

 1. Load transcripts and coding regions

 I need "Gene/Transcript relationships (tab delimited format: "gene_id(tab)transcript_id", same as used by the RSEM software)." This uses the same assembly fasta file as BLASTx.  

This is pretty quick and does not need to be a job, it can run in `interactive` mode. I need to specify the full path to the perl script it uses to make the relationship file.


`module load Trinity/2.8.4-foss-2016b`

`/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl cdhit90El.Trinity.fasta >  ELCD.fasta.gene_trans_map`


Now put that map into the sqlite database:

Still in `interactive` mode

`module load Trinotate/3.2.1`

`Trinotate Trinotate.sqlite init --gene_trans_map ELCD.fasta.gene_trans_map --transcript_fasta cdhit90El.Trinity.fasta --transdecoder_pep LORF_EL90.pep.fasta`


2. Loading BLAST homologies

**Load protein hits**:  
`Trinotate Trinotate.sqlite LOAD_swissprot_blastp ELblastp.outfmt6`


**Load transcript hits**:   
`Trinotate Trinotate.sqlite LOAD_swissprot_blastx ELblastx.outfmt6`     


3. Load Pfam domain entries (from HMMR)

`Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out`



4. Generate output of full Trinotate Anotation Report


`nano trinotate_report.sh`
```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/


module load Trinotate/3.2.1


Trinotate Trinotate.sqlite report -E 1e-3 > EL_trinotate_annotation_report.xls
```
`sbatch trinotate_report.sh`  
Submitted batch job 1133759

Main file output is: EL_trinotate_annotation_report.xls



Secure copy to computer to put on github:  
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/EL_trinotate_annotation_report.xls /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/Out_Files/`
`
