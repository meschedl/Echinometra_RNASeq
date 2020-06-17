### Mapping and Transcript Abundance Estimations

Using the [Trinity pipeline](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#alignment-free-abundance-estimation) to use Bowtie2 for read mapping and RSEM for transcript abundance estimation. This is all done in one step.

Working directory: /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment  
Make directory, move in ortholog filtered assembly file and link trimmed sample sequences:   
`mkdir alignment`  
`mv EL_ortholog_filtered.fasta alignment/`  
`ln -s /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/*.trim.fq.gz .`  

Need to make a tab-delimited file of sample names: "Please use the --samples_file parameter with the abundance estimation utility. This will organize your outputs so that each replicate will be organized in its own output directory named according to the corresponding replicate"

`nano EL-samples.txt`

29_4cell  29_4cell_rep_1  EL_29_4_C1_S49_L005_R1_001.fastq.gz.trim.fq.gz  
29_4cell  29_4cell_rep_2  EL_29_4_C2_S57_L005_R1_001.fastq.gz.trim.fq.gz  
29_4cell  29_4cell_rep_3  EL_29_4_C3_S64_L006_R1_001.fastq.gz.trim.fq.gz  
29_blast  29_blast_rep_1  EL_29_B1_S2_L001_R1_001.fastq.gz.trim.fq.gz  
29_blast  29_blast_rep_2  EL_29_B2_S10_L001_R1_001.fastq.gz.trim.fq.gz  
29_blast  29_blast_rep_3  EL_29_B3_S18_L002_R1_001.fastq.gz.trim.fq.gz  
29_gast 29_gast_rep_1 EL_29_G1_S50_L005_R1_001.fastq.gz.trim.fq.gz  
29_gast 29_gast_rep_2 EL_29_G2_S58_L005_R1_001.fastq.gz.trim.fq.gz  
29_gast 29_gast_rep_3 EL_29_G3_S65_L006_R1_001.fastq.gz.trim.fq.gz  
29_larv 29_larv_rep_1 EL_29_L1_S3_L001_R1_001.fastq.gz.trim.fq.gz  
29_larv 29_larv_rep_2 EL_29_L2_S11_L001_R1_001.fastq.gz.trim.fq.gz  
29_larv 29_larv_rep_3 EL_29_L3_S19_L002_R1_001.fastq.gz.trim.fq.gz  
33_4cell  33_4cell_rep_1  EL_33_4_C1_S25_L003_R1_001.fastq.gz.trim.fq.gz  
33_4cell  33_4cell_rep_2  EL_33_4_C2_S33_L003_R1_001.fastq.gz.trim.fq.gz  
33_4cell  33_4cell_rep_3  EL_33_4_C3_S41_L004_R1_001.fastq.gz.trim.fq.gz  
33_blast  33_blast_rep_1  EL_33_B1_S72_L007_R1_001.fastq.gz.trim.fq.gz  
33_blast  33_blast_rep_2  EL_33_B2_S79_L007_R1_001.fastq.gz.trim.fq.gz  
33_blast  33_blast_rep_3  EL_33_B3_S85_L008_R1_001.fastq.gz.trim.fq.gz  
33_gast 33_gast_rep_1 EL_33_G1_S26_L003_R1_001.fastq.gz.trim.fq.gz  
33_gast 33_gast_rep_2 EL_33_G2_S34_L003_R1_001.fastq.gz.trim.fq.gz  
33_gast 33_gast_rep_3 EL_33_G3_S42_L004_R1_001.fastq.gz.trim.fq.gz  
33_larv 33_larv_rep_1 EL_33_L1_S73_L007_R1_001.fastq.gz.trim.fq.gz  
33_larv 33_larv_rep_2 EL_33_L2_S80_L007_R1_001.fastq.gz.trim.fq.gz  
33_larv 33_larv_rep_3 EL_33_L3_S86_L008_R1_001.fastq.gz.trim.fq.gz  
eggs  eggs_rep_1  EL_EGGS_1_S1_L001_R1_001.fastq.gz.trim.fq.gz  
eggs  eggs_rep_2  EL_EGGS_2_S9_L001_R1_001.fastq.gz.trim.fq.gz  
eggs  eggs_rep_3  EL_EGGS_3_S17_L002_R1_001.fastq.gz.trim.fq.gz  



Need to load SAMtools, Trinity, and RSEM.

Need to specify full path to perl script /opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl

I had a hard time optimizing the amount of memory needed, the job would run for hours then fail. Decided on 20GB of memory.  
This in theory should be 2GB for each thread? Decided to use 10 threads because it's not asking for too much computing power but more than the default (4).

`nano read_align_transcript_quant.sh`
```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=20GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meschedl@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/

# load modules
module load Trinity/2.8.4-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load RSEM/1.3.0-foss-2016b

#Aligning reads with Bowtie2 and quantify transcripts with RSEM

/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --transcripts EL_ortholog_filtered.fasta --seqType fq --samples_file EL-samples.txt --est_method RSEM --aln_method bowtie2 --thread_count 10 --trinity_mode --prep_reference --output_dir ELrsem_outdir

```
`sbatch read_align_transcript_quant.sh`  
Submitted batch job 1216093

Ran out of memory after 7 days, only had one sample left to go. Did that one separately.

`nano eggs_align.sh`

```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=23GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meschedl@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/

# load modules
module load Trinity/2.8.4-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load RSEM/1.3.0-foss-2016b

/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --transcripts EL_ortholog_filtered.fasta --seqType fq --single EL_EGGS_3_S17_L002_R1_001.fastq.gz.trim.fq.gz --est_method RSEM --aln_method bowtie2 --thread_count 10 --trinity_mode --output_dir eggs_3

```
Submitted batch job 1307802

This made separate directories for each sample, with an RSEM.genes.results and RSEM.isoforms.results files in them with count values like FPKM and CPM.


Now to get the counts file, I do NOT want to use the Trinity abundance estimation script, because it will make non-integer count estimates which are incompatible with the DESeq2 differential expression R package.


I need to use tximport R program to do the transcript abundance estimation

I need to have all the RSEM.gene.results files on my computer.
secure copy all of them
rename all of them to "29_4cell_rep_1.RSEM.genes.results" style to be able to load them into R

```
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_4cell_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_4cell_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_4cell_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_blast_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_blast_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_blast_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_gast_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_gast_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_gast_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_larv_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_larv_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_larv_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_4cell_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_4cell_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_4cell_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_blast_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_blast_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_blast_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_gast_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_gast_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_gast_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_larv_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_larv_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_larv_rep_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/eggs_rep_1/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/eggs_rep_2/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/eggs_3/RSEM.genes.results /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/EL-DESeq2/
```
