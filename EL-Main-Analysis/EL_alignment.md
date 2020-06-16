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

ran out of memory after 7 days!!

ok got to last sample
Says this:
Expression Results are written!
Time Used for EM.cpp : 0 h 04 m 00 s

rm -rf RSEM.temp

CMD: touch RSEM.isoforms.results.ok
slurmstepd: error: Exceeded step memory limit at some point.

think I might be ok if I just re-do the last sample eggs rep 3

nano eggs_align.sh

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

ok output looks exactly right for the singleton one

all the RSEM.genes.results files are the ones I need to
Build Transcript and Gene Expression Matrices

need to use the trinity perl scrip abundance_estimates_to_matrix.pl

EL_ortholog_filtered.fasta.gene_trans_map

/home/craker/diadema/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl \
--est_method salmon --gene_trans_map none --out_prefix salmon --name_sample_by_basedir \
pH_high_rep1/quant.sf pH_high_rep2/quant.sf pH_high_rep3/quant.sf pH_high_rep4/quant.sf \
pH_high_rep5/quant.sf pH_high_rep6/quant.sf pH_high_rep7/quant.sf pH_medx_rep1/quant.sf \
pH_medx_rep2/quant.sf pH_medx_rep3/quant.sf pH_medx_rep4/quant.sf pH_medx_rep5/quant.sf \
pH_medx_rep6/quant.sf pH_lowx_rep1/quant.sf pH_lowx_rep2/quant.sf pH_lowx_rep3/quant.sf \
pH_lowx_rep4/quant.sf pH_lowx_rep5/quant.sf pH_lowx_rep6/quant.sf pH_lowx_rep7/quant.sf \
pH_lowx_rep8/quant.sf

doing gene abundances
should also do transcript abundances?
do both I guess


`nano abundence-estimates.sh`

```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=10GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meschedl@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/

# load module
module load Trinity/2.8.4-foss-2016b
module load R-bundle-Bioconductor/3.3-foss-2016b-R-3.3.1

# calculate abundance matrix for genes

/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map EL_ortholog_filtered.fasta.gene_trans_map \
--name_sample_by_basedir \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_4cell_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_4cell_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_4cell_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_blast_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_blast_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_blast_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_gast_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_gast_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_gast_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_larv_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_larv_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/29_larv_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_4cell_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_4cell_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_4cell_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_blast_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_blast_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_blast_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_gast_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_gast_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_gast_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_larv_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_larv_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/33_larv_rep_3/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/eggs_rep_1/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/eggs_rep_2/RSEM.isoforms.results \
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/eggs_3/RSEM.isoforms.results \
```
sbatch abundence-estimates.sh

ok have to specify full path to each sample?

still have to use the name samples
had a problem with that one eggs sample...
there is a results file in eggs_rep_3 directory, the error did say it wrote the file then ran out of mem
change to use thatone

Submitted batch job 1312746



ok it's having problems because some transcripts don't have values for it to caluclate
Error, no TPM value specified for transcript [TRINITY_DN100003_c0_g1_i1] of gene [TRINITY_DN100003_c0_g1] for sample 29_4cell_rep_1 at /opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl line 322.

so I don't know what happened with that


oh apparently it only works for isofom results
not a problem because "When you include the --gene_trans_map file above, it will automatically generate the gene-level count and expression matrices"


ok now haveing another problem
* Outputting combined matrix.

/net/clusterhn.cluster.com/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/run_TMM_scale_matrix.pl --matrix RSEM.isoform.TPM.not_cross_norm > RSEM.isoform.TMM.EXPR.matrixCMD: R --no-save --no-restore --no-site-file --no-init-file -q < RSEM.isoform.TPM.not_cross_norm.runTMM.R 1>&2
sh: R: command not found
Error, cmd: R --no-save --no-restore --no-site-file --no-init-file -q < RSEM.isoform.TPM.not_cross_norm.runTMM.R 1>&2  died with ret (32512)  at /net/clusterhn.cluster.com/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/run_TMM_scale_matrix.pl line 105.
Error, CMD: /net/clusterhn.cluster.com/opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/run_TMM_scale_matrix.pl --matrix RSEM.isoform.TPM.not_cross_norm > RSEM.isoform.TMM.EXPR.matrix died with ret 9728 at /opt/software/Trinity/2.8.4-foss-2016b/trinityrnaseq-Trinity-v2.8.4/util/abundance_estimates_to_matrix.pl line 385.


could be a bad script but also looks like I need to use R?? https://github.com/trinityrnaseq/trinityrnaseq/issues/345

see if works with R-bundle-Bioconductor/3.3-foss-2016b-R-3.3.1 ??

ok might be working now

yay completed really quickly!!!

"The 'counts.matrix' file is used for downstream analyses of differential expression. The TMM.EXPR.matrix file is used as the gene expression matrix in most other analyses. "

want to use DESeq2 on RSEM.gene.counts.matrix

scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/alignment/RSEM.gene.counts.matrix /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/

create file with treatment and life stage information

problem is that the files are normalized, can I run again without any normalization??
added --cross_sample_norm none \
Submitted batch job 1373320

these are still not integers... not sure what is going to happen or what to do...

maybe I shoulnd't use the gene_trans_map
"When you include the --gene_trans_map file above, it will automatically generate the gene-level count and expression matrices, using the 'scaledTPM' method as described in txImport but implemented here directly in the Trinity script. This 'scaledTPM' method for estimating gene counts accounts for differences in isoform lengths that could otherwise lead to false gene DE reporting under situations where it is differential transcript usage (DTU) as opposed to differential gene expression (DGE) occurring."
--gene_trans_map EL_ortholog_filtered.fasta.gene_trans_map \
Submitted batch job 1373371

ok need to say it but say non? --gene_trans_map none \
Submitted batch job 1373374
these are still not integers...

ugh but I wan the gene ones though not isoform!!!


need to use tximport R program to do the transcript abundance estimation

need to have all the RSEM.gene.results files on my computer.
secure copy all of them
rename all of them

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
