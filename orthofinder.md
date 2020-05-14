### Running OrthoFinder to Find Ortholog Transcripts With S. purpuratus

[OrthoFinder](https://github.com/davidemms/OrthoFinder)

Create directory:  
in: /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder  
`mkdir Orthofinder`

Copy in S. purpuratus peptide file:  
`scp /users/maggieschedl/Desktop/SPU_peptide.fasta meschedl@bluewaves:/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/Orthofinder/`

[OrthoFinder Manual](https://github.com/davidemms/OrthoFinder/blob/master/OrthoFinder-manual.pdf)

Need to put both the protein file from transdecoder and the SPU file into their own directory to call on when using OrthoFinder

Get the .pep file and rename to be specific for EL:  
`cd cdhit90El.Trinity.fasta.transdecoder_dir`  
`cp longest_orfs.pep ../LORF_EL90.pep.fasta`  
`cd..`  
Move file to Orthofinder directory and make Fasta directory inside:  
`mv LORF_EL90.pep.fasta Orthofinder/`  
`cd Orthofinder`  
`mkdir Fastas`  
`mv LORF_EL90.pep.fasta Fastas/`  
Bring in SPU file:  
`mv SPU_peptide.fasta Fastas/`

Options:  

"The '-t' option should always be used, typically with as many cores as are available. This determines how many highly-parallelisable tasks such as DIAMOND/BLAST searches, MSAs etc are run in parallel."

"RAM availability is an important consideration when using the `-a` option. Each thread loads all BLAST hits between one species and all sequences in all other species. To give some very approximate numbers, each thread might require:  
• 0.02 GB per species for small genomes (e.g. bacteria)  
• 0.04 GB per species for larger genomes (e.g. vertebrates)  
• 0.2 GB per species for even larger genomes (e.g. plants)  
I.e. running an analysis on 10 vertebrate species with 5 threads for the OrthoFinder algorithm (-a 5) might require 10 x 0.04 = 0.4 GB per thread and so 5 x 0.4 = 2 GB of RAM in total. If you have the BLAST results already then the total size of all the Blast* 0.txt files gives a good approximation of the memory requirements per thread. Additionally, the speed at which files can be read is likely to be the limiting factor when using more than 5-10 threads on current architectures so you may not see any increases in speed beyond this."

Say this is vert size, 10 threads, 0.04 * 2 * 10 is 0.8GB?

I'll just give it 4GB and 10 threads to be safe. Don't think I need to use the -a option.

Dependancies are:  
DIAMOND  
MCL  
FastME  
and potentially BLAST+

**Need to say othofinder._py_ on script for some reason specific to Bluewaves** (this is not in the manual or documentation on their site)

Make SLURM Job script:

`nano OrthoEL.sh`

```
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

```
`sbatch OrthoEL.sh`

Submitted job 1130084
