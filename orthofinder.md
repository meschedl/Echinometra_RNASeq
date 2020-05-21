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

Results are in /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/Orthofinder/Fastas/OrthoFinder/Results_May11_1

Orthologs are in Orthologs directory
data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/Orthofinder/Fastas/OrthoFinder/Results_May11_1/Orthologues/Orthologues_LORF_EL90.pep

Comparisson file between EL and SPU is called LORF_EL90.pep__v__SPU_peptide.tsv.  
Count number of lines (should be number of orthologs):  
`wc -l LORF_EL90.pep__v__SPU_peptide.tsv`  
**13012 LORF_EL90.pep__v__SPU_peptide.tsv**

There is also a SPU_peptide__v__LORF_EL90.pep.tsv file that has the same number of lines


I want to take the trinity transcripts in the 2nd column of this file and only keep those in the cut down transcriptome file (cdhit90El.Trinity.fasta)

Cut the second column into a new file.  
`cut -f2 LORF_EL90.pep__v__SPU_peptide.tsv > EL_otholog_transcripts.txt`

There are multiple transcripts in a line separated by commas, hopefully this will work for filtering the assembly file.

Copy to directory with cd-hit fasta.  
`cp EL_otholog_transcripts.txt /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/`


Using [BBMap script](https://github.com/BioInfoTools/BBMap/blob/master/sh/filterbyname.sh) filterbyname.sh to filter the assembly file to only ortholog transcripts.

Options:  
in2 and out2 are for paired reads and are optional  
include=f           Set to 'true' to include the filtered names rather than excluding them.  
substring=f         Allow one name to be a substring of the other, rather than a full match.
                         f: No substring matching.  
                         t: Bidirectional substring matching.  
                         header: Allow input read headers to be substrings of names in list.  
                         name: Allow names in list to be substrings of input read headers.    
For improved speed, add 'usejni=t' to the command line of BBMap tools which support the use of the compiled jni C code.

`nano ortho_filter.sh`

```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/

module load BBMap/38.81-foss-2018b-Java-1.8

filterbyname.sh in=cdhit90El.Trinity.fasta names=EL_otholog_transcripts.txt out=EL_ortholog_filtered.fasta include=t substring=t usejni=t
```
Submitted batch job 1133459

Once done the ortholog filtered assembly file is EL_ortholog_filtered.fasta.  
Count the number of transcripts:  
`cat EL_ortholog_filtered.fasta | grep '>' | wc -l`  
**15474**
