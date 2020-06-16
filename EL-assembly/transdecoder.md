### Running Transdecoder on CD-Hit Limited Assembly from Trinity

[Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki): Finds Coding Regions Within Transcripts

Create directory and move in CD-hit file:

`mkdir Transdecoder`  
`cp cdhit90El.Trinity.fasta Transdecoder/`  
`cd Transdecoder`


Make SLURM Job script (no other potions to optimize):

`nano EL-transdecoder.sh`

```
#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/

module load TransDecoder/5.5.0-foss-2019a

TransDecoder.LongOrfs -t cdhit90El.Trinity.fasta
```
`sbatch EL-transdecoder.sh`  


Submitted batch job 1129485

Directions in outfile: "Use file: /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/cdhit90El.Trinity.fasta.transdecoder_dir/longest_orfs.pep  for Pfam and/or BlastP searches to enable homology-based coding region identification."


Count number of clusters now:  
`cd cdhit90El.Trinity.fasta.transdecoder_dir`

`cat longest_orfs.cds | grep '>' | wc -l`  
**50583**

This is the coding sequences file. The longest_orfs.pep file has the protein sequences. 
