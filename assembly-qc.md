

Job 1125446 Finished.  Final Trinity assemblies are written to /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/trinity_out_dir.Trinity.fasta


less trinity_out_dir.Trinity.fasta
ls -l trinity_out_dir.Trinity.fasta
-rw-r--r-- 1 meschedl pradalab 328120693 May  7 07:16 trinity_out_dir.Trinity.fasta
ls -s trinity_out_dir.Trinity.
320431 trinity_out_dir.Trinity.fasta

copied to directory above
`cp trinity_out_dir.Trinity.fasta ../`

secure copied to my computer and then Google Drive
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/trinity_out_dir.Trinity.fasta /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq`


working in interactive node

Count the number of clusters from output of Trinity


cat trinity_out_dir.Trinity.fasta | grep '>' | wc -l
250831


How to get the best quality transcriptome
https://github.com/trinityrnaseq/trinity_community_codebase/wiki/Trinity-best-transcript-set

run transdecoder on output

mkdir Transdecoder
ln -s trinity_out_dir.Trinity.fasta Transdecoder/

module avail to look for transdecoder name
TransDecoder/5.5.0-foss-2019a let's see if this one will work with my transcriptome

module load TransDecoder/5.5.0-foss-2019a

TransDecoder.LongOrfs -t trinity_out_dir.Trinity.fasta

for some reason this didn't work, then it broke the link in that directory

removed link with rm trinity_out_dir.Trinity.fasta
if you use unlink it deletes the file!! why!!

cp trinity_out_dir.Trinity.fasta Transdecoder/EL.Trinity.fasta

now can I run
TransDecoder.LongOrfs -t EL.Trinity.fasta
yes it's running, might need to make a job script for this

exit out of interactive node for this

nano transdecoder.sh

```
#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/Transdecoder

module load TransDecoder/5.5.0-foss-2019a

TransDecoder.LongOrfs -t EL.Trinity.fasta
```

sbatch transdecoder.sh
Submitted batch job 1128931



looks like at some point they reduced some outputs to be only from cd-hit so how about I start with that one
