### Copying Trinity Fasta File and Running CD-Hit to Collapse Transcript Clusters

Final Trinity assemblies were written to /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/trinity_out_dir.Trinity.fasta


First check dimensions of file:  

`ls -l trinity_out_dir.Trinity.fasta`  
-rw-r--r-- 1 meschedl pradalab 328120693 May  7 07:16 trinity_out_dir.Trinity.fasta  
`ls -s trinity_out_dir.Trinity.`  
320431 trinity_out_dir.Trinity.fasta

Want to make sure this file does not get overwritten so I want to save it multiple places.  
Copy to above directory:  
`cp trinity_out_dir.Trinity.fasta ../`

Secure copied to my computer and then added to my Google Drive:  
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/trimmed-data/trinity_out_dir.Trinity.fasta /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq`

Changed to interactive mode in Bluewaves

Count the number of clusters from Trinity:  
`cat trinity_out_dir.Trinity.fasta | grep '>' | wc -l`  
**250831**

Resource on how to get best quality transcriptome from trinity [here](https://github.com/trinityrnaseq/trinity_community_codebase/wiki/Trinity-best-transcript-set), I looked at this a little.

### CD-Hit

First create separate directory and rename Trinity file to be nicer:  
`mkdir EL-CD`  
`cp trinity_out_dir.Trinity.fasta EL-CD/EL.Trinity.fasta`  
`cd EL-CD`

[CD-HIT-EST](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST)  
"CD-HIT-EST clusters a nucleotide dataset into clusters that meet a user-defined similarity threshold, usually a sequence identity. The input is a DNA/RNA dataset in fasta format and the output are two files: a fasta file of representative sequences and a text file of list of clusters. Since eukaryotic genes usually have long introns, which cause long gaps, it is difficult to make full-length alignments for these genes. So, CD-HIT-EST is good for non-intron containing sequences like EST."

Some options out of the long list to use:  

-i input file  
-o output file  
-n word size "-n 8,9    for thresholds 0.90 ~ 0.95"  
-c sequence identity threshold this is the default cd-hit's "global sequence identity" calculated as: number of identical bases in alignment divided by the full length of the shorter sequence  
-T number of threads  
-M	memory limit (in MB) for the program, default 800; 0 for unlimited  
-d	length of description in .clstr file, default 20,
 	if set to 0, it takes the fasta defline and stops at first space  
-g	1 or 0, default 0  
  	by cd-hit's default algorithm, a sequence is clustered to the first
  	cluster that meet the threshold (fast cluster). If set to 1, the program
  	will cluster it into the most similar cluster that meet the threshold
  	(accurate but slow mode)
  	but either 1 or 0 won't change the representatives of final clusters

Not sure what to set the number of threads to, first I checked some capabilities of Bluewaves:  
`lscpu`  
CPU(s):                20

Thread(s) per core:    1  
Core(s) per socket:    10  
Socket(s):             2  


CPUs = threads per core times cores per socket times sockets ( 1 X 10 X 2 = 20)

Not wanting to take up a whole node (20), I decided on 10. Guessed that 5GB is enough for memory (much more than 800MB)

Used `module avail` to check the module version I want: **CD-HIT/4.8.1-foss-2018b**  

Created SLURM Job Script:  

`nano cd-hit-90.sh`
```
#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=5G
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/

module load CD-HIT/4.8.1-foss-2018b

cd-hit-est -i EL.Trinity.fasta -o cdhit90El.Trinity.fasta -n 8 -c 0.90 -M 5000 -T 10 -d 0 -g 1
```

Submitted batch job 1128992

Output file:  
cdhit90El.Trinity.fasta  

Count number of clusters:  
`cat cdhit90El.Trinity.fasta | grep '>' | wc -l`  
**107376**

This is a dramatic decrease in clusters! 
