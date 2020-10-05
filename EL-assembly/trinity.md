## Running Trinity on EL Trimmed Sequence Files

Sequences to use:

EL_29_4_C1_S49_L005_R1_001.fastq.gz.trim.fq.gz  
EL_33_4C1_S25_L003_R1_001.fastq.gz.trim.fq.gz    
EL_29_B1_S2_L001_R1_001.fastq.gz.trim.fq.gz    
EL_33_B2_S79_L007_R1_001.fastq.gz.trim.fq.gz    
EL_29_G3_S65_L006_R1_001.fastq.gz.trim.fq.gz    
EL_33_G1_S26_L003_R1_001.fastq.gz.trim.fq.gz    
EL_29_L2_S11_L001_R1_001.fastq.gz.trim.fq.gz    
EL_33_L3_S86_L008_R1_001.fastq.gz.trim.fq.gz    
EL_EGGS_2_S9_L001_R1_001.fastq.gz.trim.fq.gz    

This is one from each temperature and each developmental stage. It's very important to include all of these.

### [Running Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)

Required options:  

--seqType <string>      : type of reads: ( fa, or fq )

--max_memory <string>      : suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc) provided in Gb of RAM, ie.  '--max_memory 10G'

If paired reads:  
  --left  <string>    : left reads, one or more file names (separated by commas, not spaces)  
  --right <string>    : right reads, one or more file names (separated by commas, not spaces)  
Or, if unpaired reads:  
  --single <string>   : single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )


I have fastq (fq) files that are single end reads.  
The maximum memory I can request on BLuewaves is 500GB.

Other important options:  
--SS_lib_type <string>     : Strand-specific RNA-Seq read orientation. If paired: RF or FR, if single: F or R.   (dUTP method = RF)  See web documentation.  
--CPU <int>                     : number of CPUs to use, default: 2  
--min_contig_length <int>       : minimum assembled contig length to report (def=200)  
--full_cleanup                  : only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta

Library time is forward read and minimum contig length should be 300bp. Going to request a full node (--exclusive flag) which has 20 CPUs.


Make SLURM Job scrip:

`nano EL-Tinity.sh`
```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meschedl@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/


module load Trinity/2.8.4-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load Jellyfish/2.2.6-foss-2016b
module load Salmon/0.10.2-foss-2016b-Python-2.7.12


Trinity --seqType fq --SS_lib_type F --min_contig_length 300 --max_memory 500G \
--single EL_29_4_C1_S49_L005_R1_001.fastq.gz.trim.fq.gz,\
EL_33_4_C1_S25_L003_R1_001.fastq.gz.trim.fq.gz,\
EL_29_B1_S2_L001_R1_001.fastq.gz.trim.fq.gz,\
EL_33_B2_S79_L007_R1_001.fastq.gz.trim.fq.gz,\
EL_29_G3_S65_L006_R1_001.fastq.gz.trim.fq.gz,\
EL_33_G1_S26_L003_R1_001.fastq.gz.trim.fq.gz,\
EL_29_L2_S11_L001_R1_001.fastq.gz.trim.fq.gz,\
EL_33_L3_S86_L008_R1_001.fastq.gz.trim.fq.gz,\
EL_EGGS_2_S9_L001_R1_001.fastq.gz.trim.fq.gz \
--CPU 20 --full_cleanup
```
`sbatch EL-Trinity.sh`  
Submitted job number 1125446
