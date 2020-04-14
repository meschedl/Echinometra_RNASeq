# Sequence QC


/data/pradalab/meschedl/Echinometra


First wanted to Get read counts in all of the files. There are 81 so doing them all in the working node will take too much computing power.

Make a batch job that takes all the

`nano read-count.sh`

```
 #!/bin/bash
 #SBATCH -t 2:00:00
 #SBATCH --nodes=1 --ntasks-per-node=1
 #SBATCH --export=NONE
 #SBATCH -D /data/pradalab/meschedl/Echinometra

for fq in *.fastq.gz
  do
  echo $fq >> read-counts.txt
  zcat $fq | echo $((`wc -l`/4)) >> read-counts.txt
  done
```

`sbatch read-count.sh`

Submitted batch job 1117192


Now I guess I want to know the data quality. I like to look at MultiQC reports, so first I will run fasqc on all the raw sequence files then compile them into a MultiQC report.

use module avail to see if it is loadable

I don't know how long this will take so at first I'll give it 3 hours?

nano fast-qc.sh

```
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra

module load FastQC/0.11.5-Java-1.8.0_92

fastqc *fastq.gz
```

Ok that wasn't long enough but it was close

nano fast-qc.sh

```
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra

module load FastQC/0.11.5-Java-1.8.0_92

fastqc HY_*
```
Submitted batch job 1117347


```
fastqc EV_EGGS_3_S66_L006_R1_001.fastq.gz
```

```
mv *.html data-QC/
mv *.zip data-QC/
cd data-qc
mkdir EL-before
mkdir EV-before
mkdir HY-before
mv EL_* EL-before/
mv EV_* EV-before/
mv HY_* HY-before/
```



```
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/data-QC/EL-before

module load MultiQC/1.7-foss-2018b-Python-2.7.15

multiqc *fastqc.zip
```
Submitted batch job 1118333

Had to use this version of the module MultiQC/1.7-foss-2018b-Python-2.7.15 because of error with matplotlib where the version's weren't compatible. Not sure exactly what the problem was but it ran! Requirement.parse('matplotlib<3.0.0,>=2.1.1'))

also was really quick so I probably don't need to make a job for them


`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EL-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

`scp -r meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EL-before/multiqc_data /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

Changed names on desktop to be EL_multiqc etc

```
cd ..
cd EV-before
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc *fastqc.zip
```
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EV-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

`scp -r meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EV-before/multiqc_data /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

Changed names on desktop to be EV_multiqc etc

```
cd ..
cd HY-before
multiqc *fastqc.zip
```

`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/HY-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

`scp -r meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/HY-before/multiqc_data /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

Changed names on desktop to be EV_multiqc etc
