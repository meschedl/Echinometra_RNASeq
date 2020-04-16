# Sequence QC


/data/pradalab/meschedl/Echinometra


First wanted to Get read counts in all of the files. There are 81 so doing them all in the working node will take too much computing power.

Make a batch job that takes all the files and counts the number of lines (reads) and then puts that information into a file.

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
Submit the job.  
`sbatch read-count.sh`

Submitted batch job 1117192


Now I guess I want to know the data quality. I like to look at MultiQC reports, so first I will run fasqc on all the raw sequence files then compile them into a MultiQC report.

Use `module avail` to see if it is loadable and what exactly to use.

I don't know how long this will take so at first I'll give it 3 hours?

`nano fast-qc.sh`

```
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra

module load FastQC/0.11.5-Java-1.8.0_92

fastqc *fastq.gz
```

Ok that wasn't long enough but it was close. Need to run all the HY files.

`nano fast-qc.sh`

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

And one last EV file.
```
fastqc EV_EGGS_3_S66_L006_R1_001.fastq.gz
```

Now move those into more organized directories.

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

Now I want to run Multiqc to make a compiled report. Do those by species because there should be some differences in those that will blind other sequence signals.

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

Had to use this version of the module MultiQC/1.7-foss-2018b-Python-2.7.15 because of error with matplotlib where the version's weren't compatible. Not sure exactly what was wrong but with tat version it ran. Error was: Requirement.parse('matplotlib<3.0.0,>=2.1.1'))

Also was really quick so I probably don't need to make a job for them

Secure copied reports to my computer to put on the repository.

`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EL-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

`scp -r meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EL-before/multiqc_data /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

Changed names on desktop to be EL_multiqc to not get confused with other files.

Ran muiltiqc on EV files.
```
cd ..
cd EV-before
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc *fastqc.zip
```
Transferred files to computer.   
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EV-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

`scp -r meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EV-before/multiqc_data /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

Changed names on desktop to be EV_multiqc etc


Ran multiqc on HY files.
```
cd ..
cd HY-before
multiqc *fastqc.zip
```
Transferred files to computer.
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/HY-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

`scp -r meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/HY-before/multiqc_data /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq/`

Changed names on desktop to be EV_multiqc etc


Now to quality trim the files. They look pretty good from the multqc reports, lots of high quality scores. Just looks like there is some adapter content and also all the bases at the beginning of the sequences are the same in all samples/reads. Trinity has the program Trimmomatic in it so I'll star with trying that.

Try on one sequence to see if it runs.
```
#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=20G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trinity/2.8.4-foss-2016b

for fastq in *.fastq.gz
do
Trinity --seqType fq --single $fastq --max_memory 20G --trimmomatic --quality_trimming_params "SLIDINGWINDOW:4:5 HEADCROP:5 LEADING:5 TRAILING:5 MINLEN:40"
gzip *.fq
done
```
Ok this didn't run, needs SAMtools installed. Maybe if I don't try running it a batch I'll be able to figure out the errors faster.

Try on one sequence

`module load Trinity/2.8.4-foss-2016b`  
`module load SAMtools/1.9-foss-2018b`

`Trinity --seqType fq --single EL_29_4_C1_S49_L005_R1_001.fastq.gz --max_memory 10G --trimmomatic --quality_trimming_params "SLIDINGWINDOW:4:5 HEADCROP:5 LEADING:5 TRAILING:5 MINLEN:40"`

Error, problem with the SAMtools I loaded I think. It needs to be a 2016b foss to match with Trinity version. `module load SAMtools/1.3.1-foss-2016b`

Another error, now need to have jellyfish loaded `module load Jellyfish/2.2.6-foss-2016b`

Another error, now I need salmon loaded `module load Salmon/0.10.2-foss-2016b-Python-2.7.12`

Now it is running. I think it might not output as a zipped file. So I will want to zip it when it's run on all the samples.


It took a long time to run and then it got stuck because it keeps going with Trinity after doing the trimming, which is not what I want because only a subset of those samples will be used for the transcriptome. And I want to be able to check the trimmed sequence quality.

Using Trimmomatic

It takes about 16 minutes for one file so I have to make a job script for testing.

`nano trim.sh`

```
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=10G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trimmomatic/0.38-Java-1.8

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE EL_29_4_C1_S49_L005_R1_001.fastq.gz EL_29_4_C1_S49_L005_R1_001_trimm_fq.gz SLIDINGWINDOW:4:5 HEADCROP:5 LEADING:5 TRAILING:5 MINLEN:40
```
That is with the recommended settings from Trinity, but minimum length set to 40bp because that's better, and cropping the first 5 bases to remove the overrepresented nucleotides at the beginning.

That ran fine so I want to check the quality with fastqc. It dropped no reads in trimmning

`module load FastQC/0.11.5-Java-1.8.0_92`

`fastqc EL_29_4_C1_S49_L005_R1_001_trimm_fq.gz`

Then copy to my computer so I can open it up in Chome.   
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/EL_29_4_C1_S49_L005_R1_001_trimm_fq_fastqc.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq`


Ok so I should have actually cut more at the beginning because it's to the first 11 bases that are all the same nucleotides actually, which you can tell in the per base sequence count.
```
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=10G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trimmomatic/0.38-Java-1.8

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE EL_29_4_C1_S49_L005_R1_001.fastq.gz EL_29_4_C1_S49_L005_R1_001_trimm_fq.gz SLIDINGWINDOW:4:5 HEADCROP:11 LEADING:5 TRAILING:5 MINLEN:40
```
`fastqc EL_29_4_C1_S49_L005_R1_001_trimm_fq.gz`

Copy to computer  
`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/data-QC/EL-before/EL_33_B2_S79_L007_R1_001_fastqc.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq`

Ok the beginning of the files are good but there is still adapter content. The fastqc tells me what type of adapters are in the files: illumina truseq adapters with specific indexes. So those need to be trimmed out. Used [Trimmomatic helpful manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) to figure out what file I needed for the adapters.

Ok I think the adapter file I need is TruSeq3-SE.fa

`wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa`

Ok try again with the adapter trimming, but don't do the head crop because that might be removed with the adapter trimming?

```
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=10G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trimmomatic/0.38-Java-1.8

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE EL_29_4_C1_S49_L005_R1_001.fastq.gz EL_29_4_C1_S49_L005_R1_001_trimmm_fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:5  LEADING:5 TRAILING:5 MINLEN:40
```

`fastqc EL_29_4_C1_S49_L005_R1_001_trimmm_fq.gz`

`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/EL_29_4_C1_S49_L005_R1_001_trimmm_fq_fastqc.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq`

Ok I still need to trim the first 11 bases, but now no longer says there is adapter in there. This makes sense because I think in the original multiqc reports it detected adapter sequences at the end not the beginning of reads.

All the files have the overrepresented nucleotides 11th bp so I think it would be fine to trim them all the same way. Just test this one file one more time to be sure.

```
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=10G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trimmomatic/0.38-Java-1.8

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE EL_29_4_C1_S49_L005_R1_001.fastq.gz EL_29_4_C1_S49_L005_R1_001_trim4_fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 HEADCROP:11 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:40
```
Submitted batch job 1119592


`fastqc EL_29_4_C1_S49_L005_R1_001_trim4_fq.gz`

`scp meschedl@bluewaves.uri.edu:/data/pradalab/meschedl/Echinometra/EL_29_4_C1_S49_L005_R1_001_trim4_fq_fastqc.html /Users/maggieschedl/Desktop/URI/Prada/Echinometra_RNASeq`

This looks good, now to run it on all  81 files.


First I'll move those testing files

```
mkdir test-trim
mv *.html test-trim/
mv *fastqc.zip test-trim/
mv *fq.gz test-trim/
```

To run all the files have to make a for loop and then have them also uniquely named.

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=20G
#SBATCH -D /data/pradalab/meschedl/Echinometra/

module load Trimmomatic/0.38-Java-1.8

for fastq in *.fastq.gz
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE $fastq $fastq.trim.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:5  LEADING:5 TRAILING:5 MINLEN:40

done
```
Submitted batch job 1119660


Now just to check I want to get a multiqc report for the trimmed files. Only the EL files are done being trimmed so I start with those. They were moved to the trimmed-data/ directory, then I made them a sub directory.
```
mkdir EL-trimmed
mv *.fq.gz EL-trimmed/
```

Make a job to run fastQC on them.

```
#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL-trimmed/

module load FastQC/0.11.5-Java-1.8.0_92

fastqc EL_*
```
