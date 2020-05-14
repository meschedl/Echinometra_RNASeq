 cd /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD

 mkdir Trinnotate

Need to build protein databases
first load module


 module load Trinotate/3.2.1

Not sure if this code will work...

TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
also does not work

no somehow I have to get it to run this command but not sure how



/opt/software/Trinotate/3.2.1/Trinotate/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate

it is a perl script so might need to load perl??

module load Perl/5.30.0-GCCcore-8.3.0

ok I jsut downloaded wget https://raw.githubusercontent.com/Trinotate/Trinotate/master/admin/Build_Trinotate_Boilerplate_SQLite_db.pl

make exicuatble?

chmod +x Build_Trinotate_Boilerplate_SQLite_db.pl
check that it's exicuable
ls -l Build_Trinotate_Boilerplate_SQLite_db.pl
-rwxr-xr-x 1 meschedl pradalab 5766 May 13 09:36 Build_Trinotate_Boilerplate_SQLite_db.pl

perl Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

hmmm this did not work

Can't locate Pipeliner.pm in @INC (you may need to install the Pipeliner module) (@INC contains: /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/../PerlLib /opt/slurm/lib64/perl5/ /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0 /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0) at Build_Trinotate_Boilerplate_SQLite_db.pl line 8.
BEGIN failed--compilation aborted at Build_Trinotate_Boilerplate_SQLite_db.pl line 8.

ok what's pipliner

ok it's in the perl library in the trinnotate
and in the script it looks for PerlLib
which is only with the full download I guess


uhhh what if I try another way

/opt/modulefiles/bio/Trinotate/Trinotate-Trinotate-v3.2.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

/opt/modulefiles/bio/Trinotate/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate


 /opt/modulefiles/bio/Trinotate


from Kevin Bryan

/opt/software/Trinotate/3.2.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

This one worked



created Pfam-A.hmm.gz  Trinotate.sqlite  uniprot_sprot.dat.gz  uniprot_sprot.pep

Now to capture BLAST Homologies

Latest BLAST
BLAST+/2.9.0-iimpi-2019a

hopefully this is right.

Want to run this on the transdecoder open reading frame file. But it can also be run on the trinity output

Need .pep file
/data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Transdecoder/cdhit90El.Trinity.fasta.transdecoder_dir

hmm I should have put this directory in the Transdecoder directory. Oh well.

use same name as before

cp longest_orfs.pep ../../LORF_EL90.pep.fasta
cd ../..
mv LORF_EL90.pep.fasta Trinnotate/


needs the uniprot_sprot.pep that I just made and the transdecoder.pep file


nano BLAST_pep.sh
 ```
#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=4GB
#SBATCH --account=pradalab
#SBATCH --export=NONE
#SBATCH -D /data/pradalab/meschedl/Echinometra/trimmed-data/EL_trimmed/EL-CD/Trinnotate/

module load BLAST+/2.9.0-iimpi-2019a


#Prepare the protein database for blast searches by:  
makeblastdb -in uniprot_sprot.pep -dbtype prot



#Search Transdecoder-predicted proteins  
blastp -query LORF_EL90.pep.fasta -db uniprot_sprot.pep -num_threads 10 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > ELblastp.outfmt6
```

Submitted batch job 1131108
