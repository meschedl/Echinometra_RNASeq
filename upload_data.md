### Uploading Data to Bluewaves

Created directory in /data/pradalab for these sequences

```
cd /data/pradalab/
mkdir meschedl
cd meschedl
mkdir Echinometra
```

Downloaded sequences to an external hard drive so they could be accessed by the command line then secure copied to the directory in Bluewaves.

```
scp -r /volumes/MES_Passport/Project_Sequences/HY_* meschedl@bluewaves:/data/pradalab/meschedl/Echinometra/
scp -r /volumes/MES_Passport/Project_Sequences/EV_* meschedl@bluewaves:/data/pradalab/meschedl/Echinometra/
scp -r /volumes/MES_Passport/Project_Sequences/EL_* meschedl@bluewaves:/data/pradalab/meschedl/Echinometra/
```


scp /volumes/MES_Passport/SPU_peptide.fasta meschedl@bluewaves:/data/pradalab/meschedl/Echinometra/
