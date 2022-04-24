# LocalHGT

## Install
```
cd LocalHGT/
make
```
## Dependency
```
Python Module: scipy, numpy, pandas, sklearn, pysam, scikit-bio, Bio, pyfaidx

Please install samtools and bwa, and add them to system path.
samtools==1.11
bwa-0.7.17
```

## Reference database
```
Please contruct the database before running. The database (a single file in fasta formate) should contain all the representative references of your concerning envirment. You can 
For human gut, we can download it from the UHGG database.

Please use samtools to index the database first by

samtools faidx ref_database.fasta

LocalHGT also needs a index file, with the first running, Local will build the index the database automatically, and it will take several hours.
```
## Running
```
For each sample, with the unzipped fastq files, we can perform LocalHGT like
bash pipeline.sh reference_database.fasta sample.1.fq sample.2.fq sample_name result_dir/ 0.1 0.08
```

