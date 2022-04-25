# LocalHGT: an ultrafast HGT detection method from large microbial communities

## Install
```
cd LocalHGT/
make
```
## Dependency
```
Python 3.8+
Python Module: scipy, numpy, pandas, sklearn, pysam, scikit-bio, Bio, pyfaidx
```
Please install samtools and bwa, and add them to system path, the version 
should be
```
samtools==1.11+
bwa-0.7.17+
```

## Reference database
Please construct the database before running. The database (a single file in fasta 
formate) should contain all the representative references of your concerning environment. 

For the human gut, we can download it from the Unified Human Gastrointestinal Genome 
([UHGG](https://www.nature.com/articles/s41587-020-0603-3)) database. 
The script build_UHGG_reference.py can be used to download the UHGG v1 genomes.

Please use samtools to index the database first by

```
samtools faidx ref_database.fasta

```
LocalHGT also needs an index file, with the first running, Local will index 
the database automatically, and it will take several hours.
## Running
For each paired-end sequencing sample, with the unzipped fastq files, we can perform
LocalHGT like

```
bash pipeline.sh reference_database.fasta sample.1.fq sample.2.fq sample_name result_dir/ 0.1 0.08
```
## Output interpretion
The HGT breakpoints would be saved in the *acc.csv file.
Each line contains a pair of breakpoints, it records the two genome's ID and 
the breakpoint positions. Here is an example:
```
from_ref,to_ref,from_pos,to_pos,from_side,to_side,if_reverse,read_seq,ref_seq,similarity
GUT_GENOME000729_19,4074,GUT_GENOME000139_8,30127,right,right,True,CAATCCTTCTATTACGTCAGGGGAAATGGTAGCCGCCCCGATTCCATACTCACAGAGGTCG,CAATCCTTCTATTACGTCAGGGGAAATGGTAGCCGCCCCGATTCCATACTCACAGAGTTCA,1.885
GUT_GENOME000729_19,39463,GUT_GENOME000139_8,30121,left,left,True,AAGCATCCCTTCTCCCATGCCGAATTCTTTCGCGGCCGAAAA,AGCATCCCTTCTCCCATGCCGAATTCTTTCGCGGCCGAAAAC,1.952
GUT_GENOME000232_7,104924,GUT_GENOME044018_8,2495,right,left,False,CATTCTCTCTGGCAGTTACACTGATTTTACCTCCTGCAGGCACATG,CATTCTCTCTGGCAGTTACACTGATTTTACCTCCTGCAGGCACATG,2.0
GUT_GENOME000232_7,108957,GUT_GENOME044018_8,2494,left,right,False,GTATGCGGTCTGTCCGTCCAAGGCGATCGTATCATCGACAAGCGTCCCATAGGA,GTATTCGGTCTGTCCGTCCAAGGCGATCGTATCATCGACAAGCGTCCCATAGGA,1.907
GUT_GENOME057980_1,853,GUT_GENOME000963_8,156831,left,right,False,CTTAATGTGAAGGTTGATGGAAGGGAAGTAAAAGAGGAGGATATTATTTTCTGCCTCCACGGTAAGAAATTTGGTG,TTGACATTGAAGGTTGATGGAAGGGAAGTAAAAGAGGAGGATATTATTTTCTGCCTCCACGGTAAGGAATTTGGTG,1.75
GUT_GENOME057980_1,851,GUT_GENOME000963_8,198483,right,left,False,ATTTCTAAGATGGAAAGCATATTAAGGCTAACTCTAATGATGAAG,CATTAACGGATGGAAAGCATATTAAGGCTAACTCTAATGATGAAG,1.644
```

