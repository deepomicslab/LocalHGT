# LocalHGT

## Install
```
cd LocalHGT/
make
```
## Dependency
```
Python Module: scipy, numpy, pandas, sklearn, pysam, scikit-bio, Bio, pyfaidx
```
Please install samtools and bwa, and add them to system path, the version 
should be
```
samtools==1.11
bwa-0.7.17
```

## Reference database
Please contruct the database before running. The database (a single file in fasta 
formate) should contain all the representative references of your concerning envirment. 

For human gut, we can download it from the UHGG database.

Please use samtools to index the database first by

```
samtools faidx ref_database.fasta

```
LocalHGT also needs a index file, with the first running, Local will build the index the database automatically, and it will take several hours.
## Running
For each paired-end sequencing sample, with the unzipped fastq files, we can perform
LocalHGT like

```
bash pipeline.sh reference_database.fasta sample.1.fq sample.2.fq sample_name result_dir/ 0.1 0.08
```
## Output interpretion
The HGT breakpoints would be saved in the *acc.csv file.
Each line contain a pair of breakpoints, it record the two genome's ID and 
the breakpint positions. Here is a example
```
from_ref,to_ref,from_pos,to_pos,from_side,to_side,if_reverse,read_seq,ref_seq,similarity
GUT_GENOME001770_6,3495,GUT_GENOME147678_1,1797790,left,right,False,TCACCTTGCCTATATGACAGGAATCTTGCCAATCAAGAAGT,TCACCTTGCCTATATGACAGGAATCTTGCCAATCAAGAAGT,2.0
GUT_GENOME001261_6,80859,GUT_GENOME140320_24,26086,right,right,True,CCAGCGAAGATCCGGCAGGTTGTGGATATCGTCAGGAAATTATA,CGAAGATCCGGCAGGTTGTGGATATCGTCAGGAAAATGAAAGAA,1.591
GUT_GENOME001770_6,96104,GUT_GENOME268486_26,113,left,right,False,TGAAATGATATGTACCCGGTTTTCTGGACACCGGTTCAAAA,ATGATTGATATGTACCCGGTTTTCTGGACACCGGTTCAAAA,1.756
GUT_GENOME001770_6,96035,GUT_GENOME000367_1,418068,left,left,True,GCTTCAGCAGTTATCGCTGCTTTGTTCCATAAATACCCTCCTTGCTGC,GCTGCTTTGTTCCATAAATACCCTCCTTGCTGCAAGTAATAAAACCAT,1.375
GUT_GENOME001770_6,97080,GUT_GENOME111111_58,2172,left,left,True,TAGAAAGGTATGGTAATTATTATGTGGAATG,GAAAGGCATGGTAATTATTATGTGGAATGGT,1.71
```

