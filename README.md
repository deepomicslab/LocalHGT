# LocalHGT: an ultrafast HGT detection method from large microbial communities

## Install
```
git clone https://github.com/deepomicslab/LocalHGT.git --depth 1
cd LocalHGT/
conda env create --prefix localhgt -f environment.yml
conda activate localhgt
make
```

## Reference database
Please construct the database before running. The database (a single file in fasta 
formate) should contain all the representative references of your concerning environment. 

For the human gut, we can download it from the Unified Human Gastrointestinal Genome 
([UHGG](https://www.nature.com/articles/s41587-020-0603-3)) database. 
The script `build_UHGG_reference.py` can be used to download the UHGG v1 genomes.

Please use samtools to index the database first by

```
samtools faidx ref_database.fasta
```
LocalHGT also needs an index file for the reference, with the first running, LocalHGT will index 
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
# the number of reads in the sample is: 41723899; Insert size is 681.
from_ref,from_pos,from_side,from_strand,to_ref,to_pos,to_side,to_strand,if_reverse,read_seq,ref_seq,similarity,from_split_reads,to_split_reads,cross_split_reads,pair_end
GUT_GENOME001745_51,3736,head,+,GUT_GENOME147854_27,1150,head,-,True,GCTGAACTAAAGGGAGTAAAACTGAAAGATTATGCAGGCACGAA,GCTGAACTAAAGGGAGTAAAACTGAAAGATTATGCAGGCACGAA,2.0,3,2,1,6371
GUT_GENOME147854_27,694,tail,-,GUT_GENOME001745_51,2933,tail,+,True,CATAACGGCACAAGAAAAGATAACCGACCTTATCGG,CATAACGGCACAAGAAAAGATAACCGACCTTATCGG,2.0,11,4,1,401
```
Interpret each column as:
| Header  | Description |
| :-------------:| :-------------: |
| from_ref  | from reference ID  |
| from_pos  | from reference breakpoint position  |
| from_side  | from reference side  |
| from_strand  | from reference strand  |
| to_ref  | to reference ID  |
| to_pos  | to reference breakpoint position  |
| to_side  | to reference side  |
| to_strand  | to reference strand  |
| if_reverse  | if the transferred gene is reverse inserted to receptor  |
| read_seq  | used reads sequence in local alignment  |
| ref_seq  | aligned reference sequence in local alignment  |
| similarity  | alignment similarity in local alignment  |
| from_split_reads  | number of split reads mapped to the from_pos  |
| to_split_reads  | number of split reads mapped to the to_pos  |
| cross_split_reads  | number of split reads supported the breakpoint pair  |
|pair_end| number of paired-end reads supported the breakpoint pair     |

## Dependency
```
Python 3.8+
Python Module: scipy, numpy, pandas, sklearn, pysam, scikit-bio, biopython, pyfaidx
```
Please install samtools and bwa, and add them to system path, the version 
should be
```
samtools==1.11+
bwa-0.7.17+
```

