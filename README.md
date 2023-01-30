# LocalHGT: an ultrafast horizontal gene transfer detection method from large microbial communities

## Install
```
git clone https://github.com/deepomicslab/LocalHGT.git --depth 1
cd LocalHGT/
conda env create --name localhgt -f environment.yml
conda activate localhgt
sudo chmod 744 ./extract_ref
python main.py -h
```

## Reference database
The database (a single `fasta` file) should contain all the representative references of your concerning bacteria. 

For example, use gut-specific representative genomes collection of the [UHGG](https://www.nature.com/articles/s41587-020-0603-3) database. 
The script `paper_results/build_UHGG_reference.py` can download the gut-specific UHGG v1 database.

Then index the database. First, build `samtools` index 
```
samtools faidx ref_database.fasta
```
Second, at the first running, `LocalHGT` will index the database automatically, and it will take several hours.

## Test
```
cd test/
sh run.sh
```
See `test/output/` for results. 

## Running
Run `python main.py` like
```
usage: main.py -h

Detect HGT events from metagenomics sequencing data.

required arguments:
  -r             <str> Reference file. (default: None)
  --fq1          <str> unzipped fastq 1 file. (default: None)
  --fq2          <str> unzipped fastq 2 file. (default: None)
  -s             <str> Sample name. (default: sample)
  -o             <str> Output folder. (default: ./)

optional arguments:
  -t             <int> number of threads (default: 5)
  --hit_ratio    <float> Minimum approximate kmer match ratio to extract a
                   reference fragment. (default: 0.1)
  --match_ratio  <float> Minimum exact kmer match ratio to extract a
                   reference fragment. (default: 0.08)
  -h, --help
```
Or, with the unzipped paired-end fastq files, run
```
bash pipeline.sh reference_database.fasta sample.1.fq sample.2.fq sample_name result_dir/ 0.1 0.08
```
A command example:
```
python ./main.py -r reference.fa --fq1 /mnt/d/breakpoints/HGT/uhgg_length//species20_snp0.01_depth50_reads100_sample_0.1.fq --fq2 /mnt/d/breakpoints/HGT/uhgg_length//species20_snp0.01_depth50_reads100_sample_0.2.fq -s species20_snp0.01_depth50_reads100_sample_0  -o test
```

## Output interpretion
The HGT breakpoints would be saved in the `*acc.csv` file. Here is an example:
```
# the number of reads in the sample is: 41723899; Insert size is 681.
from_ref,from_pos,from_side,from_strand,to_ref,to_pos,to_side,to_strand,if_reverse,read_seq,ref_seq,similarity,from_split_reads,to_split_reads,cross_split_reads,pair_end
GUT_GENOME000031_3,41994,head,+,GUT_GENOME096518_4,692224,tail,+,False,GTGTCGGGGCTTATGATAATCATATCTTATTTTTC,GTGTCGGGGCTTATGATAATCATATCTTTTTTTTC,1.857,3,2,2,5
GUT_GENOME096518_4,725079,head,+,GUT_GENOME000031_3,41992,tail,+,False,TACGCGGAGGGATTATGGGAATGCTCACGGCAATCGAAATGGGAA,CCCGCGGCGGGATTATGGGAATGCTCACGGCAATCGAAATGGGAA,1.8,1,3,1,5
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
We recommand construct the environment using `conda` with the `environment.yml`.
Users can also prepare the env as follows: 
```
Python 3.7+
Python Module: scipy, numpy, pandas, sklearn, pysam, scikit-bio, biopython, pyfaidx
```
Please install `samtools` and `bwa`, and add them to system path, the version 
should be
```
samtools==1.11+
bwa-0.7.17+
```

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).

