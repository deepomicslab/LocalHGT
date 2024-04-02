# LocalHGT: an ultrafast horizontal gene transfer detection method from large microbial communities

## Install
```
git clone https://github.com/deepomicslab/LocalHGT.git --depth 1
cd LocalHGT/
conda env create --name localhgt -f environment.yml
conda activate localhgt
make

python scripts/infer_HGT_breakpoint.py -h  # detect HGT breakpoints
python scripts/infer_HGT_event.py -h # detect complete HGT events
```

## Test
```
cd test/
sh run_BKP_detection.sh
sh run_event_detection.sh
```
See `output/test_sample.acc.csv` for breakpoint results, and see `test_event_output.csv` for event results.

## Construct reference database
Users can construct customized reference databases. The reference database (a single `fasta` file) should contain all the representative genomes of your concerning bacteria. 
For example, to analyze HGT events between bacteria A, B, and C, we should collect the representative genomes of A, B, and C and merge them into a single fasta file.

To analyze HGTs in the human gut microbiome, we use gut-specific representative genome collection of the [UHGG](https://www.nature.com/articles/s41587-020-0603-3) database. 
The script `paper_results/build_UHGG_reference.py` can download the human gut-specific UHGG v1 database. The command is
```
usage: python build_UHGG_reference.py -h

Construct the reference from the gut-specific UHGG V1 database.

optional arguments:
  -r        <str> Generated reference file. (default: uhgg_v1.rep.fasta)
  -b        <str> Folder saves all the downloaded assemblies. (default:
              genomes/)
  -m        <int> Try this number of times until all the genomes are
              downloaded. (default: 10)
  -h, --help

Example: python paper_results/build_UHGG_reference.py -r my_ref.fasta -b genomes_dir -m 4
```

Then index the database using `samtools`
```
samtools faidx ref_database.fasta
```
Also, at the first run, `LocalHGT` will index the database automatically, which will take several hours. 

Note:
- reference index file size is approx (reference size) * 4 * (number of denoted hash functions), make sure the disk has enough space.


## Run

### Refine sequencing reads
It is highly advised to refine the sequencing reads using tools like [fastp](https://github.com/OpenGene/fastp) before calling HGT.
We include `fastp` in the conda environment, and refine the reads like
```
fastp -i raw_1.fq -I raw_2.fq -o refine_1.fq -O refine_2.fq
```

Note:
- Currently, LocalHGT only supports paired-end short-read sequencing data (e.g., Illumina data).

### Detect HGT breakpoints
First, infer HGT breakpoints by running `python scripts/infer_HGT_breakpoint.py` like
```
usage: python infer_HGT_breakpoint.py -h

Detect HGT breakpoints from metagenomics sequencing data.

required arguments:
  -r             <str> reference file which contains all the representative
                   references of concerned bacteria. (default: None)
  --fq1          <str> unzipped fastq 1 file. (default: None)
  --fq2          <str> unzipped fastq 2 file. (default: None)
  -s             <str> Sample name. (default: sample)
  -o             <str> Output folder. (default: ./)

optional arguments:
  -k             <int> kmer length. (default: 32)
  -t             <int> number of threads. (default: 10)
  -e             <int> number of hash functions (1-9). (default: 3)
  -a             <0/1> 1 indicates retain reads with XA tag. (default: 1)
  -q             <int> minimum read mapping quality in BAM. (default: 20)
  --seed         <int> seed to initialize a pseudorandom number generator.
                   (default: 1)
  --use_kmer     <1/0> 1 means using kmer to extract HGT-related segment, 0
                   means using original reference. (default: 1)
  --hit_ratio    <float> minimum fuzzy kmer match ratio to extract a
                   reference fragment. (default: 0.1)
  --match_ratio  <float> minimum exact kmer match ratio to extract a
                   reference fragment. (default: 0.08)
  --max_peak     <int> maximum candidate BKP count. (default: 300000000)
  --sample       <float> down-sample in kmer counting: (0-1) means sampling
                   proportion, (>1) means sampling base count (bp). (default:
                   2000000000)
  --read_info    <0/1> 1 indicates including reads info, 0 indicates not
                   (just for evaluation). (default: 1)
  -h, --help
```
A command example:
```
python scripts/infer_HGT_breakpoint.py -r reference.fa --fq1 species20_snp0.01_depth50_reads100_sample_0.1.fq --fq2 species20_snp0.01_depth50_reads100_sample_0.2.fq -s species20_snp0.01_depth50_reads100_sample_0  -o test
```
The detected HGT breakpoints are stored in the `<sample name>.acc.csv` file within the output folder.

Note:
- With a small reference, we can skip extracting HGT-related segments by setting `--use_kmer 0`.
- With a small reference, while maintaining the extraction of HGT-related segments, we can set a small value of `-k` to reduce memory usage. 


### Detect complete HGT events
Second, infer complete HGT events by matching breakpoints after detecting HGT breakpoints for all the samples.
```
usage: python infer_HGT_event.py -h

Infer complete HGT events by matching breakpoint pairs.

required arguments:
  -r        <str> Reference file. (default: None)
  -b        <str> Folder saves all the breakpoint results from all samples.
              (default: None)
  -f        <str> Output file to save all inferred HGT events. (default:
              complete_HGT_event.csv)

optional arguments:
  -n        <int> minimum supporting split read number (default: 2)
  -m        <int> minimum transfer sequence length (default: 500)
  -h, --help
```

A command example:
```
python scripts/infer_HGT_event.py -r ref.fa -b output -f test_event_output.csv
```
The `output` folder contain all the HGT breakpoint results generated by `infer_HGT_breakpoint.py`.


Note:
- It is recommended to detect HGT breakpoints for each sample and store the results in a common output folder. Subsequently, when detecting complete HGT events, specify the output folder using the `-b` parameter. This approach allows LocalHGT to consider all the samples collectively, resulting in more reliable results for complete HGT events.

## Output interpretation

###  HGT breakpoints
The HGT breakpoints are saved in the `*acc.csv` file. Here is an example:
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
| if_reverse  | if the transferred gene is reversely inserted to receptor  |
| read_seq  | used reads sequence in local alignment  |
| ref_seq  | aligned reference sequence in local alignment  |
| similarity  | alignment similarity in local alignment  |
| from_split_reads  | number of split reads mapped to the from_pos  |
| to_split_reads  | number of split reads mapped to the to_pos  |
| cross_split_reads  | number of split reads supported the breakpoint pair  |
|pair_end| number of paired-end reads supported the breakpoint pair     |

Note:
- Before conducting the analysis, it is advisable to filter out any instances of gene transfer between different contigs of the same species that may be present in the results.
For example, GUT_GENOME000031_1 and GUT_GENOME000031_2 belong to the same species, the detected HGT breakpoint pairs or HGT events between them should be discarded.

###  HGT events
HGT event results are like
```
sample,receptor,insert_locus,donor,delete_start,delete_end,reverse_flag
species20_snp0,GUT_GENOME000149_2,369882,GUT_GENOME000537_12,40443,91965,True
species20_snp0,GUT_GENOME000189_7,22927,GUT_GENOME000534_3,78150,79127,True
species20_snp0,GUT_GENOME001261_3,336654,GUT_GENOME000202_1,72162,123917,False
```
Interpret each column as:
| Header  | Description |
| :-------------:| :-------------: |
| sample  | Sample ID  |
| receptor  | receptor genome  |
| insert_locus  | the transfer sequence is inserted at this locus |
| donor  | donor genome  |
| delete_start  | the start site of the transfer sequence on the donor genome  |
| delete_end  | the end site of the transfer sequence on the donor genome  |
| reverse_flag  | if the transferred gene is reversely inserted to receptor  |


## Dependency
We recommend constructing the environment using `conda` with the `environment.yml`.
Users can also prepare the env as follows: 
```
Python 3.7+
Python Module: scipy, numpy, pandas, sklearn, pysam, scikit-bio, biopython, pyfaidx
```
Please install `samtools` and `bwa`, and add them to the system path, the version 
should be
```
samtools==1.11+
bwa-0.7.17+
```

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).

