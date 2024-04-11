# LocalHGT: an ultrafast horizontal gene transfer detection method from large microbial communities

## Install
There are three methods to install LocalHGT:

### 1. conda install
First, install it by conda via the bioconda channel
```
conda install localhgt
```
Notably, it is recommanded to create a new conda environment and then install LocalHGT
```
conda create --name localhgt --channel wshuai294 localhgt
conda activate localhgt
```

### 2. source code install with Conda environment
Second, obtain the source code, construct the environment with conda, and then install
```
git clone https://github.com/deepomicslab/LocalHGT.git --depth 1
cd LocalHGT/
conda env create --name localhgt -f brief_env.yml
conda activate localhgt
sh install.sh
```

### 3. source code install with self-build environment
Third, obtain the source code from github, install the dependencies by yourself (which are listed at the bottom in `Dependencies` section), and then install
```
make
python setup.py install
```
Ensure your platform has `c++ compiler` (g++/clang++) and `make`.

### Run
After installation, perform LocalHGT with
```
localhgt --help
localhgt bkp --help
localhgt event --help
```
Note:
- LocalHGT only accept paired-end shotgun metagenomic sequencing data.
- LocalHGT supports Linux and MacOS platforms.
- If you meet the issue: `No module named _sysconfigdata_x86_64_conda_cos7_linux_gnu`, just run 
```
cp ${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos6_linux_gnu.py ${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos7_linux_gnu.py
```
which is referred to the [solution](https://stackoverflow.com/questions/68261254/conda-error-sysconfigdata-x86-64-conda-linux-gnu).

## Test
```
cd test/
sh run_BKP_detection.sh # test HGT breakpoint detection
sh run_event_detection.sh # test HGT event detection
```
See `output/test_sample.acc.csv` for breakpoint results, and see `test_event_output.csv` for event results.


## Basic Usage 
### Main functions
```
usage: localhgt [-h] {bkp,event} ...

LocalHGT: an ultrafast HGT detection method from large microbial communities

optional arguments:
  -h, --help   show this help message and exit

Command:
  {bkp,event}
    bkp        Detect HGT breakpoints from metagenomic sequencing data.
    event      Infer complete HGT events based on detected HGT breakpoints.

Note: To use LocalHGT, first detect HGT breakpoints with 'localhgt bkp'. After
that, detect HGT events based on the detected HGT breakpoints with 'localhgt
event'. Detailed documentation can be found at
https://github.com/deepomicslab/LocalHGT
```


### Construct reference 
LocalHGT require a reference, which contains the representative genome of your interested bacteria. The reference should be a single `fasta` file.
We have prebuilt several references for users to conveniently use: 
- [human-gut-v1 (UHGG v1)](https://doi.org/10.5281/zenodo.10908234)
- [human-oral-v1-0-1)](xx)
- [human-vaginal-v1-0](xx)

These data and related annotation information can be obtained from [MGnify](https://www.ebi.ac.uk/metagenomics/browse/genomes).
At the first run, `LocalHGT` will index the reference automatically (e.g., it will take several hours for UHGG v1). 

Note:
- reference index file size is approx (reference size) * 4 * (number of denoted hash functions), make sure the disk has enough space.
- the reference file should be uncompressed.

### Detect HGT breakpoints
First, infer HGT breakpoints by running `localhgt bkp` like
```
usage: localhgt bkp [-r] [--fq1] [--fq2] [-s] [-o] [-k] [-t]
                    [-e] [-a] [-q] [--seed] [--use_kmer]
                    [--hit_ratio] [--match_ratio] [--max_peak]
                    [--sample] [--refine_fq] [--read_info] [-h]

Detect HGT breakpoints from metagenomic sequencing data. Example: localhgt bkp
-r reference.fa --fq1 test.1.fq --fq2 test.2.fq -s test -o outdir

required arguments:
  -r             <str> Uncompressed reference file, which contains all the
                   representative references of concerned bacteria. (default:
                   None)
  --fq1          <str> Uncompressed fastq 1 file. (default: None)
  --fq2          <str> Uncompressed fastq 2 file. (default: None)
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
  --refine_fq    <0/1> 1 indicates refine the input fastq file using fastp
                   (recommended). (default: 0)
  --read_info    <0/1> 1 indicates including reads info, 0 indicates not
                   (just for evaluation). (default: 1)
  -h, --help
```
The detected HGT breakpoints are stored in the `<sample name>.acc.csv` file within the output folder.

Note:
- The fastq files should be uncompressed.
- With a small reference, we can skip extracting HGT-related segments by setting `--use_kmer 0`.
- With a small reference, while maintaining the extraction of HGT-related segments, we can set a small value of `-k` to reduce memory usage. 


### Detect complete HGT events
Second, infer complete HGT events by running `localhgt event` like
```
usage: localhgt event [-r] [-b] [-f] [-n] [-m] [-h]

Infer complete HGT events based on detected HGT breakpoints. Example: localhgt
event -r reference.fa -b outdir -f test_event.csv

required arguments:
  -r        <str> <str> Uncompressed reference file, which contains all the
              representative references of concerned bacteria. (default: None)
  -b        <str> the folder stores all the breakpoint results from all
              samples, i.e., a folder stores all the *acc.csv files generated
              by 'localhgt bkp' (default: None)
  -f        <str> Output file to save all inferred HGT events. (default:
              complete_HGT_event.csv)

optional arguments:
  -n        <int> minimum supporting split read number (default: 2)
  -m        <int> minimum transfer sequence length (default: 500)
  -h, --help
```

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
| from_ref  | first genome ID  |
| from_pos  | first genome breakpoint position  |
| from_side  | first genome side  |
| from_strand  | first genome strand  |
| to_ref  | second genome ID  |
| to_pos  | second genome breakpoint position  |
| to_side  | second genome side  |
| to_strand  | second genome strand  |
| if_reverse  | if the transferred sequence is reversely inserted to recipient  |
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
| receptor  | recipient genome  |
| insert_locus  | the transferred sequence is inserted at this locus of the recipient |
| donor  | donor genome  |
| delete_start  | the start site of the transferred sequence on the donor genome  |
| delete_end  | the end site of the transferred sequence on the donor genome  |
| reverse_flag  | if the transferred sequence is reversely inserted to recipient  |


## Dependencies
### Python modules:
```
python>=3.7.12
scikit-bio=0.5.6
networkx=2.6.3
scikit-learn
scipy
biopython
numpy
pandas
pysam
pyfaidx
```
### Third-party software
```
samtools>=1.11
bwa>=0.7.17
fastp>=0.23.2
seqkit>=2.6.1
make
cxx-compiler
```
The above tools should be installed in the system path.

### Paper results generation
The scripts to produce the results in the paper can be found in `paper_results/*`. HGT breakpoints and events detected by LocalHGT from all real samples can be seen at [Detected HGTs in real data](https://doi.org/10.5281/zenodo.10906354).

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).

