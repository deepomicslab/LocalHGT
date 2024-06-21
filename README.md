# LocalHGT: an ultrafast horizontal gene transfer detection method from large microbial communities

LocalHGT is a software tool which
- can detect HGT breakpoints from shotgun metagenome data
- can detect complete HGT events including the transferred sequence and related deletion and insertion sites
- can handle a human gut microbiome sample with ~two hours and <25G memory using ten threads

## Quick start
There are three methods to install LocalHGT:

### 1. conda install
First, install LocalHGT by conda via the [bioconda](https://anaconda.org/bioconda/localhgt) channel.
It is recommended to create a new conda environment and install LocalHGT simultaneously.
```
conda create --name localhgt -c bioconda -c conda-forge localhgt
conda activate localhgt
```
or use mamba
```
mamba create --name localhgt -c bioconda -c conda-forge localhgt
mamba activate localhgt
```

Also, you can install LocalHGT in existing conda environment. Notably, conflicts may occur with existing packages in this way.
```
conda install bioconda::localhgt
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
Third, obtain the source code from github, install the dependencies (listed at [Dependencies](#dependencies)), and then install.
Ensure your platform has `c++ compiler` (g++) and `make`.

If you have the root access, just install by:
```
make
python3 setup.py install
```
Otherwise, run 
```
make
chmod 744 scripts/*
cd scripts && ln -s localhgt.py localhgt
```
Then add the full path of `scripts/` to the system path by adding `export PATH="/full-path-to/scripts/:$PATH"` to `.bashrc`. 

### Test
```
cd test/
sh run_BKP_detection.sh # test HGT breakpoint detection
sh run_event_detection.sh # test HGT event detection
```
See `output/test_sample.acc.csv` for breakpoint results, and see `test_event_output.csv` for event results.

### Run
After installation, perform LocalHGT with
```
localhgt --help
localhgt bkp --help
localhgt event --help
```
Note:
- If you meet any issues, take a look at [Bug fix](#bug-fix).
- LocalHGT only accepts paired-end shotgun metagenomic sequencing data.
- LocalHGT supports Linux and Windows WSL platforms.


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
LocalHGT requires a reference, and LocalHGT accepts the reference by parameter `localhgt bkp -r`; please give the full path of the reference to `-r`.  The reference contains the representative genome of your interested microbes, and it should be a single `fasta` file.

We have prebuilt several references for users to conveniently use. It is recommended to choose a habitat-specific reference. For example, if you want to analyze the human oral microbiome, you can use [human-oral-v1-0-1](https://doi.org/10.5281/zenodo.10959731). The prebuilt references are:
- [human-gut-v1 (UHGG v1)](https://doi.org/10.5281/zenodo.10908234)
- [human-gut-v2 (UHGG v2)](https://doi.org/10.5281/zenodo.10959929)
- [human-oral-v1-0-1](https://doi.org/10.5281/zenodo.10959731)
- [human-vaginal-v1-0](https://doi.org/10.5281/zenodo.10952065)

If you cannot visit these links, please wait a while and try it again. The related annotation information can be obtained from [MGnify](https://www.ebi.ac.uk/metagenomics/browse/genomes). 

Moreover, [ProGenomes3](https://progenomes.embl.de/) provides several representative genome sets. Please click the below link to download it, and then unzip it. Subsequently, pass its path to `-r`. 
- [all microbial representative genomes ](https://progenomes.embl.de/data/repGenomes/progenomes3.contigs.representatives.fasta.bz2)
- [Aquatic](https://progenomes.embl.de/data/habitats/representatives.aquatic.contigs.fasta.gz)
- [Disease associated](https://progenomes.embl.de/data/habitats/representatives.disease_associated.contigs.fasta.gz)
- [Food associated](https://progenomes.embl.de/data/habitats/representatives.food_associated.contigs.fasta.gz)
- [Freshwater](https://progenomes.embl.de/data/habitats/representatives.freshwater.contigs.fasta.gz)
- [Host associated](https://progenomes.embl.de/data/habitats/representatives.host_associated.contigs.fasta.gz)
- [Host plant associated](https://progenomes.embl.de/data/habitats/representatives.host_plant_associated.contigs.fasta.gz)
- [Sediment mud](https://progenomes.embl.de/data/habitats/representatives.sediment_mud.contigs.fasta.gz)
- [Soil](https://progenomes.embl.de/data/habitats/representatives.soil.contigs.fasta.gz)


Note:
- the reference file should be uncompressed.
- given a new reference, `LocalHGT` will index it automatically, and it will take a few hours (e.g., several hours for UHGG v1). 
- reference index file size is approx (reference size) * 4 * (number of hash functions), make sure the disk has enough space. The number of hash functions is defaulted by 3, and is denoted by the parameter `-e`.


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
The detected HGT breakpoints are stored in the `<sample name>.acc.csv` file within the output folder. See [how to interpret results](#hgt-breakpoints).

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
              representative references of concerned bacteria. It should be
              the same as the reference file used in localhgt bkp -r. (default: None)
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
See [how to interpret results](#hgt-events).

Note:
- the reference file (given by `-r`) should be the same as the reference file used in `localhgt bkp` (also given by `-r`).
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
| reverse_flag  | if the transferred sequence is reversely inserted into recipient  |


## Dependencies
### Python modules:
```
python>=3.7.12
scikit-bio=0.5.6 (install this module will install scikit-learn, scipy, numpy, and pandas simultaneously)
networkx=2.6.3
typing-extensions>=4.11.0
biopython
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
The above tools should be installed in the system path. While the installation process using conda can be expedited by not specifying the version, it is important to note that this approach may potentially introduce unexpected bugs.

## Bug fix
- If you meet `No module named _sysconfigdata_x86_64_conda_cos7_linux_gnu`, according to the [solution](https://stackoverflow.com/questions/68261254/conda-error-sysconfigdata-x86-64-conda-linux-gnu), just run 
```
cp ${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos6_linux_gnu.py ${CONDA_PREFIX}/lib/python3.7/_sysconfigdata_x86_64_conda_cos7_linux_gnu.py
```
- If you meet `cannot import name TypeAlias from typing_extensions`, according to the [solution](https://github.com/alexdelorenzo/cast_control/issues/16), run
```
pip install typing-extensions --upgrade
```
- If you meet `ERROR: Could not build wheels for scikit-bio, hdmedians, which is required to install pyproject.toml-based projects` or `Python.h: No such file or directory`, according to the [solution](https://stackoverflow.com/questions/21530577/fatal-error-python-h-no-such-file-or-directory), run 
```
sudo apt-get install python3-dev
```
- If you meet `UserWarning: Signature b'\x00\xd0\xcc\xcc\xcc\xcc\xcc\xcc\xfb\xbf\x00\x00\x00\x00\x00\x00' for <class 'numpy.longdouble'> does not match any known type: falling back to type probe function.`, according to the [solution](https://github.com/nipy/nibabel/issues/1309), run
```
conda install numpy==1.24
```
- If you meet `urllib3 v2.0 only supports OpenSSL 1.1.1+, currently the 'ssl' module is compiled with LibreSSL 2.8.3.`, according to the [solution](https://stackoverflow.com/questions/76187256/importerror-urllib3-v2-0-only-supports-openssl-1-1-1-currently-the-ssl-modu), run 
```
conda install urllib3<2.0
```
- If you meet `No module named 'pexpect'`, `No module named 'decorator'`, or `No module named 'cachecontrol'`, run
```
conda install scikit-bio==0.5.6
```

## Paper results generation
The scripts to produce the results in the paper can be found in `paper_results/*`, where some scripts might require additional python modules. HGT breakpoints and events detected by LocalHGT from all real samples can be seen at [Detected HGTs in real data](https://doi.org/10.5281/zenodo.10906354).

## Citation
Shuai Wang, Yiqi Jiang, Lijia Che, Ruo Han Wang, Shuai Cheng Li, Enhancing insights into diseases through horizontal gene transfer event detection from gut microbiome, Nucleic Acids Research, 2024;, gkae515, [https://doi.org/10.1093/nar/gkae515](https://doi.org/10.1093/nar/gkae515)

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).

