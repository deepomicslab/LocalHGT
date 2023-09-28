#!/usr/bin/env python3

import os
import sys
import argparse

class Accept_Parameters:

    def __init__(self, options):

        self.reference = options.r
        self.fq1 = options.fq1
        self.fq2 = options.fq2
        self.sample_ID = options.s
        self.outdir = options.o
        self.shell_script = os.path.dirname(sys.argv[0]) + "/pipeline.sh"
        self.hit_ratio = options.hit_ratio
        self.match_ratio = options.match_ratio
        self.run_order = ''
        self.k = options.k
        self.threads = options.t

    def get_order(self):
        self.run_order = f"bash {self.shell_script} {self.reference} {self.fq1} {self.fq2} {self.sample_ID} {self.outdir} {self.hit_ratio} {self.match_ratio} {self.threads} {self.k} {options.max_peak}"
        print ("Running command:")
        print (self.run_order)

    def run(self):
        os.system(self.run_order)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Detect HGT breakpoints from metagenomics sequencing data.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-r", type=str, help="<str> Reference file.", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")

    optional.add_argument("-k", type=int, default=32, help="<int> kmer size", metavar="\b")
    optional.add_argument("-t", type=int, default=10, help="<int> number of threads", metavar="\b")
    optional.add_argument("--hit_ratio", type=float, default=0.1, help="<float> Minimum approximate kmer match ratio to extract a reference fragment.", metavar="\b")
    optional.add_argument("--match_ratio", type=float, default=0.08, help="<float> Minimum exact kmer match ratio to extract a reference fragment.", metavar="\b")
    optional.add_argument("--max_peak", type=int, default=300000000, help="<int> Maximum candidate BKP count.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv)==1:
        print (f"see python {sys.argv[0]} -h")
    else:
        acc_pa = Accept_Parameters(options)
        acc_pa.get_order()
        acc_pa.run()
