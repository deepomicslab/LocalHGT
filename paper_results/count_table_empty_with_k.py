#!/usr/bin/env python3


"""
perform k-mer-counting with different length of k-mer, based on the high complexity sample, 
generated by combining pure data and CAMI high complexity sample.
"""

import os


class Order():
    def __init__(self):
        self.all_k = [x for x in range(20, 33)]
        self.all_ratio = [13, 23, 33, 50, 75, 100]

    def generate(self, k, ratio):
        com = f"./count_diff_kmer {fq1} {fq2}  {k} {ratio}"
        return com

    def run(self):
        for k in self.all_k:
            for ratio in self.all_ratio:
                com = self.generate(k, ratio)
                print (com)
                os.system(com)

fq1 = "/mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.01_depth50_reads150_sample_0_high.1.fq"
fq2 = "/mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.01_depth50_reads150_sample_0_high.2.fq"

os.system('g++ -pthread -o count_diff_kmer src/count_diff_kmer.cpp')
order = Order()
# order.run()
os.system(order.generate(32, 13))

