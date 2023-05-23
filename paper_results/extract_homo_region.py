#!/usr/bin/env python3

import os
import numpy as np
import random
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import lzma
import csv
import scipy.special as sc
from scipy.special import comb



def extract_uniq_region(map_ref):
    genome_homo_dict = {}

    line_num = 0
    # for line in lzma.open(map_ref, mode='rt', encoding='utf-8'):
    for line in open(map_ref):
        line_num += 1
        if line_num % 1000000 == 0:
            print (line_num)
        # if line_num > 1000000:
        #     break
        array = line.split()
        if array[0] == array[1]:
            continue
        if float(array[3]) < 50: # length
            continue
        if float(array[2]) < 99: # identity
            continue
        clean_genome_1 =  "_".join(array[0].split("_")[:-1])
        clean_genome_2 =  "_".join(array[1].split("_")[:-1])
        raw_genome_2 = array[1]

        #record scaffold name
        if clean_genome_1 not in genome_homo_dict.keys():
            genome_homo_dict[clean_genome_1] = {}

        if raw_genome_2 not in genome_homo_dict[clean_genome_1]:
            genome_homo_dict[clean_genome_1][raw_genome_2] = []
        
        if float(array[8]) < float(array[9]):
            sec_s = float(array[8])
            sec_e = float(array[9])
        else:
            sec_s = float(array[9])
            sec_e = float(array[8])   
        genome_homo_dict[clean_genome_1][raw_genome_2].append([sec_s, sec_e])
    return genome_homo_dict


def is_position_in_intervals(position, intervals):
    for interval in intervals:
        if position >= interval[0] and position <= interval[1]:
            return True
    return False

        
class Acc_Bkp(object):

    def __init__(self, list):

        self.from_ref = list[0]
        self.from_bkp = int(list[1])
        self.from_side = list[2]
        self.from_strand = list[3]
        self.from_split_reads = int(list[12])

        self.to_ref = list[4]
        self.to_bkp = int(list[5])
        self.to_side = list[6]
        self.to_strand = list[7]
        self.to_split_reads = int(list[13])

        self.if_reverse = list[8]
        
        self.cross_split_reads = int(list[14])
        self.pair_end = int(list[15])

        
        self.from_ref_genome = "_".join(self.from_ref.split("_")[:-1])
        self.to_ref_genome = "_".join(self.to_ref.split("_")[:-1])
        self.score = float(list[11])


def get_split_reads_cutoff(g): # g is the number of reads in the specific sample
    # prior_a = 25 # prior probability
    prior_b = 50000000#63333330  # prior probability
    given_n = 42648185 # mean number of reads among samples
    # # n = 182534663 # max number of reads 
    # given_r = 2 # the cutoff with n reads  4
    prior_a = 20 # prior probability
    given_r = 2 # the cutoff with n reads  4

    for m in range(100):
        alpha = prior_a + m
        beta= prior_b + g - m
        p = 1
        for k in range(given_r + 1):
            choose = comb(given_n, k)
            beta1 = sc.beta(alpha + k, beta + given_n - k)
            beta2 = sc.beta(alpha, beta)
            # print (choose, beta1)
            p -= choose * beta1 / beta2
        if p > 0.9: #0.9
            break
    return m, p

def read_bkp(bkp_file, filtered_file):
    f = open(bkp_file)
    all_rows = csv.reader(f)

    out = open(filtered_file, 'w')
    csv_writer = csv.writer(out)
    final_row_num = 0

    for row in all_rows:
        if row[0][0] == "#":
            reads_num = int(row[0].split(";")[0].split(":")[1])
            min_split_num, p = get_split_reads_cutoff(reads_num)
            # print (bkp_file, reads_num, "cutoff", min_split_num)
        elif row[0] == "from_ref":
            pass
        else:
            # print (row)
            if len(row) < 13:
                continue
            eb = Acc_Bkp(row)
            if eb.cross_split_reads < min_split_num:
                continue
            if eb.from_ref_genome in genome_homo_dict:
                if eb.to_ref in genome_homo_dict[eb.from_ref_genome]:
                    intervals = genome_homo_dict[eb.from_ref_genome][eb.to_ref]
                    if is_position_in_intervals(eb.to_bkp, intervals):
                        continue
            if eb.to_ref_genome in genome_homo_dict:
                if eb.from_ref in genome_homo_dict[eb.to_ref_genome]:
                    intervals = genome_homo_dict[eb.to_ref_genome][eb.from_ref]
                    if is_position_in_intervals(eb.from_bkp, intervals):
                        continue
            final_row_num += 1
        csv_writer.writerow(row)

    f.close()
    out.close()
    if final_row_num == 0:
        os.system("rm %s"%(filtered_file))

def main(result_dir, sample_num):
    files = os.listdir(result_dir)
    # print (len(files))
    for acc_file in files:
        if not re.search("acc.csv", acc_file):
            continue
        if re.search("repeat.acc.csv", acc_file):
            continue
        # if hgt_result_dir + acc_file != "/mnt/d/breakpoints/script/analysis/hgt_results/CCIS98832363ST-4-0.acc.csv":
        #     continue
        # print (hgt_result_dir + acc_file)
        read_bkp(result_dir + acc_file, filter_hgt_result_dir + acc_file)
        sample_num += 1
        # break
    return sample_num



if __name__ == "__main__":
    blast_file = '/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.blast.out'
    genome_homo_dict_file = "/mnt/d/HGT/time_lines/ref_repeat.npy"

    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    tgs_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    wenkui_dir = "/mnt/d/breakpoints/HGT/CRC/wenkui/"

    filter_hgt_result_dir = "/mnt/d/breakpoints/script/analysis/homo_filter/"

    if os.path.isfile(genome_homo_dict_file):
        genome_homo_dict = np.load(genome_homo_dict_file, allow_pickle='TRUE').item()
        print ("loaded")
    else:
        genome_homo_dict = extract_uniq_region(blast_file)
        np.save(genome_homo_dict_file, genome_homo_dict)

    # genome_homo_dict = {}

    # intervals = genome_homo_dict["GUT_GENOME123416"]["GUT_GENOME147678_1"]
    # print (is_position_in_intervals(1690584, intervals), is_position_in_intervals(1958707, intervals))

    main(tgs_dir, 0)
    # main(hgt_result_dir, 0)
    # main(wenkui_dir, 0)