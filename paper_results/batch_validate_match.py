import sys
import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from pyfaidx import Fasta
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import DNA
import re
from Bio import pairwise2
from Bio.Seq import Seq
import pandas as pd
import subprocess




def read_event():
    hgt_event_dict = {}

    for line in open(identified_hgt):

        array = line.strip().split(',')
        if array[1] == "sample":
            continue
        sra_id = array[1]
        if sra_id not in hgt_event_dict:
            hgt_event_dict[sra_id] = []
        array[3] = int(array[3])
        array[5] = int(array[5])
        array[6] = int(array[6])
        hgt_event_dict[sra_id].append(array[2:])
    return hgt_event_dict

    
if __name__ == "__main__":

    # database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"
    # workdir = "/mnt/d/HGT/time_lines/"
    # meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    # data_pair = "/mnt/d/HGT/time_lines/SRP366030.ngs_tgs_pair.csv"
    # design_file = "/mnt/d/HGT/time_lines/sample_design.tsv"
    # result_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    # identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    # tgs_bam_dir = "/mnt/d/HGT/time_lines/tgs_bam_results"
    # verified_result = "/mnt/d/HGT/time_lines/SRP366030.verified_event.csv"

    database = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/reference/UHGG_reference.formate.fna"
    workdir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/"
    meta_data = "//mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/SRP366030.csv.txt"
    data_pair = "/mnt/disk2_workspace/wangshuai/00.strain/32.BFB/SRP366030.ngs_tgs_pair.csv"
    design_file = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid//sample_design.tsv"
    result_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/hgt/result/"
    identified_hgt = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/match/SRP366030.identified_event.csv"
    tgs_bam_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/nanopore_alignment/results/"
    verified_result = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/match/SRP366030.verified_event.csv"


    tmp_bam = workdir + "/tmp.bam"
    tmp_fastq = workdir + "/tmp.fastq"
    tmp_nano_str = workdir + "/tmp.fasta"
    reverse_tmp_ref = workdir + "/tmp.rev.fasta"
    contig_dir = workdir + "/ShastaRun/"
    contig_file = contig_dir + "/Assembly.fasta"
    standard_output = workdir + "/minimap2.log"

    hgt_event_dict = read_event()
    batch_num = 20
    h = open("work.sh", 'w')
    for i in range(batch_num):
        f = open("work_%s.sh"%(i), 'w')
        print ("nohup sh work_%s.sh &"%(i), file = h)
        f.close()
    j = 0
    for sample in hgt_event_dict:
        z = j % batch_num
        f = open("work_%s.sh"%(z), 'a')
        print ("python validate_bkp_match_for_parelle.py %s"%(sample), file = f)
        f.close()
        j += 1
    h.close()