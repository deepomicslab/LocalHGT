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


def reverse_complement_seq(sequence):
    sequence = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

def read_event():
    hgt_event_dict = {}
    f = open(identified_hgt_seq, 'w')
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
        transfer_seq = ref_fasta[array[4]][array[5]:array[6]].seq
        if array[7] == "True":
            transfer_seq = reverse_complement_seq(transfer_seq)
        print (">%s:%s-%s "%(array[4], array[5], array[6])+ line.strip(), file = f)
        print (transfer_seq, file = f)

        hgt_event_dict[sra_id].append(array[2:])
    f.close()
    return hgt_event_dict


if __name__ == "__main__":

    database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"
    identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    identified_hgt_seq = "/mnt/d/HGT/time_lines/SRP366030.identified_event.fasta"

    ref_fasta = Fasta(database)
    hgt_event_dict = read_event()

