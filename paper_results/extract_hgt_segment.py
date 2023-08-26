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
        index = array[0]
        if sra_id not in hgt_event_dict:
            hgt_event_dict[sra_id] = []
        array[3] = int(array[3])
        array[5] = int(array[5])
        array[6] = int(array[6])
        transfer_seq = ref_fasta[array[4]][array[5]:array[6]].seq
        if array[7] == "True":
            transfer_seq = reverse_complement_seq(transfer_seq)
        # print (">%s:%s-%s&%s "%(array[4], array[5], array[6], index) + line.strip(), file = f)
        print (">" + line.strip(), file = f)
        print (transfer_seq, file = f)

        hgt_event_dict[sra_id].append(array[2:])
    f.close()
    return hgt_event_dict

def get_event_index():
    event_index_dict = {}

    for line in open(identified_hgt):
        array = line.strip().split(',')
        # if array[1] == "sample":
        #     continue
        sra_id = array[1]
        index = array[0]
        event_index_dict[index] = line.strip()
    return event_index_dict
        
def read_blastn(blastn_file, fasta_file):
    contig_lengths = get_contig_lengths(fasta_file)
    # print (contig_lengths)
    map_cds = set()
    # Open the blastn results file
    with open(blastn_file, "r") as f:
        # Loop over each line in the file
        for line in f:
            # Split the line into fields
            fields = line.strip().split("\t")
            # Calculate the alignment length and proportion
            align_len = int(fields[3])
            query_name = fields[0]
            query_len = contig_lengths[query_name]
            align_prop = align_len / query_len
            # Check if the alignment proportion is larger than 0.5
            if align_prop > 0.5:
                # Extract the query name
                query_name = fields[0]
                # Print the query name
                # print(query_name)
                map_cds.add(query_name)
    print (len(map_cds), len(map_cds)/len(contig_lengths))
    return len(map_cds), len(contig_lengths)

def blast_main(fasta_file, db):
    command = f"blastn -db {db} -query {fasta_file} -outfmt 6 -out {fasta_file}.out -num_threads 8"
    # print (command)
    os.system(command)
    # mapped_num, all_num = read_blastn(f"{fasta_file}.out", fasta_file)
    # return mapped_num, all_num

def get_contig_lengths(fasta_file):   
    # Initialize an empty dictionary to store the contig lengths
    contig_lengths = {}
    
    # Parse the fasta file and loop over the records
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Store the contig ID and length in the dictionary
        contig_lengths[record.id] = len(record.seq)
    
    return contig_lengths


if __name__ == "__main__":

    database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"
    phage_db = "/mnt/d/HGT/seq_ana/BlastDB/TemPhD.fasta"
    plasmid_db = "/mnt/d/HGT/seq_ana/database/plsdb.fna"

    # identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    # identified_hgt_seq = "/mnt/d/HGT/time_lines/SRP366030.identified_event.fasta"

    identified_hgt = "/mnt/d/HGT/seq_ana/identified_event.csv"
    identified_hgt_seq = "/mnt/d/HGT/seq_ana/identified_event.fasta"

    # ref_fasta = Fasta(database)
    # hgt_event_dict = read_event()

    # event_index_dict = get_event_index()

    blast_main(identified_hgt_seq, phage_db)

