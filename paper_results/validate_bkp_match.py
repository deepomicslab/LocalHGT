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


final_total = 0
final_verified = 0

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

def get_support_reads(bam, chr, pos):
    # Open the input BAM file
    bamfile = pysam.AlignmentFile(bam, "rb")
    read_list = []
    # # Create a new BAM file for the supporting reads
    # out_bamfile = pysam.AlignmentFile("supporting_reads.bam", "wb", template=bamfile)

    # Iterate over the reads in the BAM file
    for read in bamfile.fetch(chr, pos-1, pos+1):

        if read.is_unmapped:
            continue

        # Check if the read overlaps the breakpoint position
        if read.reference_start <= pos and read.reference_end >= pos:

            # print(read.query_name, read.reference_name, read.reference_start, read.reference_end, read.cigar, sep="\t")
            # Write the read to the output BAM file
            # out_bamfile.write(read)
            read_list.append(read)

    # Close the input and output BAM files
    bamfile.close()
    # out_bamfile.close()
    return read_list


def read_meta():
    
    ngs_sample_dict = {}
    sample_tgs_dict = {}
    ngs_tgs_pair = {}

    for line in open(meta_data):
        if line.strip() == '':
            continue
        array = line.strip().split(',')
        if array[0] != 'Run':
            sra_id = array[0]
            sample_id = array[-2]
            if re.search("_", sample_id):
                sample_id = sample_id.split("_")[1]

            if re.search("PAIRED", line):
                ngs_sample_dict[sra_id] = sample_id
            elif re.search("SINGLE", line):
                sample_tgs_dict[sample_id] = sra_id
    data = []
    for ngs_sra in ngs_sample_dict:
        sample_id = ngs_sample_dict[ngs_sra]
        tgs_sra = sample_tgs_dict[sample_id]
        ngs_tgs_pair[ngs_sra] = tgs_sra
        
        data.append([sample_id, ngs_sra, tgs_sra])
        
    # df = pd.DataFrame(data, columns = ["sample", "ngs", "nanopore"])
    # df.to_csv(data_pair, sep=',')

    return ngs_tgs_pair       
        

class Map():

    def __init__(self):
        self.ref = database 
        self.ref_fasta = Fasta(self.ref)
        self.window = 200
        self.max_length = 10000

    def extract_ref_seq(self, scaffold_name, start, end):
        # self.ref_fasta = Fasta(self.ref)
        if start < 1:
            start = 1
        # print (scaffold_name, start, end)
        return self.ref_fasta[scaffold_name][start:end].seq

    def get_reverse_complement_seq(self, sequence):
        sequence = sequence[::-1]
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
        string = sequence.translate(trantab)
        return string

    def get_reads(self, hgt_event):
        insert_support_reads = get_support_reads(bamfile, hgt_event[0], hgt_event[1])
        delete_support_reads_1 = get_support_reads(bamfile, hgt_event[2], hgt_event[3])
        delete_support_reads_2 = get_support_reads(bamfile, hgt_event[2], hgt_event[4])
        reads_set = insert_support_reads + delete_support_reads_1 + delete_support_reads_2

        # templatefile = pysam.AlignmentFile(bamfile, "rb")
        # out_bamfile = pysam.AlignmentFile(tmp_bam, "wb", template=templatefile)
        # uniq_dict = {}
        # for read in reads_set:
        #     if read.query_name not in uniq_dict:
        #         out_bamfile.write(read)
        #     uniq_dict[read.query_name] = 1
        # out_bamfile.close()
        # os.system(f"samtools bam2fq {tmp_bam} > {tmp_fastq}")
        # templatefile.close()

        return reads_set

    def for_each_sample(self, sample):
        if sample not in hgt_event_dict:
            print ("No hgt for sample", sample)
        else:
            print ("Detect hgt for sample", sample)

        total_hgt_event_num = 0
        valid_hgt_event_num = 0

        for hgt_event in sorted(list(hgt_event_dict[sample])):
            # if hgt_event[0] != "GUT_GENOME095993_10":
            #     continue
            # if 'GUT_GENOME156655_37' not in hgt_event:
            #     continue
            verify_flag = False
            reads_set = self.get_reads(hgt_event)

            delete_len = hgt_event[4] - hgt_event[3]
            
            if  delete_len < self.max_length:
                insert_left = self.extract_ref_seq(hgt_event[0], hgt_event[1] - self.window, hgt_event[1])
                insert_right = self.extract_ref_seq(hgt_event[0], hgt_event[1], hgt_event[1] + self.window)
                delete_seq = self.extract_ref_seq(hgt_event[2], hgt_event[3], hgt_event[4])
                reverse_delete_seq = self.get_reverse_complement_seq(delete_seq)
                merged_seq = insert_left + delete_seq + insert_right
                rev_merged_seq = insert_left + reverse_delete_seq + insert_right
                verify_flag = self.verify(reads_set, merged_seq, rev_merged_seq)
            else:

                delete_seq = self.extract_ref_seq(hgt_event[2], hgt_event[3], hgt_event[4])
                reverse_delete_seq = self.get_reverse_complement_seq(delete_seq)

                left_ins = self.extract_ref_seq(hgt_event[0], hgt_event[1] - self.window, hgt_event[1])
                left_del = delete_seq[:self.window]
                left_junc = left_ins + left_del
                left_junc_rev = left_ins + reverse_delete_seq[:self.window]
            
                mid_locus = hgt_event[3] + self.window * 3
                mid_del = self.extract_ref_seq(hgt_event[2], mid_locus-self.window, mid_locus+self.window)
                mid_del_rev = self.get_reverse_complement_seq(mid_del)


                right_del = delete_seq[-self.window:]
                right_ins = self.extract_ref_seq(hgt_event[0], hgt_event[1], hgt_event[1] + self.window)
                right_junc = right_del + right_ins
                right_junc_rev = reverse_delete_seq[-self.window:] + right_ins

                left_flag = self.verify(reads_set, left_junc, left_junc_rev)
                mid_flag = self.verify(reads_set, mid_del, mid_del_rev)
                right_flag = self.verify(reads_set, right_junc, right_junc_rev)

                if left_flag and mid_flag and right_flag:
                    verify_flag = True
                else:
                    print (left_flag,mid_flag,right_flag )

            total_hgt_event_num += 1
            final_total += 1


            # self.genetate_fasta(tmp_ref, merged_seq)
            # self.genetate_fasta(reverse_tmp_ref, rev_merged_seq)
            del_ref_len = len(self.ref_fasta[hgt_event[2]])
            ins_ref_len = len(self.ref_fasta[hgt_event[0]])
            if verify_flag:
                print (hgt_event, "HGT event is verified.", delete_len, "genome length", del_ref_len, ins_ref_len)  
                final_data.append([sample] + hgt_event + ["Yes"])
                valid_hgt_event_num += 1
                final_verified += 1
            else:
                final_data.append([sample] + hgt_event + ["No"])
                print (hgt_event, "HGT event is not verified.", delete_len, "genome length", del_ref_len, ins_ref_len)
            print ("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
            # break
        print ("Total HGT num is %s; valid one is %s; valid ratio is %s."%(total_hgt_event_num, valid_hgt_event_num, valid_hgt_event_num/total_hgt_event_num))

    def verify(self, reads_set, merged_seq, rev_merged_seq):
        uniq_dict = {}
        max_align_length = 0
        for read in reads_set:
            if read.query_sequence == None:
                continue
            if len(read.query_sequence) < self.window:
                continue
            if read.query_name in uniq_dict:
                continue
            if countN(merged_seq)/len(merged_seq) > 0.4:
                continue

            align_len = self.align(merged_seq, rev_merged_seq, read.query_sequence)
            if align_len > max_align_length:
                max_align_length = align_len
            # print (start_end_positions, len(delete_seq), len(merged_seq))
            # print ("fr", alignment[0])
            # print ("to", alignment[1])
            uniq_dict[read.query_name] = 1

            if len(merged_seq) - max_align_length < self.window/2:
                break
        print ("aligned len", len(merged_seq), max_align_length)
        if len(merged_seq) - max_align_length < self.window/2 and max_align_length > self.window: 
            return True
        else:
            return False  

    def align(self, merged_seq, rev_merged_seq, read_seq):
        from_seq = DNA(merged_seq)
        to_seq = DNA(read_seq)
        # print (from_seq)
        # print (to_seq)
        alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq)
        align_length = start_end_positions[0][1] - start_end_positions[0][0]
        # print ("for", start_end_positions, len(merged_seq))

        from_seq = DNA(rev_merged_seq)
        alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq)
        rev_align_length = start_end_positions[0][1] - start_end_positions[0][0]
        # print ("rev", start_end_positions, len(merged_seq))

        return max([align_length, rev_align_length])


    def genetate_fasta(self, file, seq):
        f = open(file, 'w')
        print (">seq\n%s"%(seq), file = f)
        f.close()

def countN(sequence):
    # initialize a counter variable
    count = 0

    # loop through the sequence and count the number of "N" characters
    for char in sequence:
        if char.upper() == 'N':
            count += 1
    return count

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
    tmp_ref = workdir + "/tmp.fasta"
    reverse_tmp_ref = workdir + "/tmp.rev.fasta"

    hgt_event_dict = read_event()
    ngs_tgs_pair = read_meta()
    final_data = []



    for sample in hgt_event_dict:
    # sample = "SRR18490939"
        bamfile = tgs_bam_dir + "/%s.bam"%(ngs_tgs_pair[sample])
        baifile = tgs_bam_dir + "/%s.bam.bai"%(ngs_tgs_pair[sample])

        if os.path.isfile(bamfile) and os.path.isfile(baifile):
            print (sample)
            ma = Map()
            ma.for_each_sample(sample)
    print ("Total HGT num is %s; valid one is %s; valid ratio is %s."%(final_total, final_verified, final_verified/final_total))

    df = pd.DataFrame(final_data, columns = ["sample", "receptor", "insert_locus", "donor", "delete_start", "delete_end", "Verified"])
    df.to_csv(verified_result, sep=',')

