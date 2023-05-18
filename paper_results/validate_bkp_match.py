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

def get_support_reads(bam, chr, pos):
    # Open the input BAM file
    bamfile = pysam.AlignmentFile(bam, "rb")
    read_list = []
    # # Create a new BAM file for the supporting reads
    # out_bamfile = pysam.AlignmentFile("supporting_reads.bam", "wb", template=bamfile)
    near = 20000
    start = pos - near
    end = pos + near
    if start < 0:
        start = 1
    
    pos_read_num = 0

    # Iterate over the reads in the BAM file
    for read in bamfile.fetch(chr, start, end):

        if read.is_unmapped:
            continue

        # Check if the read overlaps the breakpoint position
        if read.reference_start <= pos and read.reference_end >= pos:
            pos_read_num += 1
            # print(read.query_name, read.reference_name, read.reference_start, read.reference_end, read.cigar, sep="\t")
            # Write the read to the output BAM file
            # out_bamfile.write(read)
        read_list.append(read)

    # Close the input and output BAM files
    bamfile.close()
    # out_bamfile.close()
    return read_list, pos_read_num

def extract_read_seq(bam):
    # Open the BAM file
    bamfile = pysam.AlignmentFile(bam, 'rb')

    # Create an empty dictionary to store the read sequences
    read_seqs = {}

    # Iterate over each read in the BAM file
    for read in bamfile:
        # Check if the read is mapped and primary
        if not read.is_unmapped and not read.is_secondary:
            # Get the read ID and sequence
            read_id = read.query_name
            read_seq = read.query_sequence
            # Add the sequence to the dictionary using the read ID as the key
            if read_id not in read_seqs:
                read_seqs[read_id] = read_seq
            else:
                if len(read_seq) > len(read_seqs[read_id]):
                    read_seqs[read_id]  = read_seq

    # Close the BAM file
    bamfile.close()
    return read_seqs

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
        self.window = window
        self.max_length = 8000

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
        insert_support_reads, pos_read_num_1 = get_support_reads(bamfile, hgt_event[0], hgt_event[1])
        delete_support_reads_1, pos_read_num_2 = get_support_reads(bamfile, hgt_event[2], hgt_event[3])
        delete_support_reads_2, pos_read_num_3 = get_support_reads(bamfile, hgt_event[2], hgt_event[4])
        reads_set = insert_support_reads + delete_support_reads_1 + delete_support_reads_2
        pos_read_num = pos_read_num_1 + pos_read_num_2 + pos_read_num_3


        f = open(tmp_nano_str, 'w')
        uniq_dict = {}
        for read in reads_set:
            if read.query_name in uniq_dict:
                continue
            print (f">{read.query_name}\n{read_seqs[read.query_name]}", file = f)
            uniq_dict[read.query_name] = 1
        # os.system(f"shasta --memoryBacking 2M --memoryMode filesystem --assemblyDirectory {contig_dir} --input {tmp_nano_str} --config Nanopore-May2022 2> {standard_output}")
        os.system(f"flye --threads 5 --nano-raw {tmp_nano_str} --out-dir {contig_dir} 2> {standard_output}")

        return reads_set, pos_read_num

    def for_each_sample_bk(self, sample, final_total, final_verified):
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
        return final_total, final_verified

    def for_each_sample(self, sample, final_total, final_verified, best_verified):
        if sample not in hgt_event_dict:
            print ("No hgt for sample", sample)
        else:
            print ("Detect hgt for sample", sample)

        total_hgt_event_num = 0
        valid_hgt_event_num = 0
        skip_hgt_event_num = 0

        for hgt_event in sorted(list(hgt_event_dict[sample])):
            # if hgt_event[0] != "GUT_GENOME001602_16":
            #     continue
            # if 'GUT_GENOME156655_37' not in hgt_event:
            #     continue
            verify_flag = False
            reads_set, pos_read_num = self.get_reads(hgt_event)

            if pos_read_num < 5:  # skipt the event with less than five support reads
                skip_hgt_event_num += 1
                continue

            delete_len = hgt_event[4] - hgt_event[3]
            
            insert_left = self.extract_ref_seq(hgt_event[0], hgt_event[1] - self.window, hgt_event[1])
            insert_right = self.extract_ref_seq(hgt_event[0], hgt_event[1], hgt_event[1] + self.window)
            delete_seq = self.extract_ref_seq(hgt_event[2], hgt_event[3], hgt_event[4])
            reverse_delete_seq = self.get_reverse_complement_seq(delete_seq)
            merged_seq = insert_left + delete_seq + insert_right
            rev_merged_seq = insert_left + reverse_delete_seq + insert_right
            reverse_flag = hgt_event[5]
            verify_flag, best_flag = self.verify(reads_set, merged_seq, rev_merged_seq, reverse_flag)

            total_hgt_event_num += 1
            final_total += 1

            del_ref_len = len(self.ref_fasta[hgt_event[2]])
            ins_ref_len = len(self.ref_fasta[hgt_event[0]])
            if verify_flag:
                print (hgt_event, "HGT event is verified. delete length:", delete_len, "long read num", len(reads_set), "if totally coverred:", best_flag)  
                final_data.append([sample] + hgt_event + ["Yes"])
                valid_hgt_event_num += 1
                final_verified += 1
            else:
                final_data.append([sample] + hgt_event + ["No"])
                print (hgt_event, "HGT event is not verified. delete length:", delete_len, "long read num", len(reads_set), "if totally coverred:", best_flag)
            if best_flag:
                best_verified += 1
            print ("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
            # break
        print ("#Sample %s, Total HGT num is %s; valid one is %s; valid ratio is %s; skip num is %s."%(sample, total_hgt_event_num, valid_hgt_event_num, valid_hgt_event_num/total_hgt_event_num, skip_hgt_event_num))
        return final_total, final_verified, best_verified

    def verify_bk(self, reads_set, merged_seq, rev_merged_seq):
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

            if len(merged_seq) - max_align_length < self.window/2 and max_align_length > self.window:
                break
        print ("aligned len", len(merged_seq), max_align_length)
        if len(merged_seq) - max_align_length < self.window/2 and max_align_length > self.window: 
            return True
        else:
            return False  

    def align(self, merged_seq, rev_merged_seq, read_seq):
        from_seq = DNA(merged_seq)
        to_seq = DNA(read_seq)
        print ("1", merged_seq)
        print ("2", read_seq)
        # print (from_seq)
        # print (to_seq)
        alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq, gap_open_penalty=3, gap_extend_penalty=2,mismatch_score=-3 )
        align_length = start_end_positions[0][1] - start_end_positions[0][0]
        # print ("for", start_end_positions, len(merged_seq))
        # print ("for del",alignment)

        from_seq = DNA(rev_merged_seq)
        alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq, gap_open_penalty=3, gap_extend_penalty=2,mismatch_score=-3)
        rev_align_length = start_end_positions[0][1] - start_end_positions[0][0]
        # print ("rev", start_end_positions, len(merged_seq))
        # print ("rev del",alignment)
    
        from_seq = DNA(merged_seq)
        to_seq = DNA(self.get_reverse_complement_seq(read_seq))
        alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq, gap_open_penalty=3, gap_extend_penalty=2,mismatch_score=-3)
        rev_read_align_length = start_end_positions[0][1] - start_end_positions[0][0]
        # print ("for", start_end_positions, len(merged_seq))
        # print ("for del, for read",alignment)

        from_seq = DNA(rev_merged_seq)
        alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq, gap_open_penalty=3, gap_extend_penalty=2,mismatch_score=-3)
        rev_read_rev_align_length = start_end_positions[0][1] - start_end_positions[0][0]
        # print ("rev", start_end_positions, len(merged_seq))
        # print ("rev del, rev read",alignment)

        # return max([align_length, rev_align_length])
        return max([align_length, rev_align_length, rev_read_align_length, rev_read_rev_align_length])

    def verify_bk2(self, reads_set, merged_seq, rev_merged_seq, reverse_flag):
        uniq_dict = {}
        max_align_length = 0
        verify = False
        best_flag = False
        all_start_flag = False
        all_end_flag = False

        start_interval = [float("inf"), 0] # the mapped region of the reads mapped to start junctions
        end_interval = [float("inf"), 0] # the mapped region of the reads mapped to end junctions

        for read in reads_set:
            read.query_sequence = read_seqs[read.query_name]
            # print (read.query_name, read.query_sequence)
            if read.query_sequence == None:
                continue
            if len(read.query_sequence) < self.window:
                continue
            if read.query_name in uniq_dict:
                continue
            if countN(merged_seq)/len(merged_seq) > 0.4:
                continue
            
            if reverse_flag == "True":
                start_flag, end_flag, best_flag, start_match_interval, end_match_interval = minimap2_align(rev_merged_seq, read.query_sequence)
            else:
                start_flag, end_flag, best_flag, start_match_interval, end_match_interval = minimap2_align(merged_seq, read.query_sequence)


            all_start_flag = all_start_flag or  start_flag
            all_end_flag = all_end_flag or end_flag

            if start_flag and end_flag:
                verify = True

            if start_flag:
                if start_match_interval < start_interval[0]:
                    start_interval[0] = start_match_interval
                if end_match_interval > start_interval[1]:
                    start_interval[1] = end_match_interval
            if end_flag:
                if start_match_interval < end_interval[0]:
                    end_interval[0] = start_match_interval
                if end_match_interval > end_interval[1]:
                    end_interval[1] = end_match_interval

            if len(merged_seq) > self.max_length:
                if all_start_flag and all_end_flag:
                    if end_interval[0] < start_interval[1]:
                        verify = True
            # print ("<<<<<<<<<<<<<", start_interval, end_interval, start_flag, end_flag)
            # else:

            if best_flag:
                break

            uniq_dict[read.query_name] = 1

        # if len(merged_seq) < self.max_length: # check the whole inserted segment if short
        #     verify = best_flag

        return verify, best_flag 

    def verify(self, reads_set, merged_seq, rev_merged_seq, reverse_flag):
        uniq_dict = {}
        max_align_length = 0
        verify = False
        best_flag = False
        all_start_flag = False
        all_end_flag = False
        
        

        start_interval = [float("inf"), 0] # the mapped region of the reads mapped to start junctions
        end_interval = [float("inf"), 0] # the mapped region of the reads mapped to end junctions

        sequence_list = []
        if os.path.isfile(contig_file):
            contigs = Fasta(contig_file)
            for name in contigs.keys():
                sequence = str(contigs[name])
                sequence_list.append([name, sequence])
        for read in reads_set:
            if read.query_name in uniq_dict:
                continue
            sequence = read_seqs[read.query_name]
            sequence_list.append([read.query_name, sequence])
            uniq_dict[read.query_name] = 1

        for name, sequence in sequence_list:
            if sequence == None:
                continue
            if len(sequence) < self.window:
                continue
            # if read.query_name in uniq_dict:
            #     continue
            if countN(merged_seq)/len(merged_seq) > 0.4:
                continue
            
            if reverse_flag == "True":
                start_flag, end_flag, best_flag, start_match_interval, end_match_interval = minimap2_align(rev_merged_seq, sequence, name)
            else:
                start_flag, end_flag, best_flag, start_match_interval, end_match_interval = minimap2_align(merged_seq, sequence, name)


            all_start_flag = all_start_flag or  start_flag
            all_end_flag = all_end_flag or end_flag

            if start_flag and end_flag:
                verify = True

            if start_flag:
                if start_match_interval < start_interval[0]:
                    start_interval[0] = start_match_interval
                if end_match_interval > start_interval[1]:
                    start_interval[1] = end_match_interval
            if end_flag:
                if start_match_interval < end_interval[0]:
                    end_interval[0] = start_match_interval
                if end_match_interval > end_interval[1]:
                    end_interval[1] = end_match_interval

            if len(merged_seq) > self.max_length:
                if all_start_flag and all_end_flag:
                    if end_interval[0] < start_interval[1]:
                        verify = True
            # print ("<<<<<<<<<<<<<", start_interval, end_interval, start_flag, end_flag)
            # else:

            if best_flag:
                break
            # uniq_dict[read.query_name] = 1
        # if len(merged_seq) < self.max_length: # check the whole inserted segment if short
        #     verify = best_flag
        os.system(f"rm -rf {contig_dir}")
        return verify, best_flag 

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

def write_fasta_file(filename, seq_name, sequence):
    """Write a sequence to a FASTA file"""
    with open(filename, "w") as f:
        f.write(f">{seq_name}\n{sequence}\n")

def minimap2_align(seq1, seq2, name):
    # Example sequences
    seq1_file = workdir + "/seq1.fasta"
    seq2_file = workdir + "/seq2.fasta"
    paf = workdir + "/minimap2.paf"

    # Write sequences to FASTA files
    write_fasta_file(seq1_file, "seq1", seq1)
    write_fasta_file(seq2_file, "seq2", seq2)
    # sam = "/mnt/d/HGT/time_lines/minimap2.sam"
    

    # cmd = f"minimap2 {seq1_file} {seq2_file} -a -o {sam}"
    cmd = f"minimap2 -t 10 -x map-ont {seq2_file} {seq1_file} -o {paf} 2> {standard_output}"
    os.system(cmd)
    # os.system(f"cat {paf}|cut -f 1-4")
    return parse_paf(paf, name)

def parse_paf(paf_file, name):
    """Parse a PAF file and extract the Target start and end positions"""
    map_flag = False
    start_flag = False
    end_flag = False
    best_flag = False
    start_match_interval = float('inf') # the mapped interval of the reads 
    end_match_interval = 0 # the mapped interval of the reads 
    junc_len = 30
    with open(paf_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            target_start = int(fields[2])
            target_end = int(fields[3])
            target_len = int(fields[1])
            read_len = int(fields[6])

            # start_junc = [window/2, 1.5*window] 
            # end_junc = [target_len - 1.5*window, target_len-window/2]

            start_junc = [window-junc_len, window+junc_len] 
            end_junc = [target_len - window - junc_len, target_len-window + junc_len]
            
            # start_junc = [window/2, target_len/2] 
            # if start_junc[1] - start_junc[0] > max_seg:
            #     start_junc[1] = start_junc[0] + max_seg
            # end_junc = [target_len/2, target_len-window/2]
            # if end_junc[1] - end_junc[0] > max_seg:
            #     end_junc[0] = end_junc[1] - max_seg

            if start_junc[0] >= target_start and start_junc[1] <= target_end:
                start_flag = True
            if end_junc[0] >= target_start and end_junc[1] <= target_end:
                end_flag = True
            if start_flag and end_flag:
                map_flag = True
            if target_start < start_match_interval:
                start_match_interval = target_start
            if target_end > end_match_interval:
                end_match_interval = target_end
            # if map_flag:
            #     break
            # print (name, read_len, target_start, target_end, "|", fields[7], fields[8], start_junc, end_junc, sep = "\t")
            if start_junc[0] >= target_start and end_junc[1] <= target_end:
                best_flag = True
            if best_flag:
                break
    # print (start_flag, end_flag)
    # print (name, "<<<<<<<<<<<<<<<<<<<<<<")
    return start_flag, end_flag, best_flag, start_match_interval, end_match_interval

if __name__ == "__main__":

    database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"
    workdir = "/mnt/d/HGT/time_lines/"
    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    data_pair = "/mnt/d/HGT/time_lines/SRP366030.ngs_tgs_pair.csv"
    design_file = "/mnt/d/HGT/time_lines/sample_design.tsv"
    result_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    tgs_bam_dir = "/mnt/d/HGT/time_lines/tgs_bam_results"
    verified_result = "/mnt/d/HGT/time_lines/SRP366030.verified_event.csv"

    # database = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/reference/UHGG_reference.formate.fna"
    # workdir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/"
    # meta_data = "//mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/SRP366030.csv.txt"
    # data_pair = "/mnt/disk2_workspace/wangshuai/00.strain/32.BFB/SRP366030.ngs_tgs_pair.csv"
    # design_file = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid//sample_design.tsv"
    # result_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/hgt/result/"
    # identified_hgt = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/match/SRP366030.identified_event.csv"
    # tgs_bam_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/nanopore_alignment/results/"
    # verified_result = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/match/SRP366030.verified_event.csv"


    tmp_bam = workdir + "/tmp.bam"
    tmp_fastq = workdir + "/tmp.fastq"
    tmp_nano_str = workdir + "/tmp.fasta"
    reverse_tmp_ref = workdir + "/tmp.rev.fasta"
    contig_dir = workdir + "/flye/"
    contig_file = contig_dir + "/assembly.fasta"
    standard_output = workdir + "/minimap2.log"

    hgt_event_dict = read_event()
    ngs_tgs_pair = read_meta()
    final_data = []

    final_total = 0
    final_verified = 0
    best_verified = 0
    window = 200
    max_seg = 5000

    # minimap2_align()

    for sample in hgt_event_dict:
        # if sample != "SRR18490939":
        #     continue
        bamfile = tgs_bam_dir + "/%s.bam"%(ngs_tgs_pair[sample])
        baifile = tgs_bam_dir + "/%s.bam.bai"%(ngs_tgs_pair[sample])

        if os.path.isfile(bamfile) and os.path.isfile(baifile):
            print (sample)
            read_seqs = extract_read_seq(bamfile)
            ma = Map()
            final_total, final_verified, best_verified = ma.for_each_sample(sample, final_total, final_verified, best_verified)
    print ("Total HGT num is %s; valid one is %s; valid ratio is %s, best valid ratio is %s."%(final_total, final_verified, final_verified/final_total, best_verified/final_total))

    df = pd.DataFrame(final_data, columns = ["sample", "receptor", "insert_locus", "donor", "delete_start", "delete_end", "reverse_flag", "Verified"])
    df.to_csv(verified_result, sep=',')

