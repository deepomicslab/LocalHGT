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

        self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        self.abundance = None
        self.split_abundance = None

        self.read = None

def read_bkp(bkp_file):
    my_bkps = []
    f = open(bkp_file)
    all_rows = csv.reader(f)
    total_HGT_split_num = 0
    for row in all_rows:
        if row[0][0] == "#":
            reads_num = int(row[0].split(";")[0].split(":")[1])
            # print (row, self.reads_num)
            pass
        elif row[0] == "from_ref":
            pass
        else:
            if reads_num == 0:
                print ("old bkp", bkp_file)
                break
            eb = Acc_Bkp(row)
            eb.abundance = eb.cross_split_reads/reads_num
            if eb.from_ref_genome == eb.to_ref_genome:
                continue
            if eb.cross_split_reads < split_cutoff:
                continue
            if eb.abundance < abun_cutoff:
                continue
            total_HGT_split_num += eb.cross_split_reads
            # if eb.hgt_tag not in self.all_hgt:
            my_bkps.append(eb)
    for eb in my_bkps:
        eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
    return my_bkps

def find_chr_segment_name(bed_file):
    # bed_file =  "/mnt/d/breakpoints/HGT/test_5_11/species20_snp0.01_depth30_reads150_sample_1.interval.txt.bed"  
    chr_segments = {} 
    for line in open(bed_file):   
        segment = line.strip()
        my_chr = segment.split(":")[0]
        start = int(segment.split(":")[1].split("-")[0])
        end = int(segment.split(":")[1].split("-")[1])
        if my_chr not in chr_segments:
            chr_segments[my_chr] = []
        chr_segments[my_chr].append([start, end])
    return chr_segments

def convert_chr2_segment(ref, pos):
    tolerate_gap = 150
    for interval in chr_segments[ref]:
        if pos >= interval[0] - tolerate_gap and pos <= interval[1] + tolerate_gap:
            new_pos = pos - interval[0] 
            segment_name = "%s:%s-%s"%(ref, interval[0], interval[1])
            if new_pos < 1:
                new_pos = 1
            return segment_name, new_pos
    print ("Can't find corresponding for", ref, pos)
    return "NA", 0

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
            # Write the read to the output BAM file
            # out_bamfile.write(read)
            read_list.append(read)

    # Close the input and output BAM files
    bamfile.close()
    # out_bamfile.close()
    return read_list

def get_tandem_repeat():
    tandem_repeat_dict = {}
    f = open(tandem_repeat)
    for line in f:   
        if line[0] == "#":
            continue
        array = line.split()
        if array[2] == "region":
            continue
        chrom = array[0]
        repeat_start = int(array[3])
        repeat_end = int(array[4])
        if chrom not in tandem_repeat_dict:
            tandem_repeat_dict[chrom] = []
        tandem_repeat_dict[chrom].append([repeat_start, repeat_end])
    f.close()
    # print (len(tandem_repeat_dict))
    return tandem_repeat_dict

def get_TEI():
    TEI_dict = {}
    f = open(TEI)
    for line in f:   
        if line[0] == "#":
            continue
        array = line.split()
        if array[2] == "region":
            continue
        chrom = array[0]
        repeat_start = int(array[1])
        repeat_end = int(array[2])
        if chrom not in TEI_dict:
            TEI_dict[chrom] = []
        TEI_dict[chrom].append([repeat_start, repeat_end])
    f.close()
    return TEI_dict

def is_in_intervals(value, intervals):
    for interval in intervals:
        if value >= interval[0] and value <= interval[1]:
            return True
    return False

class Mechanism():

    def __init__(self, bkp_object, ref_fasta):
        self.bkp = bkp_object
        # self.ref = database_dir + "/UHGG_reference.formate.fna"
        self.ref_fasta = ref_fasta 
        self.cutoff = 100
        

    def extract_ref_seq(self, scaffold_name, start, end):
        if start < 1:
            start = 1
        # print (scaffold_name, start, end)
        return self.ref_fasta[scaffold_name][start:end].seq

    def compare_seq_ins(self, chrom, pos, strand):
        # chrom = self.bkp.from_ref
        # pos = self.bkp.from_bkp
        segment_name, segment_pos = convert_chr2_segment(chrom, pos)
        read_list = get_support_reads(bam, segment_name, segment_pos)
        for read in read_list:
            read_seq = read.query_sequence
            if len(read_seq) < 50:
                continue
            ref_seq = self.extract_ref_seq(chrom, pos - len(read_seq), pos)
            if re.search("N", ref_seq) or re.search("N", read_seq):
                continue
            if strand == "-":
                ref_seq = self.get_reverse_complement_seq(ref_seq)   
            read_seq = DNA(read_seq)
            ref_seq = DNA(ref_seq)

            alignment, score, start_end_positions = global_pairwise_align_nucleotide(read_seq, ref_seq)
            # print (read.query_name, alignment[0])
            # print (read.query_name, alignment[1])
            insertion_list = extract_insertion(str(alignment[1]))
            if len(insertion_list) == 0:
                return 0
            else:
                return max(insertion_list)
        return 0

    def compare_seq_homo(self):
        self.bkp.from_bkp -= 1
        from_seq = self.extract_ref_seq(self.bkp.from_ref, self.bkp.from_bkp-self.cutoff, self.bkp.from_bkp+self.cutoff)
        to_seq = self.extract_ref_seq(self.bkp.to_ref, self.bkp.to_bkp-self.cutoff, self.bkp.to_bkp+self.cutoff)
        if self.bkp.from_strand == "-":
            from_seq = self.get_reverse_complement_seq(from_seq)      
        if self.bkp.to_strand == "-":
            to_seq = self.get_reverse_complement_seq(to_seq)
        if re.search(">", from_seq) or re.search(">", to_seq) or re.search("N", from_seq) or re.search("N", to_seq):
            return 0
        from_seq = DNA(from_seq)
        to_seq = DNA(to_seq)
        alignment, score, start_end_positions = global_pairwise_align_nucleotide(from_seq, to_seq)
        # print ("fr", alignment[0])
        # print ("to", alignment[1])
        max_homology_len = extract_homology(alignment)
        return max_homology_len

    def get_reverse_complement_seq(self, sequence):
        sequence = sequence[::-1]
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
        string = sequence.translate(trantab)
        return string

    def check_tandem(self, scaffold_name, pos):
        # check if the breakpoint locates in the tandem repeat region
        if scaffold_name not in tandem_repeat_dict:
            return False
        else:
            return is_in_intervals(pos, tandem_repeat_dict[scaffold_name])

    def check_TEI(self, scaffold_name, pos):
        # check if the breakpoint locates in the tandem repeat region
        if scaffold_name not in TEI_dict:
            return False
        else:
            return is_in_intervals(pos, TEI_dict[scaffold_name])

    def main(self, from_type, to_type): 
        from_tandem_repeat_flag = self.check_tandem(self.bkp.from_ref, self.bkp.from_bkp)
        from_TEI_flag = self.check_TEI(self.bkp.from_ref, self.bkp.from_bkp)
        to_tandem_repeat_flag = self.check_tandem(self.bkp.to_ref, self.bkp.to_bkp)
        to_TEI_flag = self.check_TEI(self.bkp.to_ref, self.bkp.to_bkp)

        from_ins =  self.compare_seq_ins(self.bkp.from_ref, self.bkp.from_bkp, self.bkp.from_strand)
        to_ins =  self.compare_seq_ins(self.bkp.to_ref, self.bkp.to_bkp, self.bkp.to_strand)

        homo = self.compare_seq_homo()
        # print (tandem_repeat_flag, from_ins, to_ins, homo)
        from_mechanism = self.classify(from_type, from_tandem_repeat_flag, from_TEI_flag, from_ins, homo)
        to_mechanism = self.classify(to_type, to_tandem_repeat_flag, to_TEI_flag, to_ins, homo)

        if from_type == "del":
            return from_mechanism, self.bkp.from_ref, self.bkp.from_bkp
        elif to_type == "del":
            return to_mechanism, self.bkp.to_ref, self.bkp.to_bkp
        # print (from_mechanism, to_mechanism)
    
    def classify(self, break_type, tandem_repeat_flag, TEI_flag, ins_num, homo_num):
        if break_type == "ins":
            if TEI_flag:
                mechanism_type = "TEI"
            elif tandem_repeat_flag:
                mechanism_type = "VNTR"
            else:
                mechanism_type = "NA"
        elif break_type == "del":
            if TEI_flag:
                mechanism_type = "TEI"
            elif tandem_repeat_flag:
                mechanism_type = "VNTR"
            else:
                if ins_num > 0:
                    if ins_num > 10:
                        mechanism_type = "FoSTeS/MMBIR"
                    else:
                        mechanism_type = "NHEJ"
                else:
                    if homo_num > 100:
                        mechanism_type = "NAHR"
                    elif homo_num >= 2:
                         mechanism_type = "alt-EJ"
                    else:
                        mechanism_type = "NHEJ"
        return mechanism_type

def read_verified(verified_result):
    verified_HGT_dict = {}
    record_HGT = {}
    for line in open(verified_result):
        array = line.strip().split(",")
        if array[1] == "sample":
            continue
        if array[-1] != "Yes":
            continue
        sample = array[1]

        if sample not in verified_HGT_dict:
            verified_HGT_dict[sample] = {}
            record_HGT[sample] = {}
        insert = array[2] + "&" + array[3]
        delete1 = array[4] + "&" + array[5]
        delete2 = array[4] + "&" + array[6]
        verified_HGT_dict[sample][insert] = "ins"
        verified_HGT_dict[sample][delete1] = "del"
        verified_HGT_dict[sample][delete2] = "del"

        record_HGT[sample][array[2] + "&" + array[3] + "&" + array[4] + "&" + array[5]] = 1
        record_HGT[sample][array[4] + "&" + array[5] + "&" + array[2] + "&" + array[3]] = 1

        record_HGT[sample][array[2] + "&" + array[3] + "&" + array[4] + "&" + array[6]] = 1
        record_HGT[sample][array[4] + "&" + array[6] + "&" + array[2] + "&" + array[3]] = 1

    return verified_HGT_dict, record_HGT


def extract_insertion(ref_alignment):
    insertion_list = []
    insertion_length = 0
    insertion_flag = False
    for i in range(1, len(ref_alignment)-1):
        if ref_alignment[i -1] != "-" and ref_alignment[i] == "-":
            insertion_length = 1
            insertion_flag = True
        elif insertion_flag and ref_alignment[i] == "-":
            insertion_length += 1
        elif insertion_flag and ref_alignment[i] != "-":
            insertion_list.append(insertion_length)    
            insertion_flag = False
            insertion_length = 0
    # print (insertion_list)
    return insertion_list

def extract_homology(alignment):
    homology_list = []
    
    homology_flag = False
    homology_seq = ''
    for i in range(len(alignment[0])):
        if str(alignment[0][i]) != "-" and str(alignment[1][i]) != "-":
            # print (i, alignment[0][i], alignment[1][i])
            if homology_flag == False:
                homology_flag = True
                homology_seq += str(alignment[0][i])
            else:
                homology_seq += str(alignment[0][i])
        else:
            if len(homology_seq) > 0:
                homology_list.append(homology_seq)
            homology_flag = False
            homology_seq = ''
    if len(homology_seq) > 0:
        homology_list.append(homology_seq)  
    homology_len_list = [len(x) for x in homology_list]
    # if len(homology_len_list) > 0 and max(homology_len_list) > 20:
    #     print (homology_list)
    #     print (alignment)
    # print (homology_list)
    # return homology_list
    if len(homology_len_list) > 0:
        return max(homology_len_list)
    else:
        return 0
    
    
    
if __name__ == "__main__":
    bin_size = 500
    split_cutoff = 0  #10
    sample_cutoff = 8  # 8
    abun_cutoff = 1e-7  #1e-7

    # result_dir = "/mnt/d/HGT/time_lines/"
    # bkp_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    # database_dir = "/mnt/d/HGT/UHGG/"
    # verified_result = "/mnt/d/HGT/time_lines/SRP366030.verified_event.csv"
    

    result_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/hgt/result/"
    bkp_dir = result_dir
    database_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/reference/"
    verified_result = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/match/SRP366030.verified_event.csv"

    mechanism_result = result_dir + "/mechanism_result.txt"
    fout = open(mechanism_result, 'w')

    tandem_repeat = database_dir + "/UHGG_reference.formate.tandem_repeat.gff"  # Source-version MISA 2.1
    TEI = database_dir + "/TEI_elements.txt" # repeatmasker: LTR, LINE, SINE
    # sample = "SRR18491235"
    tandem_repeat_dict = get_tandem_repeat()
    TEI_dict = get_TEI()
    verified_HGT_dict, record_HGT = read_verified(verified_result)
    mechanism_freq_dict = {}


    ref_fasta = Fasta(database_dir + "/UHGG_reference.formate.fna")
    for sample in verified_HGT_dict:
        bam = "%s/%s.unique.bam"%(result_dir, sample)
        bed = "%s/%s.interval.txt.bed"%(result_dir, sample)
        bkp_file = bkp_dir + "/%s.acc.csv"%(sample)
        if not (os.path.isfile(bam) and os.path.isfile(bed) and os.path.isfile(bkp_file)):
            continue
        chr_segments = find_chr_segment_name(bed)
        my_bkps = read_bkp(bkp_file)
        print (sample, "N.O. of bkp:", len(my_bkps))
        for bkp in my_bkps:          
            # if bkp.from_ref + "&" + str(bkp.from_bkp) in verified_HGT_dict[sample] and bkp.to_ref + "&" + str(bkp.to_bkp) in verified_HGT_dict[sample]:
            if bkp.from_ref + "&" + str(bkp.from_bkp) + "&" + bkp.to_ref + "&" + str(bkp.to_bkp) in record_HGT[sample]:
                # print (bkp.from_ref, bkp.from_bkp)
                mac = Mechanism(bkp, ref_fasta)
                from_type = verified_HGT_dict[sample][bkp.from_ref + "&" + str(bkp.from_bkp)]
                to_type = verified_HGT_dict[sample][bkp.to_ref + "&" + str(bkp.to_bkp)]
                mechanism, bkp_chrom, bkp_pos = mac.main(from_type, to_type)
                print (sample, bkp_chrom, bkp_pos, mechanism)
                print (sample, bkp_chrom, bkp_pos, mechanism, file = fout)
                if mechanism not in mechanism_freq_dict:
                    mechanism_freq_dict[mechanism] = 0
                mechanism_freq_dict[mechanism] += 1
                # break

    total = sum(list(mechanism_freq_dict.values()))
    for mechanism in mechanism_freq_dict:
        print ("#", mechanism, mechanism_freq_dict[mechanism], mechanism_freq_dict[mechanism]/total, file = fout)
        print ("#", mechanism, mechanism_freq_dict[mechanism], mechanism_freq_dict[mechanism]/total)
    
    fout.close()