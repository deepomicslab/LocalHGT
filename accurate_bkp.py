#!/usr/bin/env python3

import csv
import pysam
from pyfaidx import Fasta
import skbio
from skbio import DNA, TabularMSA
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import StripedSmithWaterman
import numpy as np
import argparse
import sys
import re


from get_raw_bkp import getInsertSize, readFilter

min_match_score = 0.8 #1.6  #0.8
min_seq_len = 15
cigar_dict = {0:'M',1:'M',2:'M',3:'M',4:'N',5:'N'}
tolerate_read_mismatch_num = 20
bkp2end = 15

def compute_scores(dna1, dna2):
    # StripedSmithWaterman docs:
    # http://scikit-bio.org/docs/0.4.2/generated/skbio.alignment.StripedSmithWaterman.html
    ssw1 = StripedSmithWaterman(dna1)  #match_score=1, mismatch_score=0
    # AlignmentStructure docs:
    # http://scikit-bio.org/docs/0.4.2/generated/skbio.alignment.AlignmentStructure.html
    # https://github.com/biocore/scikit-bio/blob/9dc60b4248912a4804c90d0132888d6979a62d51/skbio/alignment/_lib/ssw.c
    align = ssw1(dna2)  #the map score is equal to the match base number
    return align.optimal_alignment_score

def get_reverse_complement_seq(sequence):
    sequence = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

class Each_Bkp(object):

    def __init__(self, one_bkp):
        self.ref1 = one_bkp[0].strip()       
        self.ref1_positions = [int(one_bkp[1]), int(one_bkp[2]), int(one_bkp[3])] 
        self.ref2 = one_bkp[4].strip()       
        self.ref2_positions = [int(one_bkp[5]), int(one_bkp[6]), int(one_bkp[7])] 
        self.split_reads_num = one_bkp[8]
        self.direction =  one_bkp[9].strip() 
        self.support_reads = []
        self.pos1 = 0
        self.pos2 = 0

    def add_support_reads(self, read_obj):
        self.support_reads.append(read_obj)

    def sort_support_reads(self):
        mean_pos1 = np.median(self.ref1_positions)
        mean_pos2 = np.median(self.ref2_positions)
        dist = {}
        record = {}
        # print (self.support_reads)
        for read_obj in self.support_reads:
            record[read_obj.qname] = read_obj
            dist[read_obj.qname] = abs(read_obj.pos1 - mean_pos1) + abs(read_obj.pos2 - mean_pos2)
        new_dist = sorted(dist.items(), key=lambda item: item[1])
        # print (new_dist)
        new_reads = []
        for new_d in new_dist:
            qname = new_d[0]
            new_reads.append(record[qname])
        self.support_reads = new_reads


    def reverse_bkp(self):
        a = self.ref1
        b = self.ref1_positions
        self.ref1 = self.ref2
        self.ref1_positions = self.ref2_positions
        self.ref2 = a
        self.ref2_positions = b
        return self

def key_name(a, b):
    return str(a) + '=' + str(b)

def reverse_name(name):
    diff = name.split('=')
    return str(diff[1]) + '=' + str(diff[0])

class Read_Raw_Bkp():
    def __init__(self, raw_bkp_file):
        self.raw_bkp_file = raw_bkp_file
        self.raw_bkps = []
        self.read_bkp()
        self.raw_bkps_cluster = {}
        self.max_dist = 50  #rlen

    def read_bkp(self):
        f = open(self.raw_bkp_file)
        all_rows = csv.reader(f)
        for row in all_rows:
            eb = Each_Bkp(row)
            self.raw_bkps.append(eb)
        print ('raw bkp num is', len(self.raw_bkps))

    def cluster_bkp(self):
        for bkp in self.raw_bkps:
            if key_name(bkp.ref1, bkp.ref2) in self.raw_bkps_cluster:
                self.raw_bkps_cluster[key_name(bkp.ref1, bkp.ref2)] = \
                self.update_cluster(self.raw_bkps_cluster[key_name(bkp.ref1, bkp.ref2)], bkp)
            elif key_name(bkp.ref2, bkp.ref1) in self.raw_bkps_cluster:
                bkp.reverse_bkp()
                self.raw_bkps_cluster[key_name(bkp.ref1, bkp.ref2)] = \
                self.update_cluster(self.raw_bkps_cluster[key_name(bkp.ref1, bkp.ref2)], bkp)            
            else:
                self.raw_bkps_cluster[key_name(bkp.ref1, bkp.ref2)] = [bkp]
        print ('gene pair number is', len(self.raw_bkps_cluster))

    def update_cluster(self, xy_cluster, bkp):
        flag = False
        for cluster in xy_cluster:
            if bkp.ref1 == cluster.ref1 and bkp.ref2 == cluster.ref2 and bkp.direction == cluster.direction:
                if abs(bkp.ref1_positions[0] - cluster.ref1_positions[0]) < self.max_dist and\
                    abs(bkp.ref2_positions[0] - cluster.ref2_positions[0]) < self.max_dist:
                    cluster.ref1_positions += bkp.ref1_positions
                    cluster.ref2_positions += bkp.ref2_positions
                    flag = True
                    # break
            elif bkp.ref1 == cluster.ref2 and bkp.ref2 == cluster.ref1 and bkp.direction == cluster.direction:
                if abs(bkp.ref1_positions[0] - cluster.ref2_positions[0]) < self.max_dist and\
                    abs(bkp.ref2_positions[0] - cluster.ref1_positions[0]) < self.max_dist:
                    cluster.ref1_positions += bkp.ref2_positions
                    cluster.ref2_positions += bkp.ref1_positions
                    flag = True
                    # break
        if flag == False:
            xy_cluster.append(bkp)
        return xy_cluster

    def sort_cluster(self):
        for my_key in self.raw_bkps_cluster:
            for cluster in self.raw_bkps_cluster[my_key]:
                cluster.ref1_positions = sorted(list(set(cluster.ref1_positions)))
                cluster.ref2_positions = sorted(list(set(cluster.ref2_positions)))

class Each_Split_Read(object):
    def __init__(self, read):
        self.clipped_direction = ''
        
        self.ref1 = read.reference_name
        self.ref2 = read.get_tag('SA').split(',')[0]
        self.raw_ref1 = self.ref1
        self.raw_ref2 = self.ref2
        self.ref2_cigar = read.get_tag('SA').split(',')[3]
        self.pos1 = read.reference_start
        self.pos2 = int(read.get_tag('SA').split(',')[1])
        self.qname = read.qname
        self.hard_clipped = 0
        self.end_point = False

        if self.qname not in reads_mapped_len:
            reads_mapped_len[self.qname] = 0


        self.ref2_clipped_direction = ''  #covert ref2 pos
        self.direction_confict = False
        self.get_ref2_clipped_direction()

        m = self.map_length(read)
        if self.clipped_direction == 'right':
            self.pos1 += m
        if self.clipped_direction == 'left':
            self.mapped_len = rlen - m
        else:
            self.mapped_len = m

        if args['n'] == 1:
            self.update_pos()
        self.real_ref = self.ref1
        
        self.clipped = 2 #indicate this reads is clipped in ref2, so only use it to find acc in ref2.
        if len(read.query_sequence) < rlen:
            self.seq1 = '' #seq1 save the clipped sequence for ref1
            self.seq2 = '' #seq1 save the clipped sequence for ref2
        else:
            if self.clipped_direction == 'right':
                self.seq1 = read.query_sequence[:m]
                self.seq2 = read.query_sequence[m:]
            else:
                self.seq1 = read.query_sequence[m:]  
                self.seq2 = read.query_sequence[:m]
        # if self.qname == "GUT_GENOME036853_4-1172_2":
        #     print (self.qname, self.ref1, self.pos1, self.ref2, self.pos2, self.clipped_direction, len(self.seq1), len(self.seq2))

    def get_ref2_clipped_direction(self):
        left_re = re.search('^(\d+)([SH])', self.ref2_cigar)
        right_re = re.search('(\d+)([SH])$', self.ref2_cigar)
        left, right = 0, 0
        if left_re:
            
            left = int(left_re.group(0)[:-1])
            
        if right_re:
            # print ("right", self.ref2_cigar)
            right = int(right_re.group(0)[:-1])
            
        if left > right:
            self.ref2_clipped_direction = "left"
        else:
            self.ref2_clipped_direction = "right"
            self.pos2 += (rlen-right)
            
        if self.ref2_clipped_direction == self.clipped_direction:
            self.direction_confict = True

    def update_pos(self):
        seg1_len = int(self.ref1.split(':')[1].split('-')[1]) - int(self.ref1.split(':')[1].split('-')[0])
        seg2_len = int(self.ref2.split(':')[1].split('-')[1]) - int(self.ref2.split(':')[1].split('-')[0])
        seg1_start = int(self.ref1.split(':')[1].split('-')[0])
        seg2_start = int(self.ref2.split(':')[1].split('-')[0])
        if (seg1_start > 100 and self.pos1 < bkp2end) or (seg2_start > 100 and self.pos2 < bkp2end) \
            or seg1_len - self.pos1 < bkp2end or seg2_len - self.pos2 < bkp2end:
            self.end_point = True
        self.pos1 += int(self.ref1.split(':')[1].split('-')[0])
        self.pos2 += int(self.ref2.split(':')[1].split('-')[0])

        self.ref1 = self.ref1.split(':')[0]
        self.ref2 = self.ref2.split(':')[0]

    def map_length(self, read):
        l, r = 0, 0
        cigar_len = len(read.cigar)

        for ci in read.cigar: # record left clipped length
            if cigar_dict[ci[0]] != 'M':
                l += ci[1]
            else:
                break
        for i in range(cigar_len): # record right clipped length
            ci = read.cigar[cigar_len-1-i]
            if cigar_dict[ci[0]] != 'M':
                r += ci[1]
            else:
                break

        for ci in read.cigar: # record mapped_length
            if cigar_dict[ci[0]] == 'M':
                reads_mapped_len[self.qname] += ci[1]

        if r > l :
            self.clipped_direction = 'right'
            return rlen - r
        else:
            self.clipped_direction = 'left'
            return l

    def reverse_read_obj(self):
        a = self.ref1
        b = self.pos1
        c = self.seq1
        self.ref1 = self.ref2
        self.pos1 = self.pos2
        self.seq1 = self.seq2
        self.ref2 = a
        self.pos2 = b
        self.seq2 = c
        if self.clipped == 2:
            self.clipped = 1

def read_split_bam(split_bam_name):
    used_reads_num = 0
    split_bamfile = pysam.AlignmentFile(filename = split_bam_name, mode = 'r')
    for read in split_bamfile:
        if not read.has_tag('SA'): 
            continue   
        read_obj = Each_Split_Read(read)
        if read_obj.ref1 == read_obj.ref2:
            continue
        if len(read_obj.seq1) == 0 and len(read_obj.seq2) == 0:
            continue
        xy_name = key_name(read_obj.ref1, read_obj.ref2)
        if xy_name in rrm.raw_bkps_cluster.keys():
            for cluster in rrm.raw_bkps_cluster[xy_name]: #for a R1 and R2
                add_support_split_reads(cluster, read_obj)
        elif reverse_name(xy_name) in rrm.raw_bkps_cluster.keys():
            read_obj.reverse_read_obj()
            for cluster in rrm.raw_bkps_cluster[reverse_name(xy_name)]: #for a R1 and R2
                add_support_split_reads(cluster, read_obj)
                used_reads_num += 1
        else:
            continue
    # print (used_reads_num)
        # print (read.reference_name, read.cigar, read.get_tag('SA'))

def add_support_split_reads(cluster, read_obj):
    max_dist_support_read = insert_size #rlen #
    flag = False
    for i in range(len(cluster.ref1_positions)):
        for j in range(len(cluster.ref2_positions)):
            if abs(read_obj.pos1 - cluster.ref1_positions[i]) < max_dist_support_read and \
            abs(read_obj.pos2 - cluster.ref2_positions[j]) < max_dist_support_read:   
                cluster.add_support_reads(read_obj)
                # if read_obj.qname == "GUT_GENOME036853_4-1172_2":
                #     print (read_obj.qname, cluster.ref1, cluster.ref1_positions, cluster.ref2, cluster.ref2_positions)
                flag = True
                break
        if flag:
            break
    
def extract_ref_seq(scaffold_name, start, end):
    if start < 1:
        start = 1
    return ref_fasta[scaffold_name][start:end].seq

def find_accurate_bkp():
    bkp_num_support = 0
    bkp_num = 0
    for species_pair in rrm.raw_bkps_cluster:
        raw_bkp_clusters = rrm.raw_bkps_cluster[species_pair]
        # print (len(raw_bkp_clusters))
        for cluster in raw_bkp_clusters:
            # if cluster.ref1 in ['GUT_GENOME000634_11', '25913', 'GUT_GENOME001340_8', '45847']:
            #     print (cluster.ref1, cluster.ref1_positions, cluster.ref2, cluster.ref2_positions, len(cluster.support_reads))
            if len(cluster.support_reads) == 0: # ignore the bkp not supported by split reads
                # print ("no reads", cluster.ref1, cluster.ref1_positions, cluster.ref2, cluster.ref2_positions)
                continue
            choose_acc_from_cluster(cluster)
            bkp_num_support += 1
        # break
        # print (bkp_num_support)
    print ('number of bkp with support reads is %s.'%(bkp_num_support), bkp_num)

def choose_acc_from_cluster(cluster):
    flag = False
    locus_score = {}
    score1 = 0
    score2 = 0
    inte = 2 * rlen # search with a larger interval
    cluster.sort_support_reads()
    for readobj in cluster.support_reads:
        # if readobj.qname == "GUT_GENOME036853_4-1172_2":
        #     print ("*******",readobj.qname, reads_mapped_len[readobj.qname], readobj.end_point)
        # if reads_mapped_len[readobj.qname] < rlen - tolerate_read_mismatch_num:
        #     continue
        if readobj.end_point == True: #the pos is near the segment end, so may be false positive
            continue
        score1 = 0
        my_pos1 = 0
        score2 = 0
        my_pos2 = 0
        from_side = ''
        to_side = ''
        
        if cluster.direction == 'True':
            extract_ref_direction = 'right'
        else:
            extract_ref_direction = 'left'
        #for right clipped seq, if the seg is reverse-complement, extract seq from left to the breakpoint.
        #else, we extract seq from the breakpoint to right
        # test = ['GUT_GENOME000634_11', '25913', 'GUT_GENOME001340_8', '45847']
        # if readobj.ref1 in test and readobj.ref2 in test:

        #     print (readobj.qname, readobj.ref1, readobj.pos1, readobj.ref2, readobj.pos2,\
        #         readobj.mapped_len, readobj.clipped_direction, readobj.end_point, readobj.clipped,\
        #          readobj.direction_confict, readobj.ref2_clipped_direction)

        read_seq = readobj.seq1
        read_seq_len = len(read_seq)
        if read_seq_len > min_seq_len and readobj.clipped == 1:
            for possible_bkp in range(cluster.ref1_positions[0]-inte, cluster.ref1_positions[-1] + inte):
                if readobj.clipped_direction == extract_ref_direction:
                    ref_seq = extract_ref_seq(cluster.ref1, possible_bkp-read_seq_len, possible_bkp)
                else:
                    ref_seq = extract_ref_seq(cluster.ref1, possible_bkp, possible_bkp + read_seq_len)
                if cluster.direction == 'True':
                    ref_seq = get_reverse_complement_seq(ref_seq)
                if readobj.clipped_direction == 'right':
                    to_side = 'left'
                    if cluster.direction == 'True':
                        from_side = 'left'
                    else:
                        from_side = 'right'
                else:
                    to_side = 'right'
                    if cluster.direction == 'True':
                        from_side = 'right'
                    else:
                        from_side = 'left'
                matches = compute_scores(read_seq, ref_seq)/read_seq_len
                if matches > score1 and matches > min_match_score: #ssw alignment
                    score1 = matches
                    cluster.pos1 =  possible_bkp  
                    if readobj.real_ref == cluster.ref2:
                        cluster.pos2 =  readobj.pos2 
                    acc1 = Acc_Bkp(cluster, from_side, to_side, read_seq, ref_seq, score1)
        read_seq = readobj.seq2
        read_seq_len = len(read_seq)
        if read_seq_len > min_seq_len and readobj.clipped == 2:
            for possible_bkp in range(cluster.ref2_positions[0]-inte, cluster.ref2_positions[-1]+inte): #cluster.ref2_positions: 
                if readobj.clipped_direction == extract_ref_direction:
                    ref_seq = extract_ref_seq(cluster.ref2, possible_bkp-read_seq_len, possible_bkp)
                else:
                    ref_seq = extract_ref_seq(cluster.ref2, possible_bkp, possible_bkp + read_seq_len)
                if cluster.direction == 'True':
                    ref_seq = get_reverse_complement_seq(ref_seq)

                if readobj.clipped_direction == 'right':
                    from_side = 'left'
                    if cluster.direction == 'True':
                        to_side = 'left'
                    else:
                        to_side = 'right'
                else:
                    from_side = 'right'
                    if cluster.direction == 'True':
                        to_side = 'right'
                    else:
                        to_side = 'left'

                matches = compute_scores(read_seq, ref_seq)/read_seq_len
                if matches > score2 and matches > min_match_score: #ssw alignment
                    score2 = matches
                    cluster.pos2 = possible_bkp
                    if readobj.real_ref == cluster.ref1:
                        cluster.pos1 =  readobj.pos1   
                    acc2 = Acc_Bkp(cluster, from_side, to_side, read_seq, ref_seq, score2)

        if cluster.pos1 > 0 and cluster.pos2 > 0:
            if score1 > min_match_score and acc1.recheck():
                acc_bkp_list.append(acc1)
            elif score2 > min_match_score and acc2.recheck():
                acc_bkp_list.append(acc2)
            # if readobj.ref1 in test and readobj.ref2 in test:
            #     print (readobj.qname, readobj.ref1, readobj.pos1, readobj.ref2, readobj.pos2,\
            #         readobj.mapped_len, readobj.clipped_direction, score1, score2, read_seq,\
            #          ref_seq, read_seq_len, reads_mapped_len[readobj.qname])
            break #keep searching accurate bkp until the acc pos is found.

class Acc_Bkp(object):
    def __init__(self, cluster, from_side, to_side, read_seq, ref_seq, score):
        self.from_ref = cluster.ref1
        self.to_ref = cluster.ref2
        self.from_bkp = cluster.pos1
        self.to_bkp = cluster.pos2
        self.if_reverse = cluster.direction
        self.from_side = from_side
        self.to_side = to_side
        self.read_str = read_seq
        self.ref_str = ref_seq
        self.similarity = round(score, 3)
        self.refs_sim = 0
        self.max_refs_sim = 0.4

    def print_out(self):
        print (self.from_ref, self.from_bkp, self.to_ref, self.to_bkp, self.from_side,\
         self.to_side, self.if_reverse, self.similarity)

    def write_out(self, writer):
        writer.writerow ([self.from_ref, self.from_bkp, self.to_ref, self.to_bkp, \
        self.from_side, self.to_side, self.if_reverse, self.read_str, self.ref_str, self.similarity])

    def compare_two_refs(self):
        check_len = 50
        from_ref_seq = extract_ref_seq(self.from_ref, self.from_bkp-check_len, self.from_bkp + check_len)         
        to_ref_seq = extract_ref_seq(self.to_ref, self.to_bkp-check_len, self.to_bkp + check_len)   

        matches1 = compute_scores(from_ref_seq, to_ref_seq)/len(from_ref_seq)
        from_ref_seq = get_reverse_complement_seq(from_ref_seq)
        matches2 = compute_scores(from_ref_seq, to_ref_seq)/len(from_ref_seq)
        if matches1 > matches2:
            score = matches1
        else:
            score = matches2
        self.refs_sim = round(score, 3)

    def recheck(self):
        self.compare_two_refs()
        if self.refs_sim > self.max_refs_sim:
            return False
        else:
            return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get accurate hgt breakpoints", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-g", type=str, help="<str> Original Metagenomic reference", metavar="\b")
    required.add_argument("-a", type=str, help="<str> Raw breakpoints file.", metavar="\b")
    required.add_argument("-o", type=str, help="<str> output file of accurate breakpoints.", metavar="\b")
    required.add_argument("-u", type=str, help="<str> unique reads bam file.", metavar="\b")
    required.add_argument("-s", type=str, help="<str> split reads bam file.", metavar="\b")
    optional.add_argument("-t", type=int, default=5, help="<int> number of threads", metavar="\b")
    optional.add_argument("-n", type=int, default=1, help="<0/1> 1 indicates the aligned-ref is extracted.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())

    unique_bam_name = args["u"]
    unique_bamfile = pysam.AlignmentFile(filename = unique_bam_name, mode = 'r')
    mean, sdev, rlen = getInsertSize(unique_bamfile)
    insert_size = int(mean + 2*sdev)
    rlen = int(rlen)

    raw_bkp = args["a"]
    split_bam_name = args["s"]
    ref = args["g"]
    output_acc_bkp_file = args["o"]


    if len(sys.argv) == 1:
        print (Usage%{'prog':sys.argv[0]})
    else:
        reads_mapped_len = {}
        acc_bkp_list = []
        rrm = Read_Raw_Bkp(raw_bkp)
        ref_fasta = Fasta(ref)
        rrm.cluster_bkp()
        rrm.sort_cluster()


        read_split_bam(split_bam_name)
        find_accurate_bkp()

        f = open(output_acc_bkp_file, 'w', newline='')
        writer = csv.writer(f)
        header = ['from_ref','to_ref','from_pos','to_pos','from_side','to_side',\
        'if_reverse','read_seq','ref_seq','similarity']
        writer.writerow(header)
        for acc in acc_bkp_list:
            # acc.print_out()
            acc.write_out(writer)
        f.close()
        print ('Final bkp num is %s'%(len(acc_bkp_list)))


        #edit the cluster raw bkp parameter.
        #edit the reads supporting parameter in add_support_split_reads