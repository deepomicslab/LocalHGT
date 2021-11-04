"""
To solve FPR problem, extract reads that support breakpoints, map the reads to the breakpoint region, and get precise results.
"""

import sys
import os
import pysam
from sklearn.cluster import DBSCAN
import numpy as np
import re
import logging


cigar_dict = {0:'M',1:'M',2:'M',3:'M',4:'N',5:'N'}

def read_BP_and_reads(raw_bp_file):
    print ('start read')
    break_point_pair_dict = {}
    chrom_names = set()

    locus_list = []
    reads_ID_list = []
    f = open(raw_bp_file)
    extract_region = 2000
    index = 0
    for line in f:
        array = line.strip().split()
        if array[0] == 'Point':
        # reads_ID_list += array[:-5]
            break_point_pair_dict[index] = [array[1], int(array[2]), array[1], int(array[2])]
            chrom_names.add(array[1])
            index += 1
        elif array[0] == 'Reads':
            reads_ID_list.append(int(array[1])) 

    point_list, no_use = genome_coordinate(break_point_pair_dict, list(chrom_names))
    print ('raw break points clustering is done.', len(point_list))
    # print (point_list)
    chrom_points_dict = {}
    index = 0
    for p in point_list:
        if index % 2 == 1:
            index += 1
            continue
        # if p[0] == 'NZ_NQMH01000001.1':
        #     print (p[1])
        if p[0] in chrom_points_dict.keys():
            chrom_points_dict[p[0]].append(int(p[1]))
        else:
            chrom_points_dict[p[0]] = [int(p[1])]
        index += 1

    # point_list = []
    for chrom in chrom_points_dict.keys():
        locus = sorted(chrom_points_dict[chrom])
        locus_region = []
        for point in locus:
            pos1 = int(point)
            start1 = pos1 - extract_region
            if start1 <= 0 :
                start1 = 1
            locus_region.append([start1, pos1 + extract_region])
            if len(locus_region) > 1:
                flag = True
                while flag:
                    flag = False
                    for i in range(len(locus_region)-1):
                        if locus_region[i+1][0] >= locus_region[i][0] and locus_region[i+1][0] <= locus_region[i][1]:
                            locus_region[i][1] = locus_region[i+1][1]
                            del locus_region[i+1]
        for loc in locus_region:
            locus_list.append('%s:%s-%s'%(chrom, loc[0], loc[1]))
    print (len(locus_list))
        # print (chrom,locus_region )
    #     print (locus) 
    #     # if len(locus) == 1:
    #     start1 = locus[0] - extract_region
    #     if start1 <= 0 :
    #         start1 = 1  
    #     locus_list.append('%s:%s-%s'%(chrom, start1, locus[0] + extract_region))   
    #     past_end = locus[0] + extract_region
    #     if len(locus) > 1:
    #         for i in range(1, len(locus)):
    #             if locus[i] 
    #             start1 = locus[i] - extract_region
    #             if start1 <= 0 :
    #                 start1 = 1  
    #             if locus[i] - locus[i - 1] <= extract_region:
    #                 locus_list.append('%s:%s-%s'%(chrom, start1, locus[i] + extract_region))   
    #             else:
    #                 locus_list.append('%s:%s-%s'%(chrom, start1, locus[i-1] + extract_region)) 


    # print (point_list)
    # for point in point_list:
    #     ###
    #     pos1 = int(point[1])
    #     start1 = pos1 - extract_region
    #     if start1 <= 0 :
    #         start1 = 1
    #     locus_list.append('%s:%s-%s'%(point[0], start1, pos1 + extract_region))
    return locus_list, reads_ID_list  
   
def reads_ID_filter(reads_ID_list):
    new_reads_ID_list = set()
    for ID in reads_ID_list:
        new_reads_ID_list.add(int(int(ID-1)/2))
    return (new_reads_ID_list)

def extract_reads_by_ID(fq1, fq2, reads_ID_list, new_fq1, new_fq2):
    new_reads_ID_list = reads_ID_filter(reads_ID_list)
    # print (new_reads_ID_list)

    out1 = open(new_fq1, 'w')
    out2 = open(new_fq2, 'w')

    reads_ID = 0
    fq1_f = False
    fq2_f = False
    with open(fq1, 'r') as f1, open(fq2, 'r') as f2:    
        i = 0
        for line1 in f1:
            line2 = f2.readline()
            if i % 4 == 0:
                #fq1
                if reads_ID in new_reads_ID_list:
                    fq1_f = True
                reads_ID += 1
            if fq1_f:
                print (line1, end = '', file = out1)
                print (line2, end = '', file = out2)
            if i % 4 == 3:
                fq1_f = False
            i += 1 
    out1.close()
    out2.close()

def extract_ref(ref, locus_list, extracted_ref, ID):
    # ID = 'test'  
    os.system('samtools faidx %s'%(ref))
    max_extract = 1000 
    if len(locus_list) < max_extract:
        order = 'samtools faidx %s '%(ref)
        for locus in locus_list:
            order += locus
            order += ' '
        order += '>%s'%(extracted_ref)
        # print (order)
        os.system(order)
    else:
        
        for i in range(int(len(locus_list)/max_extract) + 1):
            order = 'samtools faidx %s '%(ref)
            tmp_file = extracted_ref + '.tmp%s.txt'%(i)
            for j in range(max_extract * i, max_extract * (i+1)):
                if j >= len(locus_list):
                    continue
                locus = locus_list[j]
                order += locus
                order += ' '
            order += '>%s'%(tmp_file)
            os.system(order)
            print ('%s fragments were extracted.'%(i*max_extract+max_extract))
        #merge fragments
        os.system('cat %s.tmp*.txt >%s'%(extracted_ref, extracted_ref))
        os.system('rm %s.tmp*.txt'%(extracted_ref))
    print ('Reference extraction is done.')
    
def bwa_map(extracted_ref, new_fq1, new_fq2, bam):
    #-F 4
    order = """
    ref=%s
    #ref='/mnt/d/breakpoints/meta_simu/01.refs/640_reference.fa'
    fq1=%s
    fq2=%s
    bam=%s
    bwa index $ref
    bwa mem $ref $fq1 $fq2| samtools view -bS |samtools sort -o $bam
    samtools index $bam  
    """%(extracted_ref, new_fq1, new_fq2, bam)
    os.system(order)

def check_junction_pair(bam, bkp_pair_list):
    support_reads = []
    for i in range(len(bkp_pair_list)):
        support_reads.append([])
    samfile = pysam.AlignmentFile(bam, "rb")
    near = 1000
    print ('###################start junction.')
    for read in samfile:
        if read.is_unmapped == False and read.mate_is_unmapped == False and read.flag < 2048 and not read.is_duplicate: 
            # if read.qname == 'SRR5024277.1.20723410':
            #     print (read.qname, read.cigar, read.mapping_quality)
            if read.reference_name == read.next_reference_name and abs(read.reference_start - read.next_reference_start) < 1000:
                continue
            if read.has_tag('XA') or read.mapping_quality < 20:
                continue
            
            no_clip = True
            for ci in read.cigar:
                if ci[0] == 4 or ci[0] == 5:
                    no_clip = False
            if no_clip == False:
                continue
            
            # print (read, len(read.query_sequence), read.reference_length, read.query_alignment_length)
                # print (read, read.qname, read.reference_name, read.next_reference_name, read.next_reference_start, read.get_tag('XA'))
            chrom_name = read.reference_name.split(':')[0]
            if len(read.reference_name.split(':')) == 2:
                chrom_pos = int(read.reference_name.split(':')[1].split('-')[0]) + read.reference_start
            else:
                chrom_pos = read.reference_start
            # chrom_pos = int(read.reference_name.split(':')[1].split('-')[0]) + read.reference_start
            next_name = read.next_reference_name.split(':')[0]
            # next_pos = int(read.next_reference_name.split(':')[1].split('-')[0]) + read.next_reference_start
            if len(read.next_reference_name.split(':')) == 2:
                next_pos = int(read.next_reference_name.split(':')[1].split('-')[0]) + read.next_reference_start
            else:
                next_pos = read.next_reference_start
            for i in range(len(bkp_pair_list)):
                bk_pair = bkp_pair_list[i]
                if bk_pair[0] == chrom_name and bk_pair[2] == next_name:
                    if abs(bk_pair[1] - chrom_pos) < near and abs(bk_pair[3] - next_pos) < near:
                        if read.qname not in support_reads[i]:
                            support_reads[i].append(read.qname)
                elif bk_pair[2] == chrom_name and bk_pair[0] == next_name:
                    if abs(bk_pair[3] - chrom_pos) < near and abs(bk_pair[1] - next_pos) < near:
                        if read.qname not in support_reads[i]:
                            support_reads[i].append(read.qname)
            # print (chrom_name, chrom_pos, next_name, next_pos, read)
    return support_reads

def blast_map2_self(a, b):
    makedb = 'makeblastdb -in %s -dbtype nucl -out tmp_ref'%(a)
    blastn = 'blastn -query %s -outfmt 6 -out tmp_blast.out -db tmp_ref'%(b)
    os.system(makedb)
    os.system(blastn)

def check_spilt_reads(bam, split_bam, filter_sv, original_ref):
    # ref = '/mnt/d/breakpoints/big/gut.reference.filter.fa'
    filer_f = open(filter_sv, 'w')
    chrom_names = set()
    break_point_pair_dict = {}
    break_point_reads_locus_dict = {}
    map_region = 500
    # reads_len = 250
    samfile = pysam.AlignmentFile(split_bam)
    read_num = 0
    # samfile = pysam.AlignmentFile(bam, "rb")
    for read in samfile.fetch():
        read_num += 1
        if read.is_unmapped or read.mate_is_unmapped or read.is_duplicate:
            continue
        # if read.mapping_quality < 10:
        #     continue
        if not read.has_tag('SA'): #only include split read
            continue   
        if len(read.get_tag('SA').split(';')) > 2:   #make sure only one supplementary alignment
            continue
  
        # uncoverd_length = split_complementary(read.cigarstring, read.get_tag('SA').split(',')[3])
        # if uncoverd_length > 50: #20
        #     continue

        ignore = False
        bp_f = False
        bP_pos = 0
        real_locus_pos, real_pos = 0, 0
        # print (read.reference_name.split(':')[0], read.next_reference_name.split(':')[0], read, read.cigar)
        if len(read.cigar) == 0:
            continue
        start = cigar_dict[read.cigar[0][0]]
        for ci in read.cigar:
            if cigar_dict[ci[0]] != start:
                real_pos = bP_pos # to record the break locus of the reads
                if start == 'M':
                    real_locus_pos = bP_pos  #to find the exact break point locus
                else:
                    real_locus_pos = 0               
                if bp_f == False:
                    bp_f = True
                else:
                    ignore = True
                start = cigar_dict[ci[0]]
                # break
            # if start == 'M':
            bP_pos += ci[1]
        # if bp_f and ignore == False:  
        if True:        
            break_point_pos = read.reference_start + real_locus_pos
            chrom_name = read.reference_name.split(':')[0]
            if len(read.reference_name.split(':')) == 2:
                chrom_pos = int(read.reference_name.split(':')[1].split('-')[0]) + break_point_pos
            else:
                chrom_pos = break_point_pos
            chrom_names.add(chrom_name)
            #update read name
            # full_name = read.query_name
            if read.is_read1:
                full_name = read.query_name + '/1'
            else:
                full_name = read.query_name + '/2'
            if full_name in break_point_pair_dict.keys():
                break_point_pair_dict[full_name] += [chrom_name, chrom_pos]
                break_point_reads_locus_dict[full_name] += [real_pos]
            else:
                break_point_pair_dict[full_name] = [chrom_name, chrom_pos]
                break_point_reads_locus_dict[full_name] = [real_pos]

    print ('Number of raw split reads pair:', len(break_point_reads_locus_dict), 'read num is:', read_num)
    new_break_point_pair_dict = {}

    for read_name in break_point_pair_dict:
        if len(break_point_reads_locus_dict[read_name]) < 2:
            # print ('hi', read_name, break_point_reads_locus_dict[read_name], break_point_pair_dict[read_name])
            continue
        elif len(break_point_reads_locus_dict[read_name]) > 2:
            continue
        else:
            new_break_point_pair_dict[read_name] = break_point_pair_dict[read_name]
            #print (read_name, break_point_pair_dict[read_name], break_point_reads_locus_dict[read_name])
    print ('Number of split reads:', len(new_break_point_pair_dict))
    # print (new_break_point_pair_dict)
    logging.info("Number of considered split reads is %s."%(len(new_break_point_pair_dict)))
    break_point_pair_dict = new_break_point_pair_dict
    point_list, bkp_pair_list, support_split_reads_num_list = genome_coordinate(break_point_pair_dict, list(chrom_names))
    # print ('DBSCAN clustering is done.', bkp_pair_list)
    support_reads = check_junction_pair(bam, bkp_pair_list)
    species_dict = scaffold_species_name(original_ref)  #The scaffold's corresponding species name.
    result_dict = {}
    species_index = {}
    j = 0
    filter_support_reads = 0
    filter_hit = 0
    filter_same_species = 0
    final_num = 0 
    print ('raw break point num is:', len(bkp_pair_list))
    logging.info("Number of raw bkp is %s."%(len(bkp_pair_list)))
    for i in range(len(bkp_pair_list)):
    # for bkp_pair in bkp_pair_list:
        bkp_pair = bkp_pair_list[i]
        support_split_reads_num = support_split_reads_num_list[i]
        # print (bkp_pair[0], bkp_pair[1], bkp_pair[2], bkp_pair[3])
        if bkp_pair[0] == bkp_pair[2] and abs(bkp_pair[1] - bkp_pair[3]) < 500:
            filter_same_species += 1
            continue
        if bkp_pair[0] == bkp_pair[2]:
            filter_same_species += 1
            continue
        species_1 = species_dict[bkp_pair[0]]
        species_2 = species_dict[bkp_pair[2]]
        if species_1 == species_2:
            filter_same_species += 1
            continue
        
        # print ('raw bkp:', bkp_pair, support_split_reads_num, len(support_reads[i]))
        if support_split_reads_num < 1 or len(support_reads[i]) < 2:
            # print ('hi', support_split_reads_num, len(support_reads[i]))
            filter_support_reads += 1
            continue
        """
        for j in range(2):
            start = bkp_pair[2*j+1] - map_region
            end = bkp_pair[2*j+1] + map_region
            if start < 0:
                start = 1
            os.system('samtools faidx %s %s:%s-%s>tmp_blast.%s.fa'%(original_ref,bkp_pair[2*j],start,end,j))
        blast_map2_self('tmp_blast.0.fa', 'tmp_blast.1.fa')
        hit_count = len(open('tmp_blast.out','rU').readlines())
        if hit_count > 0:
            filter_hit += 1
            # print (bkp_pair)
            continue
        """
        final_num += 1
        ##sort the results
        if species_1 not in species_index:
            species_index[species_1] = j
            j += 1
        if species_2 not in species_index:
            species_index[species_2] = j
            j += 1   
        if species_index[species_1] > species_index[species_2]:
            species_key = species_2 + '@' + species_1
            save = bkp_pair[0]
            bkp_pair[0]=bkp_pair[2]
            bkp_pair[2]=save
            save = bkp_pair[1]
            bkp_pair[1]=bkp_pair[3]
            bkp_pair[3]=save
        else:
            species_key = species_1 + '@' + species_2
        result = bkp_pair[0]+'\t'+str(bkp_pair[1])+'\t'+bkp_pair[2]+'\t'+str(bkp_pair[3])+'\t'+str(support_split_reads_num)+'\t'+str(len(support_reads[i]))+'\t'+species_key#+'\t'+str(support_reads[i])#+'\t'+species_1+'\t'+species_2
        
        if species_key not in result_dict:
            result_dict[species_key] = [result]
        else:
            result_dict[species_key].append(result)
        # print (bkp_pair[0], bkp_pair[1], bkp_pair[2], bkp_pair[3], support_split_reads_num, len(support_reads[i]), hit_count, species_1, species_2)
    for species_key in result_dict.keys():
        for result in result_dict[species_key]:
            print (result, file = filer_f)
    logging.info("Remove %s bkp by supporting reads."%(filter_support_reads))
    logging.info("Remove %s bkp by blast."%(filter_hit))
    logging.info("Remove %s bkp by same species."%(filter_same_species))
    logging.info("We get %s bkp finally."%(final_num))
    print (filter_support_reads, filter_hit, final_num)
    filer_f.close()
    
def genome_coordinate(break_point_pair_dict, chrom_names):
    
    # chrom_names = []
    # for read_name in break_point_pair_dict:
    #     pair = break_point_pair_dict[read_name]
    #     chrom_names.append(pair[0])
    # chrom_names = list(break_point_pair_dict.keys())
    sort_chrom_name = {}
    for i in range(len(chrom_names)):
        sort_chrom_name[chrom_names[i]] = i
    genome_coordinate_points = {}
    for read_name in break_point_pair_dict:
        bkp_pair = break_point_pair_dict[read_name]
        if len(bkp_pair) < 4:
            continue
        # print (bkp_pair)

        #sort the break point pair
        if sort_chrom_name[bkp_pair[0]] > sort_chrom_name[bkp_pair[2]]:
            bkp_pair = [bkp_pair[2], bkp_pair[3], bkp_pair[0],bkp_pair[1]]
        elif sort_chrom_name[bkp_pair[0]] == sort_chrom_name[bkp_pair[2]]:
            if bkp_pair[1] > bkp_pair[3]:
                 bkp_pair = [bkp_pair[0], bkp_pair[3], bkp_pair[2],bkp_pair[1]]
            
        coordinate = str(bkp_pair[0]) + '@' + str(bkp_pair[2])
        if coordinate not in genome_coordinate_points:
            genome_coordinate_points[coordinate] = [[bkp_pair[1], bkp_pair[3]]] 
        else:
            genome_coordinate_points[coordinate] += [[bkp_pair[1], bkp_pair[3]]] 
    # print (read_name, break_point_pair_dict[read_name])

    point_list = []
    bkp_pair_list = []
    support_split_reads_num_list = []
    for co in genome_coordinate_points:        
        # print (np.array(genome_coordinate_points[co]))
        cluster_dict = {}
        clustering = DBSCAN(eps=200, min_samples=1).fit(np.array(genome_coordinate_points[co]))
        labels = clustering.labels_
        for i in range(len(labels)):
            if labels[i] not in cluster_dict:
                cluster_dict[labels[i]] = [genome_coordinate_points[co][i]]
            else:
                cluster_dict[labels[i]] += [genome_coordinate_points[co][i]]
        chrom_name = co.split('@')
        if chrom_name[0] == chrom_name[1]:
            continue
        for clu in cluster_dict:
            if len(cluster_dict[clu]) < 1:
                # print (cluster_dict[clu])
                continue
            center = np.median(cluster_dict[clu], axis=0)
            # print (clu, cluster_dict[clu], np.median(cluster_dict[clu], axis=0))
            # print (chrom_name[0], int(center[0]), chrom_name[1], int(center[1]))
            bkp_pair_list.append([chrom_name[0], int(center[0]), chrom_name[1], int(center[1])])
            support_split_reads_num_list.append(len(cluster_dict[clu]))
            point_list.append([chrom_name[0], int(center[0])])
            point_list.append([chrom_name[1], int(center[1])])
        # print (co, genome_coordinate_points[co], clustering.labels_)
    
    return point_list, bkp_pair_list, support_split_reads_num_list

def main(ID, ref):

    fq1 = '/mnt/d/breakpoints/meta_simu/03.samples//%s_HGT.1.fq'%(ID)
    fq2 = '/mnt/d/breakpoints/meta_simu/03.samples//%s_HGT.2.fq'%(ID)


    # origin_map = """bwa index %s\nbwa mem -t 10 %s %s %s| samtools view -bS -F 4 |samtools sort -o /mnt/d/breakpoints/meta_simu/03.samples/%s.full.bam"""%(ref, ref, fq1, fq2, ID)
    # os.system(origin_map)

    new_fq1 = '/mnt/d/breakpoints/meta_simu/03.samples//%s_HGT.1.extracted.fq'%(ID)
    new_fq2 = '/mnt/d/breakpoints/meta_simu/03.samples//%s_HGT.2.extracted.fq'%(ID)
    bam='/mnt/d/breakpoints/meta_simu/03.samples/%s.extracted.bam'%(ID)
    filter_sv='/mnt/d/breakpoints/meta_simu/03.samples/%s.result.filter.txt'%(ID)
    # 

    raw_bp_file = '/mnt/d/breakpoints/meta_simu/03.samples/%s.our.result.txt'%(ID)
    locus_list, reads_ID_list = read_BP_and_reads(raw_bp_file)
    

    # locus_list = ['NZ_FLOC01000012.1:65274-68274']
    extracted_ref = '/mnt/d/breakpoints/meta_simu/01.refs/%s_extracted.fasta'%(ID)
    extract_ref(ref, locus_list, extracted_ref, ID)
    extract_reads_by_ID(fq1, fq2, reads_ID_list, new_fq1, new_fq2)
    bwa_map(extracted_ref, new_fq1, new_fq2, bam)
    print ('finish preocessing, start identify HGT:')   
    check_spilt_reads(bam, filter_sv)
    print ('done')

def precise_bp(bam, split_bam, filter_sv, original_ref):
    # check_junction_pair(bam)
    check_spilt_reads(bam, split_bam, filter_sv, original_ref)
    print ('it is done')

def scaffold_species_name(ref):
    name = ref + '.name'
    # name = '/mnt/d/breakpoints/big/gut.reference.filter.fa.name'
    if not os.path.isfile(name):
        f = open(ref)
        out=open(name, 'w')
        for line in f:   
            if line[0] != '>':
                continue
            line = line.strip()
            array =line[1:].split()
            # print (line, array)
            scaffold = array[0]
            species = array[1] + '_' + array[2]   
            print (scaffold, species, file = out)
        out.close()
        f.close()

    species_dict = {}
    f = open(name)
    for line in f:
        array = line.strip().split()
        species_dict[array[0]] = array[1]
    f.close()
    print ('Extract scaffold species name done.')
    return species_dict

def split_complementary(a, b):
    # a = '106S104M3D41M'
    # b = '106M145S'
    r = r"\d+\D"
    a_cigar = re.findall(r, a)
    b_cigar = re.findall(r, b)
    start = 0
    a_M = []
    for ac in a_cigar:
        if ac[-1] == 'D':
            continue
        if ac[-1] != 'S' and  ac[-1] != 'H':
            a_M.append([start, start + int(ac[:-1])])
        start = start + int(ac[:-1])
    # print (a_M)
    start = 0
    b_M = []
    for bc in b_cigar:
        if ac[-1] == 'D':
            continue
        if bc[-1] != 'S' and  bc[-1] != 'H':
            b_M.append([start, start + int(bc[:-1])])
        start = start + int(bc[:-1])
    # print (b_M)
    read_len = start
    ab_M = a_M + b_M
    # sorted_ab_M = []
    # print (ab_M)
    flag = True
    while flag:
        flag = False
        for i in range(len(ab_M) - 1):
            if ab_M[i][0] > ab_M[i+1][0]:
                m = ab_M[i]
                ab_M[i] = ab_M[i+1]
                ab_M[i+1] = m
                flag = True    
    # print (ab_M)
    all_cover = merge_interval(ab_M)
    return read_len - all_cover

def merge_interval(intervals):
    # intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
    all_cover = 0
    for me in merged:
        all_cover += (me[1] - me[0])
    # print (merged, all_cover)
    return all_cover

if __name__ == "__main__":
    # main(sys.argv[1], sys.argv[2])
    # scaffold_species_name('/mnt/d/breakpoints/meta_simu/01.refs/gut.reference.fa')
    split_complementary()
