#!/usr/bin/env python3

import os
import numpy as np
import random
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

def DNA_complement2(sequence):
    sequence = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

def add_snp(sequence, rate):
    locus_set = {}
    allele = ['A', 'T', 'C', 'G']
    ref_len = len(sequence)
    num = int(ref_len* rate)
    i = 0
    #'''
    while i < num:
        random_locus = np.random.randint(ref_len - 1 )
        if random_locus not in locus_set.keys():
            ref = sequence[random_locus]
            while True:
                allele_index = np.random.randint(4)
                alt = allele[allele_index]
                if alt != ref:
                    break
            # print (random_locus, ref, alt)
            sequence = sequence[:random_locus] + alt + sequence[random_locus+1:]
            # locus_set.append(random_locus)
            locus_set[random_locus] = 1
            i += 1
    return sequence

def add_indel(sequence, rate):
    locus_set = {}
    allele = ['A', 'T', 'C', 'G']
    ref_len = len(sequence)
    num = int(ref_len* rate)
    i = 0
    #'''
    while i < num:
        random_locus = np.random.randint(ref_len - 1 )
        if random_locus not in locus_set.keys():
            if np.random.random() < 0.5:
                #50% for deletion
                sequence = sequence[:random_locus] + sequence[random_locus+1:]
            else:
                #50% for insertion
                allele_index = np.random.randint(4)
                insert_allele = allele[allele_index]
                sequence = sequence[:random_locus] + insert_allele + sequence[random_locus:]
            # locus_set.append(random_locus)
            locus_set[random_locus] = 1
            i += 1
    return sequence       

def extract_uniq_region(map_ref):
    remove_dict = {}
    all_scaffold = {}
    line_num = 0
    for line in open(map_ref, 'r'):
        line_num += 1
        if line_num % 1000000 == 0:
            print (line_num)
        # if line_num > 1000000:
        #     break
        array = line.split()

        #record scaffold name
        if array[0] not in all_scaffold.keys():
            all_scaffold[array[0]] = 1
        #     continue
        
        if float(array[6]) < float(array[7]):
            fir_s = float(array[6])
            fir_e = float(array[7])
        else:
            fir_s = float(array[7])
            fir_e = float(array[6])     
        if float(array[8]) < float(array[9]):
            sec_s = float(array[8])
            sec_e = float(array[9])
        else:
            sec_s = float(array[9])
            sec_e = float(array[8])   

        if array[0] == array[1]:
            if fir_s >= sec_s and fir_s <= sec_e:
                continue
            elif fir_e >= sec_s and fir_e <= sec_e:
                continue

        if float(array[3]) < 50: # length
            continue
        if array[0] not in remove_dict.keys():
            remove_dict[array[0]] = []
        qstart = fir_s#float(array[6])
        qend = fir_e#float(array[7])
        overlap = False
        for loci in remove_dict[array[0]]:
            if qstart < loci[0] and qend > loci[0]:
                loci[0] = qstart
                overlap = True
            elif qstart < loci[1] and qend > loci[1]:
                loci[1] = qend
                overlap = True
            elif qstart >= loci[0] and qend <= loci[1]:
                overlap = True
        if overlap == False:
            remove_dict[array[0]].append([qstart, qend])
    print ("start sort segs...", len(all_scaffold.keys()))
    i = 0
    for scaffold in list(all_scaffold.keys()):
        #sometimes there is no repeat region in the scaffold, it not in remove_dict's keys.
        if scaffold in remove_dict.keys():
            remove_dict[scaffold] = sort_remove_loci(remove_dict[scaffold], scaffold)
        else:
            remove_dict[scaffold] = [[1,'end']]
        i += 1
        if i % 10000 == 0:
            print ("sort", i)
    uniq_segs_loci = remove_dict
    return uniq_segs_loci

def remove_overlap(sorted_remove_loci):
    flag = True
    while flag:
        flag = False
        new_sorted_remove_loci = []
        i = 0
        while i < len(sorted_remove_loci):
            if i + 1 < len(sorted_remove_loci) and sorted_remove_loci[i+1][0] >= sorted_remove_loci[i][0] and sorted_remove_loci[i+1][0] <= sorted_remove_loci[i][1]:
                flag = True
                if sorted_remove_loci[i][1] >= sorted_remove_loci[i+1][1]:
                    new_sorted_remove_loci.append(sorted_remove_loci[i])
                else:
                    new_sorted_remove_loci.append([sorted_remove_loci[i][0], sorted_remove_loci[i+1][1]])
                i += 2
            else:
                new_sorted_remove_loci.append(sorted_remove_loci[i])
                i += 1
        sorted_remove_loci = new_sorted_remove_loci[:]
    return new_sorted_remove_loci

def sort_remove_loci(remove_loci, scaffold):
    sorted_remove_loci = []
    for loci in remove_loci:
        sorted_remove_loci.append(loci)
    flag = True
    while flag:
        flag = False
        for i in range(len(sorted_remove_loci) - 1):
            if sorted_remove_loci[i][0] > sorted_remove_loci[i+1][0]:
                flag = True
                m = sorted_remove_loci[i]
                sorted_remove_loci[i] = sorted_remove_loci[i+1]
                sorted_remove_loci[i+1] = m
    # print ('s',sorted_remove_loci)
    sorted_remove_loci = remove_overlap(sorted_remove_loci)
    start = 1
    uniq_segs_loci = []
    for loci in sorted_remove_loci:
        if int(loci[0]) - start > 1000:
            uniq_segs_loci.append([start, int(loci[0])])
        start = int(loci[1])
    uniq_segs_loci.append([start, 'end'])
    return uniq_segs_loci

def read_fasta(file):
    seq_dict = {}
    fasta_sequences = SeqIO.parse(open(file),'fasta')
    for record in fasta_sequences:
        seq_dict[str(record.id)] = str(str(record.seq))
    #     scaffold_name.append(str(record.id))
    return seq_dict

def blast_map2_self(file, map_ref):
    makedb = 'makeblastdb -in %s -dbtype nucl -out ref'%(file)
    blastn = 'blastn -query %s -outfmt 6 -out %s -db ref'%(file, map_ref)
    os.system(makedb)
    os.system(blastn)

def random_HGT(pa):
    truth_HGT = open(pa.outdir +'/%s.true.sv.txt'%(pa.sample), 'w')

    new_seq_dict = pa.seq_dict.copy()
    chrom_list = list(pa.seq_dict.keys())
    chrom_num = len(chrom_list)
    length_change_chroms = {}
    i = 0
    while i < pa.HGT_num:
        first_chrom = np.random.randint(chrom_num)
        second_chrom = np.random.randint(chrom_num)

        # print (first_chrom, second_chrom)
        while first_chrom == second_chrom or pa.species_dict[chrom_list[first_chrom]] == pa.species_dict[chrom_list[second_chrom]]:
            second_chrom = np.random.randint(chrom_num)
        first_chrom = chrom_list[first_chrom]
        second_chrom = chrom_list[second_chrom]

        if first_chrom in length_change_chroms.keys() or second_chrom in length_change_chroms.keys():  #one chrom only envolved in one HGT.
            continue

        first_seq = new_seq_dict[first_chrom]
        second_seq = new_seq_dict[second_chrom]
        if len(first_seq) < 100000 or len(second_seq) < 100000:
            continue

        #make sure no N near the four points.
        first_uniq_region_s, first_uniq_region_e =  uniq_seg(len(pa.seq_dict[first_chrom]), pa.uniq_segs_loci[first_chrom])
        second_uniq_region_s, second_uniq_region_e =  uniq_seg(len(pa.seq_dict[second_chrom]), pa.uniq_segs_loci[second_chrom])

        if first_uniq_region_s + first_uniq_region_e == 0 or second_uniq_region_s + second_uniq_region_e == 0: #one more time
            return 0

        length_change_chroms[first_chrom] = 1
        length_change_chroms[second_chrom] = 1
        if np.random.random() < 0.5: #reverse the seq
            reverse_flag = True
            new_seq_dict[first_chrom] = first_seq[:first_uniq_region_s] + DNA_complement2(second_seq[second_uniq_region_s:second_uniq_region_e])\
             + first_seq[first_uniq_region_s:]
        else:
            reverse_flag = False
            new_seq_dict[first_chrom] = first_seq[:first_uniq_region_s] + second_seq[second_uniq_region_s:second_uniq_region_e]\
             + first_seq[first_uniq_region_s:]  #insert seg from chrom2 to chrom1
        # new_seq_dict[second_chrom] = second_seq[:second_uniq_region_s] + second_seq[second_uniq_region_e:]
        if pa.donor_in_flag == False:
            del new_seq_dict[second_chrom]

        print (i, first_chrom, first_uniq_region_s, second_chrom, second_uniq_region_s, \
        second_uniq_region_e, i, chrom_num)
        print (first_chrom, first_uniq_region_s, second_chrom, second_uniq_region_s,\
         second_uniq_region_e, reverse_flag, file = truth_HGT)

        i += 1

    #we should add mutations to each scaffold.
    t0 = time.time()
    for chrom in new_seq_dict.keys():
        new_seq_dict[chrom] = add_snp(new_seq_dict[chrom], pa.snp_rate)
        new_seq_dict[chrom] = add_indel(new_seq_dict[chrom], pa.indel_rate)
    t1 = time.time()
    generate_fastq(new_seq_dict, truth_HGT, pa)
    print (t1-t0, time.time()-t1)
    truth_HGT.close()
    return 1

def generate_fastq(new_seq_dict, truth_HGT, pa):
    fasta_file = pa.outdir+'/%s.true.fasta'%(pa.sample)
    fasta = open(fasta_file, 'w')
    i = 0
    sequencing_method = {75:'NS50', 100:'HS20', 150:'HS25'}
    for chrom in new_seq_dict.keys():
        print ('>%s'%(chrom), file = fasta)
        print (new_seq_dict[chrom], file = fasta)

        tmp_file = pa.outdir + '/%s.HGT.tmp.fasta'%(pa.sample) 
        tmp_f = open(tmp_file, 'w')
        print ('>%s'%(chrom), file = tmp_f)
        print (new_seq_dict[chrom], file = tmp_f)   
        tmp_f.close()
        ######### random depth ###########
        # dp = np.random.randint(10, 50)

        fq = 'art_illumina -ss %s --noALN -p -i %s -l %s -m 350 -s 10 --fcov %s -o %s/%s_HGT_%s_tmp.'\
        %(sequencing_method[pa.reads_len], tmp_file, pa.reads_len, pa.depth, pa.outdir, pa.sample, i)
        # fq = 'wgsim -e 0 -N %s -r 0.00 -S 1 -1 150 -2 150 %s /mnt/d/breakpoints/meta_simu/03.samples//%s_HGT_%s_tmp.1.fq /mnt/d/breakpoints/meta_simu/03.samples//%s_HGT_%s_tmp.2.fq'%(int(len(new_seq_dict[chrom])/6), tmp_file, ID, i, ID, i)
        os.system(fq)
        i += 1
    fasta.close()
    os.system('cat %s/%s_HGT_*_tmp.1.fq >%s/%s.1.fq'%(pa.outdir, pa.sample, pa.outdir, pa.sample))
    os.system('cat %s/%s_HGT_*_tmp.2.fq >%s/%s.2.fq'%(pa.outdir, pa.sample, pa.outdir, pa.sample))
    os.system('rm %s/%s_HGT_*_tmp.*.fq'%(pa.outdir, pa.sample))
    os.system('rm %s'%(tmp_file))

def uniq_seg(chrom_len, chrom_uniq_segs):
    flag = True
    i = 0
    max_HGT_len = 55000
    while flag:
        # print (chrom_len)
        HGT_len = np.random.randint(500,max_HGT_len)
        if HGT_len >= chrom_len:
            print ('HGT_len >= chrom_len')
        locus = np.random.randint(500, chrom_len - HGT_len)
        
        left_flag = False
        right_flag = False
        for segs in chrom_uniq_segs:
            #print (locus, segs)
            if segs[1] == 'end':
                segs[1] = chrom_len 

            if locus > segs[0] + 200 and locus < segs[1] - 200:
                left_flag = True
            if locus + HGT_len > segs[0] + 200 and locus + HGT_len < segs[1] - 200:
                right_flag = True

            if left_flag and right_flag:
                return locus, locus + HGT_len    
        i += 1
        if i > 1000000:
            print ('iterate too many times!')
            # break
            return 0, 0

def UHGG_snp(uniq_segs_loci): 
    species_dict = {} 
    pa = Parameters()
    pa.get_dir("/mnt/d/breakpoints/HGT/uhgg_snp/")
    pa.add_segs(uniq_segs_loci)
    pa.get_uniq_len()

    for snp_rate in pa.snp_level:
    # for snp_rate in [0.06, 0.07, 0.08, 0.09]:
        pa.change_snp_rate(snp_rate)
        for index in range(pa.iteration_times):
            pa.get_ID(index)
            if_success = 0
            while if_success == 0:
                ############random select scaffold###########
                all_ref = pa.outdir + '/%s.fa'%(pa.sample)
                fasta_sequences = SeqIO.parse(open(pa.origin_ref),'fasta')       
                f = open(all_ref, 'w')
                select_num = 0
                for record in fasta_sequences:
                    if len(record.seq) < pa.min_genome:
                        continue
                    if pa.uniq_len[str(record.id)] < pa.min_uniq_len:
                        continue
                    if np.random.random() < pa.random_rate:
                        rec1 = SeqRecord(record.seq, id=str(record.id), description="simulation")
                        SeqIO.write(rec1, f, "fasta") 
                        # uniq_segs_loci[str(record.id)] = [[1, len(record.seq)]]
                        species_dict[str(record.id)] = str(record.id)
                        select_num += 1
                    if select_num == pa.scaffold_num + pa.HGT_num:
                        break

                f.close()
                print ('%s scaffolds were extracted.'%(select_num))
                if select_num ==pa.scaffold_num + pa.HGT_num:
                    seq_dict = read_fasta(all_ref)  
                    pa.add_species(species_dict, seq_dict)                  
                    if_success = random_HGT(pa)

def UHGG_cami(): 
    pa = Parameters()
    pa.get_dir("/mnt/d/breakpoints/HGT/uhgg_snp/")

    for snp_rate in pa.snp_level: #[0.02, 0.04]:
        pa.change_snp_rate(snp_rate)
        index = 0
        pa.get_ID(index)
        for level in pa.complexity_level:
            cami_ID = pa.sample + '_' + level
            for j in range(1, 3):
                combine = "cat %s/%s.%s.fq %s/%s.fq/%s.%s.fq >%s/%s.%s.fq"%(pa.outdir, pa.sample, j, \
                pa.cami_dir, pa.cami_data[level], pa.cami_data[level], j, pa.outdir, cami_ID, j)
                print (combine)
                os.system(combine)

def UHGG_cami2(): 
    pa = Parameters()
    pa.get_dir("/mnt/d/breakpoints/HGT/uhgg_snp/")
    com = Complexity()
    for snp_rate in [0.02, 0.04]:
        pa.change_snp_rate(snp_rate)
        index = 0
        pa.get_ID(index)
        for level in pa.complexity_level:
            for j in range(1, 3):
                combine = f"cat {pa.outdir}/{pa.sample}.{j}.fq {com.complexity_dir}/{level}.{j}.fq >{pa.outdir}/{pa.sample}_{level}.{j}.fq"
                print (combine)
                os.system(combine)

class Parameters():
    def __init__(self):
        self.HGT_num = 20
        self.scaffold_num = self.HGT_num
        self.snp_rate = 0.01
        self.indel_rate = self.snp_rate * 0.1  
        self.depth = 100
        self.reads_len = 150
        self.iteration_times = 1
        self.donor_in_flag = False
        self.min_genome = 100000
        self.min_uniq_len = 20000
        self.random_rate = 0.01
        self.cami_data = {'low':'RL_S001__insert_270', 'medium':'RM2_S001__insert_270',\
         'high':'RH_S001__insert_270'}
        self.complexity_level = ['high', 'low', 'medium']
        self.snp_level = [round(x*0.01,2) for x in range(1,6)]
        
        self.uniq_segs_loci = {}
        self.species_dict = {}
        self.sample = ''
        self.outdir = ''
        self.origin_ref = '/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna'
        self.cami_dir = '/mnt/d/breakpoints/HGT/CAMI/'
        self.seq_dict = {}
        self.seq_len = {}
        self.cal_genome_len()
        self.uniq_len = {}

    def add_segs(self, uniq_segs_loci):
        self.uniq_segs_loci = uniq_segs_loci
        
    def add_species(self, species_dict, seq_dict):
        self.species_dict = species_dict
        self.seq_dict = seq_dict

    def change_snp_rate(self, snp_rate):
        self.snp_rate = snp_rate
        self.indel_rate = snp_rate * 0.1

    def get_ID(self, index):
        self.sample = 'species%s_snp%s_depth%s_reads%s_sample_%s'%(self.scaffold_num, \
        self.snp_rate, self.depth, self.reads_len, index)

    def get_dir(self, outdir):
        self.outdir = outdir

    def cal_genome_len(self):
        for line in open(self.origin_ref+'.fai'):
            array = line.split()
            self.seq_len[array[0]] = int(array[1])

    def get_uniq_len(self):
        for sca in self.uniq_segs_loci.keys():
            sca_len = 0
            for interval in self.uniq_segs_loci[sca]:
                start = interval[0]
                end = interval[1]
                if end == 'end':
                    end = self.seq_len[sca]
                sca_len += (end- start) 
            self.uniq_len[sca] = sca_len

class Complexity():
    def __init__(self):
        self.levels = ['low', 'medium', 'high']
        self.size = {'high':1000000000, 'medium':500000000, 'low':100000000}
        self.depths = {'high':10, 'medium':20, 'low':100}
        self.origin_ref = '/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna'
        self.complexity_dir = '/mnt/d/breakpoints/HGT/complexity/'
        self.random_rate = 0.1
        self.read_pair_num = int(10000000000/300)

    def select_genome(self, level):
        all_ref = self.complexity_dir + '/%s.fa'%(level)
              
        f = open(all_ref, 'w')
        select_size = 0
        max_size = self.size[level]
        while select_size < max_size:
            fasta_sequences = SeqIO.parse(open(self.origin_ref),'fasta') 
            for record in fasta_sequences:
                if np.random.random() < self.random_rate:
                    rec1 = SeqRecord(record.seq, id=str(record.id), description="simulation")
                    SeqIO.write(rec1, f, "fasta") 
                    select_size += len(record.seq)
                    if select_size >  max_size:
                        break
        print ("genomes selecting done.")
        self.fastq(all_ref, level)

    def fastq(self, genome, level):
        order = f"wgsim -1 150 -2 150 -r 0.001 -N {self.read_pair_num} {genome} \
        {self.complexity_dir}/{level}.1.fq {self.complexity_dir}/{level}.2.fq"
        print (order)
        os.system(order)

    def run(self):
        # for level in self.levels:
        for level in ['high']:    
            print (level)
            self.select_genome(level)

def generate_complexity():
    com = Complexity()
    com.run()

if __name__ == "__main__":

    t0 = time.time()
    pa = Parameters()
    uniq_segs_file = "/mnt/d/breakpoints/HGT/UHGG/uniq_region_uhgg.npy"
    blast_file = '/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.blast.out'

    # generate_complexity()


    # if os.path.isfile(uniq_segs_file):
    #     uniq_segs_loci = np.load(uniq_segs_file, allow_pickle='TRUE').item()
    # else:
    #     uniq_segs_loci = extract_uniq_region(blast_file)  
    #     np.save(uniq_segs_file, uniq_segs_loci) 
    # t1 = time.time()
    # print ('Uniq extraction is done.', t1 - t0)
    # print ("genome num:", len(uniq_segs_loci))
    # UHGG_snp(uniq_segs_loci)


    UHGG_cami()
