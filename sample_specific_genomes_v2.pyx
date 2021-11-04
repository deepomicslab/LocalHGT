#python3.7
import numpy as np
import time
cimport cython
from libc.string cimport memset
import array
from libc.stdlib cimport malloc, free
import itertools
import re
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
import logging
import dnaio


cdef extern from "Python.h":
    char* PyUnicode_AsUTF8(object unicode)

def DNA_complement2(sequence):
    # sequence = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

cdef class Index_Ref():
    cdef unsigned char [:] AG_ref
    cdef int kmer
    cdef int[:,:] coder_array
    cdef int[:,:] to_num
    cdef int[:] complementary
    cdef int[:] valid_base
    cdef long[:] base
    cdef long[:] single_array_size
    cdef int coder_num
    cdef long number
    cdef long single_kmer_table_size
    cdef long my_scaffold, my_ref_len
    cdef str index_dir
    def __cinit__(self, int k, int coder_num, str index_dir):
        #cdef long number = coder_num * 2 ** k 
        self.index_dir = index_dir
        self.number = coder_num * 2 ** k 
        self.single_kmer_table_size = 2 ** k 
        cdef int[:,:] permutation
        cdef int i, j, z
        cdef int[:] random_code
        self.coder_num = coder_num
        self.single_array_size = np.zeros(coder_num, dtype = 'int64')
        for i in range(coder_num):
            self.single_array_size[i] = i*2 ** k
        #self.single_array_size = np.array([0, 2 ** k, 2*2**k])
        self.my_scaffold = 0
        self.my_ref_len = 0
        self.complementary = np.zeros(85, dtype = 'int32')
        
        self.complementary[65] = 84         # A (65) - T(84) # C(67) - G(71)
        self.complementary[84] = 65
        self.complementary[67] = 71
        self.complementary[71] = 67
        self.kmer = k
        p = itertools.permutations([0, 1, 2])
        permutation = np.array([list(x) for x in p], dtype = 'int32') #np.random.randint(6, size = k, dtype = 'int32')#
        self.coder_array = np.zeros((k, coder_num), dtype = 'int32')
        
        for i in range(coder_num):
            if i % 3 == 0:
                random_code = np.random.randint(6, size = k, dtype = 'int32')
            z = i % 3
            for j in range(k):
                self.coder_array[j][i] = permutation[random_code[j]][z]
        np.save(self.index_dir +'/'+'coder', self.coder_array)

        self.to_num = np.zeros((90, 3), dtype = 'int32')

        self.to_num[84][1] = 1
        self.to_num[84][2] = 1
        self.to_num[67][0] = 1
        self.to_num[67][2] = 1
        self.to_num[71][0] = 1
        self.to_num[71][1] = 1

        self.valid_base = np.zeros(256, dtype = 'int32')
        self.valid_base[65] = 1
        self.valid_base[84] = 1
        self.valid_base[67] = 1
        self.valid_base[71] = 1

        self.base = np.zeros(k, dtype = 'int64')
        for z in range(k):
            self.base[z] = 2 ** (k-1-z)

    def filter_ref(self, reffile):
        t0=time.time()
        cdef int ref_len
        cdef int i, z
        cdef str line
        cdef int plasmid = 0        
        cdef int scanned_scaffold = 1
        f = open(self.index_dir +'/species_len.txt', 'w')
        cdef str scaffold_seq = ''
        for line in open(reffile, 'r'):
            line = line.strip()
            if line == '':
                continue
            if line[0] == '>':
                if scaffold_seq == '':
                    pass
                else:
                    ###for each scaffold
                    # if len(scaffold_seq) > 10000 and plasmid == 0:
                    self.choose_chr(scaffold_seq,scaffold_name,scanned_scaffold,f)
                    scanned_scaffold += 1
                scaffold_name = line[1:].split()[0]
                if scanned_scaffold % 1000 == 0:
                    print ('new scanned_scaffold',scanned_scaffold, time.time()-t0,scaffold_name)                
                scaffold_seq = ''
                continue
            scaffold_seq += line.upper()
        ###for each scaffold
        self.choose_chr(scaffold_seq, scaffold_name,scanned_scaffold,f)
        f.close()
        print ('ref index done.', time.time()-t0)
    def choose_chr(self, sequence, scaffold_name, scanned_scaffold, f):
        cdef long first_kmer_index, second_kmer_index, kmer_index, real_index
        cdef int i, j, z, n
        cdef int ref_len = len(sequence) 
        cdef unsigned char *c_ref
        cdef no_zero = 0
        cdef unsigned long[:,:] save_index = np.zeros((ref_len, self.coder_num), dtype = 'uint64')
        cdef bint N_flag

        c_ref = <unsigned char *> malloc((ref_len) * sizeof(unsigned char)) 
        for i in range(ref_len):
            c_ref[i] = ord(sequence[i]) 

        # print (ref_len)
        f.write(scaffold_name + '\t' + str(ref_len) + '\n') 
                
        for i in range(self.coder_num):
            for j in range(ref_len - self.kmer + 1):
                kmer_index = 0
                first_kmer_index = 0
                second_kmer_index = 0
                N_flag = False
                for z in range(self.kmer):
                    n = c_ref[j+z]
                    if self.valid_base[n] == 0: #not valid base
                        N_flag = True
                        break
                    first_kmer_index += self.base[z]*self.to_num[n][self.coder_array[z][i]]
                    n = self.complementary[n]
                    second_kmer_index += self.base[self.kmer-1-z] * self.to_num[n][self.coder_array[self.kmer-1-z][i]] 

                if first_kmer_index <= second_kmer_index:
                    kmer_index = first_kmer_index
                else:
                    kmer_index = second_kmer_index
            
                if N_flag: #kmer index is zero if the kmer contain N.
                    kmer_index = 0
                else:
                    real_index = kmer_index + self.single_array_size[i]
                    save_index[j][i] = real_index
        # new_depth_ref = [x for x in depth_ref[:ref_len]]
        # np.save(self.index_dir +'/'+scaffold_name, save_index)
        free(c_ref)

cdef class Extract_frag():
    cdef unsigned char *fastq_array
    cdef int reads_len, kmer, normal_kmer
    cdef int[:,:] coder_array
    cdef int[:,:] to_num
    cdef int[:] valid_base
    cdef int[:] complementary
    cdef long[:] base
    cdef long[:] single_array_size
    cdef int coder_num, window
    cdef int least_depth
    cdef str ID, index_dir
    cdef long number
    cdef long single_kmer_table_size
    cdef long my_scaffold, my_ref_len, fragment_num
    cdef float coverage_ratio, downsampling_ratio, coverage_k_ratio 

    def __cinit__(self, int coder_num, int reads_len, ID, coverage_ratio, coverage_k_ratio,\
     least_depth, window, downsampling_ratio, index_dir):
        self.index_dir = index_dir
        self.coder_array = np.load(self.index_dir + '/coder.npy')
        self.kmer = len(self.coder_array)
        print ('In reference index, kmer length is %s; the number of coder is %s.'%(self.kmer, len(self.coder_array[0])))    
        self.number = coder_num * 2 ** self.kmer
        self.single_kmer_table_size = 2 ** self.kmer
        cdef int[:,:] permutation
        cdef int i, j, z
        cdef int[:] random_code
        self.coder_num = coder_num
        self.window = window
        # allocate some memory (uninitialised, may contain arbitrary data)
        self.fastq_array = <unsigned char *> malloc(self.number * sizeof(unsigned char))
        memset(self.fastq_array, 0, self.number)
        self.single_array_size = np.zeros(coder_num, dtype = 'int64')
        for i in range(coder_num):
            self.single_array_size[i] = i*2 ** self.kmer
        self.coverage_ratio = coverage_ratio
        self.downsampling_ratio = downsampling_ratio 
        self.ID = ID
        self.least_depth = least_depth
        self.my_scaffold = 0
        self.my_ref_len = 0
        self.complementary = np.zeros(85, dtype = 'int32')
        self.complementary[65] = 84
        self.complementary[84] = 65
        self.complementary[67] = 71
        self.complementary[71] = 67
        self.reads_len = reads_len #150
        self.coverage_k_ratio = coverage_k_ratio
        self.fragment_num = 0


        self.to_num = np.zeros((90, 3), dtype = 'int32')
        self.to_num[84][1] = 1
        self.to_num[84][2] = 1
        self.to_num[67][0] = 1
        self.to_num[67][2] = 1
        self.to_num[71][0] = 1
        self.to_num[71][1] = 1

        self.valid_base = np.zeros(256, dtype = 'int32')
        self.valid_base[65] = 1
        self.valid_base[84] = 1
        self.valid_base[67] = 1
        self.valid_base[71] = 1

        self.base = np.zeros(self.kmer, dtype = 'int64')
        for z in range(self.kmer):
            self.base[z] = 2 ** (self.kmer-1-z)

    def reads2kmer(self, str reads):  
        cdef long first_kmer_index, second_kmer_index, kmer_index, real_index
        cdef int i, j, z, n
        cdef int[:] c_reads = np.empty(self.reads_len, dtype = 'int32')
        cdef bint N_flag


        for i in range(self.reads_len):
            c_reads[i] = ord(reads[i])
        for i in range(self.coder_num):            
            for j in range(self.reads_len - self.kmer + 1):               
                first_kmer_index = 0
                second_kmer_index = 0
                N_flag = False
                for z in range(self.kmer):
                    n = c_reads[j+z]
                    if self.valid_base[n] == 0: #not valid base
                        N_flag = True
                        break
                    first_kmer_index += self.base[z] * self.to_num[n][self.coder_array[z][i]]
                    n = self.complementary[n]
                    second_kmer_index += self.base[self.kmer-1-z] * self.to_num[n][self.coder_array[self.kmer-1-z][i]] 
                if first_kmer_index <= second_kmer_index:
                    kmer_index = first_kmer_index
                else:
                    kmer_index = second_kmer_index  
                if N_flag:
                    kmer_index = 0
                else:
                    real_index = kmer_index + self.single_array_size[i]
                    if self.fastq_array[real_index] < self.least_depth:
                        self.fastq_array[real_index] += 1

    def filter_ref(self, str fq1, str fq2, str reffile, str outdir):
        t0=time.time()
        cdef int ref_len
        cdef long i, z, reads_num
        cdef int scanned_scaffold = 1
        cdef str scaffold_name, line
        extracted_ref_interval_file = outdir + '%s.extract.ref.interval.txt'%(self.ID)
        # use dnaio module
        with dnaio.open(fq1) as f1, dnaio.open(fq2) as f2:   
            i = 0
            for record1 in f1:
                if np.random.random() < self.downsampling_ratio: 
                    self.reads2kmer(record1.sequence)
                    # break
                i += 1 
            for record1 in f2:
                if np.random.random() < self.downsampling_ratio: 
                    self.reads2kmer(record1.sequence)
                    # break
                i += 1 
                           
        reads_num = i/2
        f1.close()
        f2.close()
        logging.info("Reading fastq into table takes %s"%(time.time()-t0))
        print("Reading fastq into table takes %s"%(time.time()-t0))
        
        t1 = time.time()      
        f = open(extracted_ref_interval_file, "w")
        f.close()
        f = open(extracted_ref_interval_file, "a")
        for line in open(self.index_dir + '/species_len.txt', 'r'):
            array = line.strip().split()
            scaffold_name = array[0]
            # print ('scanned_scaffold',scanned_scaffold,scaffold_name) 
            self.choose_chr(scaffold_name, f, scanned_scaffold)
            if scanned_scaffold % 10000 == 0:
                print ('scanned_scaffold',scanned_scaffold,scaffold_name,self.my_ref_len, time.time()-t1) 
            scanned_scaffold += 1      
        f.close()
        free(self.fastq_array)
        t2 = time.time()
        print ('%s of the reference is extracted.'%(self.my_ref_len))
        logging.info("%s scaffolds were selected."%(self.my_scaffold))
        logging.info("%s fragments were selected."%(self.fragment_num))
        logging.info("Extracted referece size is %s bp."%(self.my_ref_len))
        logging.info("Sliding the reference costs %s s."%(t2-t1))  
        self.extract_seq(reffile, extracted_ref_interval_file,outdir)
        logging.info("Extracting sequence costs %s s."%(time.time() - t2))  
        return self.my_ref_len
    def choose_chr(self, str scaffold_name, f, int scanned_scaffold): 
        cdef long real_index
        cdef int i, j
        cdef unsigned long[:,:] save_index = np.load(self.index_dir + '/' +  scaffold_name+'.npy')
        cdef int ref_len = len(save_index)
        cdef unsigned char[:] depth_ref = np.zeros((ref_len + 1), dtype = 'uint8')
        cdef unsigned char[:] kmer_cover = np.zeros((ref_len + 1), dtype = 'uint8')
        cdef unsigned char perfect_depth
        cdef int start = 0
        cdef int end=0
        cdef int all_windows = 0
        cdef int kmer_cover_windows = 0
        cdef int min_cov = int(self.window * self.coverage_ratio)
        cdef int min_k_cov = int(self.window * self.coverage_k_ratio)
        cdef int frag_len = 0
        cdef int merge_close_fragment = 200
        cdef int[:] my_fragments = np.empty(20000, dtype = 'int32')
        cdef int frag_index = 0
        cdef bint flag = False
        cdef bint good_window = False
        cdef int windows_base_num = 0
        
        for j in range(ref_len - self.kmer + 1):
            perfect_depth = 0
            for i in range(self.coder_num):
                real_index = save_index[j][i]  
                if self.fastq_array[real_index] == self.least_depth:
                    if depth_ref[j] == 0:
                        depth_ref[j] = 1
                    perfect_depth += 1              
            if perfect_depth == self.coder_num:
                kmer_cover[j] = 1
      
        for j in range(ref_len - self.kmer + 1):
            if windows_base_num < self.window:
                all_windows += depth_ref[j]
                kmer_cover_windows += kmer_cover[j]
                windows_base_num += 1
            else:
                all_windows = all_windows - depth_ref[j-self.window] + depth_ref[j]
                kmer_cover_windows = kmer_cover_windows - kmer_cover[j-self.window] + kmer_cover[j]

            
            if all_windows >= min_cov and kmer_cover_windows > min_k_cov:
                good_window = True
            else:
                good_window = False

            if flag == False and good_window:
                start = j - 2 * self.window
                if start < 1:
                    start = 1
                flag = True
            if flag == True and good_window == False:
                end = j + self.window
                if frag_index > 0 and start - my_fragments[2*frag_index-1] < merge_close_fragment:
                    my_fragments[2*frag_index-1] = end
                else:
                    my_fragments[2 * frag_index] = start
                    my_fragments[2 * frag_index+1] = end
                    frag_index += 1
                flag = False

        if flag == True and good_window == True:
            end = ref_len
            if frag_index > 0 and start - my_fragments[2*frag_index-1] < merge_close_fragment:
                my_fragments[2*frag_index-1] = end
            else:
                my_fragments[2 * frag_index] = start
                my_fragments[2 * frag_index+1] = end
                frag_index += 1

        if frag_index > 0:
            self.my_scaffold += 1       
            for i in range(frag_index):
                self.fragment_num += 1
                self.my_ref_len += (my_fragments[2 * i+1] - my_fragments[2 * i])  
                f.write(scaffold_name + '\t' + str(my_fragments[2 * i]) + '\t' + str(my_fragments[2 * i+1]) + '\n')
                # print (scanned_scaffold, self.my_ref_len, frag_index, scaffold_name , str(my_fragments[2 * i]) , str(my_fragments[2 * i+1]))

    def extract_seq(self, reffile, extracted_ref_interval_file,outdir):
        # extracted_ref = reffile + '.%s.mid.fasta'%(self.ID)
        extracted_ref = outdir + '/%s.extract.ref.mid.fasta'%(self.ID)
        f = open(extracted_ref_interval_file, "r")
        locus_list = []
        for line in f:
            array = line.strip().split()
            locus_list.append('%s:%s-%s'%(array[0], array[1], array[2]))
        print ('find the extracted region, start samtools extraction.')
        samtools_extract_ref(reffile, locus_list, extracted_ref)
        f.close()

def samtools_extract_ref(ref, locus_list, extracted_ref):
    # ID = 'test'  
    os.system('samtools faidx %s'%(ref))
    max_extract = 300 
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

def extraction(fq1, fq2, coder_num, ID, reads_len, ref, coverage_ratio, coverage_k_ratio, \
    least_depth, window, downsampling_ratio, outdir, index_dir):  
    t0=time.time()
    # print (ID, 'start')
    ch = Extract_frag(coder_num, reads_len, ID, coverage_ratio, coverage_k_ratio, least_depth, window, downsampling_ratio, index_dir)
    my_ref_len = ch.filter_ref(fq1, fq2, ref, outdir)
    ref = ref + '.%s.mid.fasta'%(ID)
    print ('%s of the reference is extracted in %s.'%(my_ref_len, ID), time.time()-t0)
    return my_ref_len

def index(k, coder_num, ref, index_dir):
    os.system('mkdir %s'%(index_dir))
    ind = Index_Ref(k, coder_num, index_dir)
    ind.filter_ref(ref)










