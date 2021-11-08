import sys
import os
from pyfaidx import Fasta

def extract_ref_seq(scaffold_name, start, end):
    return ref_fasta[scaffold_name][start:end].seq

def index2chr():
    convert_dict = {}
    f = open(reffile)
    index = 0
    for line in f:
        if line[0] == '>':
            chr_name = line[1:].split()[0]
            convert_dict[index] = chr_name
            index += 1
    return convert_dict

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

def extract_seq():   
    f = open(extracted_ref_interval_file, "r")
    locus_list = []
    for line in f:
        array = line.strip().split()
        locus_list.append('%s:%s-%s'%(array[0], array[1], array[2]))
    print ('start samtools extraction.')
    samtools_extract_ref(reffile, locus_list, extracted_ref)
    f.close()

def extract():
    genes = Fasta(reffile)
    outf = open(extracted_ref, 'w')
    f = open(extracted_ref_interval_file, "r")
    h = open(extracted_ref_interval_file+'.eva.txt', "w")
    for line in f:
        array = line.strip().split()
        # locus_list.append('%s:%s-%s'%(array[0], array[1], array[2]))
        chr_index=int(array[0])-1
        start = int(array[1])
        end = int(array[2])
        print (">%s:%s-%s"%(genes[chr_index].name, start, end), file = outf)
        print (genes[chr_index].name, start, end, file = h)
        print (genes[chr_index][start:end], file = outf)
    outf.close()
    f.close()
    h.close()

def find_chr_name():
    genes = Fasta(reffile)
    f = open(extracted_ref_interval_file, "r")
    h = open(extracted_ref_interval_file+'.eva.txt', "w")
    for line in f:
        array = line.strip().split()
        print (index2name_dict[int(array[0])], array[1], array[2], file = h)
    f.close()
    h.close()

def get_name():
    name_file = reffile + '.name'
    if not os.path.isfile(name_file):
        f = open(name_file, 'w')
        for line in open(reffile):
            if line[0] != ">":
                continue
            name = line[1:].split()[0]
            print (name, file = f)
        f.close()

def index2name():
    index2name_dict = {}
    name_file = reffile + '.name'
    ref_index = 1
    for line in open(name_file):
        index2name_dict[ref_index] = line.strip()
        ref_index += 1
    return index2name_dict

if __name__ == "__main__":
    reffile = sys.argv[1]
    extracted_ref = sys.argv[2]
    extracted_ref_interval_file = sys.argv[3]


    get_name()
    index2name_dict = index2name()
    find_chr_name()
    extracted_ref_interval_file = extracted_ref_interval_file+'.eva.txt'
    extract_seq()

