import sys
import time
import os


def find_chr_name():
    f = open(extracted_ref_interval_file, "r")
    h = open(extracted_ref_interval_file+'.bed', "w")
    
    for line in f:
        array = line.strip().split()
        print ('%s:%s-%s'%(index2name_dict[int(array[0])], array[1], array[2]), file = h)
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
    extracted_ref_interval_file = sys.argv[2]


    get_name()
    index2name_dict = index2name()
    find_chr_name()

