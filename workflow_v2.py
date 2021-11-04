from sample_specific_genomes_v2 import extraction, index
import time
import numpy as np
import precise_map_by_bwa
import os
import linecache
import sys
import logging
from datetime import datetime


########for evaluation########
def read_interval(interval_file, true_locus):
    gap = 50
    cover_flag = False
    f = open(interval_file)
    for line in f:
        array = line.strip().split()
        if array[0] == true_locus[0] and int(true_locus[1]) > int(array[1]) + gap and int(true_locus[1]) < int(array[2]) - gap:
            cover_flag = True
            # print ('#', true_locus, array)
    return cover_flag

def check_if_bkp_in_extracted_ref(true, interval_file):
    true_bkp = read_true_loci(true)
    i = 0
    for true_locus in true_bkp:
        if read_interval(interval_file, true_locus):
            i += 1
        # else:
        #     print (true_locus)
    print (len(true_bkp), float(i)/len(true_bkp))
    return float(i)/len(true_bkp)

def read_true_loci(true):
        true_bkp = []
        for line in open(true):
                array = line.strip().split()
                true_bkp.append([array[0], array[1]])
                true_bkp.append([array[2], array[3]])
                true_bkp.append([array[2], array[4]])
        return true_bkp
########for evaluation########


def blast_map2_self(file, map_ref):
    makedb = 'makeblastdb -in %s -dbtype nucl -out ref'%(file)
    blastn = 'blastn -query %s -outfmt 6 -out %s -db ref'%(file, map_ref)
    os.system(makedb)
    os.system(blastn)

def count_read_length(fq1):
    i = 0
    for line in open(fq1):
        if i == 1:
            return len(line.strip())
        i += 1

def iter_count(file_name):
    from itertools import (takewhile, repeat)
    buffer = 1024 * 1024
    with open(file_name) as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        return sum(buf.count('\n') for buf in buf_gen)/4

def map_and_call(extracted_ref, ID, downsampling_ratio, fq1, fq2, output_file, outdir, bwa_thread, original_ref):
    t0 = time.time()
    os.system('bwa index %s'%(extracted_ref))
    t1 = time.time()
    order = """
    sample=%s
    ref=%s
    fq1=%s
    fq2=%s
    dir=%s
    # bwa index $ref
    samtools faidx $ref
    echo sample specific reference is indexed. 

    bwa mem -M -t %s -R "@RG\\tID:id\\tSM:sample\\tLB:lib" $ref $fq1 $fq2 > $dir/$sample.sam
    python3 %s/extractSplitReads_BwaMem.py -i $dir/$sample.sam > $dir/$sample.split.sam

    # Sort bam file
    #samtools sort -o $dir/$sample.bam $dir/$sample.unsort.bam
    #rm $dir/$sample.unsort.bam
    #samtools index $dir/$sample.bam   
    """%(ID, extracted_ref, fq1, fq2, outdir, bwa_thread, sys.path[0])
    os.system(order)
    t2 = time.time()
    logging.info("Alignment costs %s s; Index cost %s s."%(t2-t0, t1-t0))
    bam = '%s/%s.sam'%(outdir, ID)   
    split_bam = '%s/%s.split.sam'%(outdir, ID)  
    precise_map_by_bwa.precise_bp(bam, split_bam, output_file, original_ref)
    os.system('rm %s*'%(bam))
    t3 = time.time()
    logging.info("Finding bkp from split reads costs %s s."%(t3-t2))

def main(options):
    original_ref = options.ref
    ID = options.sample_ID
    fq1 = options.fq1
    fq2 = options.fq2
    outdir = options.outdir
    max_size = options.sampling_size*1000000
    k = options.kmer_size
    coverage_ratio = options.coverage_ratio
    coder_num = options.coder_num
    least_depth = options.least_depth
    window = options.window

    if not os.path.isfile(fq1):       
        if os.path.isfile(fq1+'.gz'):
            logging.warning("The fastq file not unzipped! Unzipping...")
            os.system('gzip -d %s.gz'%(fq1))
            os.system('gzip -d %s.gz'%(fq2))
        else:
            logging.error("The fastq file not found, check if unzipped.")
    t0 = time.time()
    reads_len = count_read_length(fq1)
    pair_reads_num = iter_count(fq1)
    datasize = round(pair_reads_num * reads_len * 2, 2)   
    downsampling_ratio = round(float(max_size)/datasize, 6)
    logging.info("Sample ID is %s."%(ID))
    logging.info("Data size is %s bp."%(datasize))
    logging.info("Sampling size is %s bp."%(max_size))
    logging.info("Sampling ratio is %s."%(downsampling_ratio))

    extracted_ref = outdir + '/%s.extract.ref.mid.fasta'%(ID)
    # if not os.path.isfile(extracted_ref):
    my_ref_len = extraction(fq1, fq2, coder_num, ID, reads_len, original_ref, coverage_ratio,\
    options.coverage_k_ratio, least_depth, window, downsampling_ratio, outdir, options.ref_index_dir)  
    logging.info("Extracted referece size is %s bp."%(my_ref_len))
    t1 = time.time()
    logging.info("Reference Extraction costs %s s."%(t1-t0))  
    # else:
    #     print ('The reference has been extracted already!')
    t1 = time.time()
    os.system('bash /mnt/d/breakpoints/script/pipeline.sh %s %s %s %s %s %s'%(extracted_ref, fq1, fq2, ID, outdir, original_ref))
    logging.info("Calling break points costs %s s."%(time.time()-t1))  
    logging.info("HGT Detection Finished.")
    logging.info("The HGT detection process costs %s s."%(time.time()-t0))

# print (ID, 'used time:', time.time() - t0)
# x = np.load('coder.npy')
# print (x)
if __name__ == "__main__":
    # ID = sys.argv[1]
    # print (ID, 'start call HGT...')
    # fq1=sys.argv[2]
    # fq2=sys.argv[3]
    # ref = sys.argv[4]  #'/mnt/d/breakpoints/big/gut.reference.filter.fa' 
    # outdir = sys.argv[5]
    # main(ref, fq1, fq2, ID, outdir)
    from argparse import ArgumentParser

    parser = ArgumentParser(description='HGT detection')
    parser.add_argument('-i', '--sample_ID', help='Sample ID', required=True)
    parser.add_argument('-1', '--fq1', help='unzipped fastq 1 file.', required=True)
    parser.add_argument('-2', '--fq2', help='unzipped fastq 2 file.', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', required=True)   
    parser.add_argument('-r', '--ref', help='ref', required=True)
    parser.add_argument('-x', '--ref_index_dir', help='ref index dir', required=True)
    #parameters
    parser.add_argument('-k', '--kmer_size', help='kmer_size, must same with that used in indexing ref.', required=False, default = 30, type=int)  
    parser.add_argument('-c', '--coverage_ratio', help='coverage_ratio', required=False, default = 0.3, type=float)
    parser.add_argument('-a', '--coverage_k_ratio', help='coverage_k_ratio', required=False, default = 0.01, type=float)
    parser.add_argument('-e', '--coder_num', help='coder_num', required=False, default = 3, type=int)
    parser.add_argument('-d', '--least_depth', help='least_depth', required=False, default = 1, type=int)
    parser.add_argument('-w', '--window', help='window size', required=False, default = 500, type=int)  
    parser.add_argument('-t', '--bwa_thread', help='bwa threads', required=False, default = 2, type=int)  
    parser.add_argument('-s', '--sampling_size', help='sampling_size (M). You can choose a very high value to avoid downsampling.', required=False, default = 500, type=int)
    options = parser.parse_args()
    give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
    logging.basicConfig(filename = options.outdir +'/'+ options.sample_ID +'_' + give_time + '.log',format='[%(asctime)s-%(filename)s-%(levelname)s:%(message)s]', level = logging.DEBUG,filemode='w')
    print (options.sample_ID, 'start call HGT...')
    logging.info("Running parameters: %s"%(' '.join(sys.argv)))
    logging.info("Start call HGT...")
    main(options)
    true = '/mnt/d/breakpoints/HGT/snp_fq/%s.HGT.true.sv.txt'%('_'.join(options.sample_ID.split('_')[:-1]))
    interval_file = options.outdir + '/%s.extract.ref.interval.txt'%(options.sample_ID)
    ref_accuracy = check_if_bkp_in_extracted_ref(true, interval_file)
    logging.info("Reference extraction accuracy is: %s"%(ref_accuracy))




    """
    ref_index_dir = '/mnt/e/HGT/ref_index_10_12'
    #index(32, 3, '/mnt/e/HGT/test.fa', ref_index_dir)
    index(32, 3, '/mnt/d/breakpoints/big/gut.reference.filter.fa', ref_index_dir)
    """

    """
    arser = ArgumentParser(description='HGT detection')
    parser.add_argument("function", 
                        nargs="?",
                        choices=['ind', 'run'],
                        default='function1',
                        help='Choose sub-script, set ind to index the reference, set run to detect HGT.',
                        )
    args, sub_args = parser.parse_known_args()

    if args.function == "ind":
        parser = ArgumentParser()
        parser.add_argument('-k', '--kmer_size', help='kmer_size, must same with that used in indexing ref.', required=False, default = 30, type=int)
        parser.add_argument('-e', '--coder_num', help='coder_num', required=False, default = 3, type=int)
        parser.add_argument('-r', '--ref', help='ref', required=True)
        options = parser.parse_args()
        index(options.kmer_size, options.coder_num, options.ref)
    elif args.function == "run":
        parser.add_argument('-i', '--sample_ID', help='Sample ID', required=True)
        parser.add_argument('-1', '--fq1', help='unzipped fastq 1 file.', required=True)
        parser.add_argument('-2', '--fq2', help='unzipped fastq 2 file.', required=True)
        parser.add_argument('-o', '--outdir', help='outdir', required=True)   
        parser.add_argument('-r', '--ref', help='ref', required=True)
        #parameters
        parser.add_argument('-k', '--kmer_size', help='kmer_size, must same with that used in indexing ref.', required=False, default = 30, type=int)  
        parser.add_argument('-c', '--coverage_ratio', help='coverage_ratio', required=False, default = 0.7, type=float)
        parser.add_argument('-e', '--coder_num', help='coder_num', required=False, default = 3, type=int)
        parser.add_argument('-d', '--least_depth', help='least_depth', required=False, default = 1, type=int)
        parser.add_argument('-w', '--window', help='window size', required=False, default = 500, type=int)  
        parser.add_argument('-t', '--bwa_thread', help='bwa threads', required=False, default = 2, type=int)  
        parser.add_argument('-s', '--sampling_size', help='sampling_size (M). You can choose a very high value to avoid downsampling.', required=False, default = 500, type=int)
        options = parser.parse_args()
        logging.basicConfig(filename = options.outdir +'/'+ options.sample_ID + '.log',format='[%(asctime)s-%(filename)s-%(levelname)s:%(message)s]', level = logging.DEBUG,filemode='w')
        print (options.sample_ID, 'start call HGT...')
        logging.info("Start call HGT...")
        main(options)  
    else:
        print ('Please select one sub-script, ind or run, set ind to index the reference, set run to detect HGT.')
    """

