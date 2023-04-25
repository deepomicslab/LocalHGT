
import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mechanism_taxonomy import Taxonomy
import pandas as pd

level_list = ["phylum", "class", "order", "family", "genus", "species", "genome"]

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

class Basic_count():

    def __init__(self):
        self.cohort_data = {}
        
    def read_samples(self):

        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        os.system(f"ls {tgs_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        os.system(f"ls {wenkui_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            if len(my_bkps) > 0:
                self.cohort_data[sra_id] = my_bkps
        print ("data loaded")
        
    def read_bkp(self, bkp_file):
        my_bkps = []
        f = open(bkp_file)
        all_rows = csv.reader(f)
        sample_dict = {}
        total_HGT_split_num = 0
        reads_num = 0
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
                    # os.system("rm %s"%(bkp_file))
                    break
                eb = Acc_Bkp(row)
                eb.abundance = eb.cross_split_reads/reads_num
                if eb.from_ref_genome == eb.to_ref_genome:
                    continue
                if eb.abundance < abun_cutoff:
                    continue
                # if self.get_genome_taxa(eb.from_ref_genome) == self.get_genome_taxa(eb.to_ref_genome): # classify intra-genus and inter-genus HGTs
                #     continue
                total_HGT_split_num += eb.cross_split_reads
                # if eb.hgt_tag not in self.all_hgt:
                sample_dict[eb.hgt_tag] = 1
                my_bkps.append(eb)
        for eb in my_bkps:
            eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
        f.close()
        return my_bkps


    def count(self):
        # bkp number distribution across all the samples
        # count mean and median value of a sample
        bkp_count_list = []
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            bkp_count_list.append(len(sample_bkp_list))
        sorted_bkp_count_list = sorted(bkp_count_list)
        print ("sample num is %s, mean bkp count is %s, median bkp count is %s, minimum bkp count is %s,\
         max bkp count is %s."%(len(self.cohort_data), np.mean(bkp_count_list), np.median(bkp_count_list),\
          sorted_bkp_count_list[0], sorted_bkp_count_list[-1]))
    
    def compare_intra_inter(self, level):
        # count the intra-taxa HGT frequency in each sample
        inter_freq = []
        intra_freq = []
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            inter_num, intra_num = 0, 0
            for bkp in sample_bkp_list:
                s1 = get_genome_taxa(bkp.from_ref_genome, level)
                s2 = get_genome_taxa(bkp.to_ref_genome, level)
                if s1[1:] == '__' or  s2[1:] == '__':
                    continue
                if s1 == s2 :
                    intra_num += 1
                else:
                    inter_num += 1
            inter_freq.append(inter_num/(intra_num + inter_num))
            intra_freq.append(intra_num/(intra_num + inter_num))
        U1, p = mannwhitneyu(inter_freq, intra_freq)
        print (level, p, np.mean(intra_freq), np.median(intra_freq))

    def sort_taxa_by_freq(self, level):
        # sort taxa by its mean frequency in each sample
        total_freq_dict = {}
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            sample_count = {}
            for bkp in sample_bkp_list:
                s1 = get_genome_taxa(bkp.from_ref_genome, level)
                s2 = get_genome_taxa(bkp.to_ref_genome, level)
                if s1[1:] == '__' or  s2[1:] == '__':
                    continue
                if s1 not in sample_count:
                    sample_count[s1] = 0
                if s2 not in sample_count:
                    sample_count[s2] = 0
                sample_count[s1] += 1
                sample_count[s2] += 1

            sample_freq = get_freq(sample_count)
            for taxa in sample_freq:
                if taxa not in total_freq_dict:
                    total_freq_dict[taxa] = []
                total_freq_dict[taxa].append(sample_freq[taxa])
        sample_num = len(self.cohort_data)
        sorted_dict = count_mean_freq(total_freq_dict, sample_num)
        return sorted_dict

    def count_inter_taxa_HGT(self):
        select_num = 5
        for level in range(1, 7):
            sorted_count = self.sort_taxa_by_freq(level)
            class_list = [x[0] for x in sorted_count][:select_num] + ['Others']     #list(count_dict.keys())
            index_dict = {}
            for i in range(len(class_list)):
                index_dict[class_list[i]] = i
            for co in sorted_count:
                if co[0] not in index_dict:
                    index_dict[co[0]] = len(class_list)-1

            head_name = [x[0][3:] for x in sorted_count][:select_num] + ['Others']
            data = []
            for i in range(len(class_list)):
                one_class = [0]*len(class_list)
                for sample in self.cohort_data:
                    for bkp in self.cohort_data[sample]:
                        from_tax = get_genome_taxa(bkp.from_ref_genome, level)
                        to_tax = get_genome_taxa(bkp.to_ref_genome, level)

                        if from_tax not in index_dict:
                            index_dict[from_tax] = len(class_list)-1
                        if to_tax not in index_dict:
                            index_dict[to_tax] = len(class_list)-1

                        if index_dict[from_tax] == i:
                            j = index_dict[to_tax]
                            if j < i:
                                one_class[j] += 1
                        if index_dict[to_tax] == i:
                            j = index_dict[from_tax]
                            if j < i:
                                one_class[j] += 1
                data.append(one_class)
            df = pd.DataFrame(data, columns = head_name)
            df.index = head_name
            df.to_csv('/mnt/d/R_script_files/inter_taxa_files/inter_taxa_%s.csv'%(level_list[level-1]), sep='\t')
            print (level_list[level-1], "done")       
            
def get_genome_taxa(genome, level):
    # g1 = get_pure_genome(genome)
    taxa = taxonomy.taxonomy_dict[genome].split(";")[level]
    return taxa

def get_freq(count_dict):
    total_count = sum(list(count_dict.values()))
    freq_dict = {}
    for taxa in count_dict:
        freq_dict[taxa] = count_dict[taxa]/total_count
    return freq_dict

def count_mean_freq(raw_dict, sample_num):
    mean_freq_dict = {}
    for taxa in raw_dict:
        if len(raw_dict[taxa])/sample_num < 0.1: # the taxa must exists in more than 10% of samples
            continue
        mean_freq_dict[taxa] = np.mean(raw_dict[taxa])
    sorted_dict = sorted(mean_freq_dict.items(), key=lambda item: item[1], reverse = True)
    print ("<<<<<<<<<<<<<<<<", len(sorted_dict))
    for i in range(10):
        print (i, sorted_dict[i][0], sorted_dict[i][1])
    return sorted_dict
    

if __name__ == "__main__":
    bin_size = 100
    abun_cutoff = 1e-7  #1e-7

    taxonomy = Taxonomy()
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    # hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"
    tgs_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    wenkui_dir = "/mnt/d/breakpoints/HGT/CRC/wenkui/"
    ba = Basic_count()
    ba.read_samples()
    # ba.count()
    ba.count_inter_taxa_HGT()
    # for level in range(1, 7):
    #     # ba.compare_intra_inter(level)
    #     # ba.sort_taxa_by_freq(level)
    #     ba.sort_inter_taxa_by_freq(level)
    