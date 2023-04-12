import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
 
def read_meta():
    
    sra_sample_dict = {}

    for line in open(meta_data):
        if line.strip() == '':
            continue
        array = line.strip().split(',')
        if array[0] != 'Run':
            sra_id = array[0]
            sample_id = array[-2]
            if re.search("_", sample_id):
                sample_id = sample_id.split("_")[1]
            sra_sample_dict[sra_id] = sample_id
    return sra_sample_dict

def read_design():
    
    sample_individual_dict = {}
    sample_time_point = {}

    for line in open(design_file):
        if line.strip() == '':
            continue
        array = line.strip().split()
        if array[0] != 'Sample':
            sample_id = array[0]
            individual = array[3]
            sample_individual_dict[sample_id] = individual
            sample_time_point[sample_id] = int(array[4])
    return sample_individual_dict, sample_time_point


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



class Time_Lines():

    def __init__(self):
        self.cohort_data = {}
        self.all_hgt = {}
        self.hgt_count = {} # record the hgt exists in how many samples
        self.sample_array_dict = {}
        self.HGT_list = []
        
    def read_samples(self):

        all_acc_file = result_dir + "/acc.list"
        os.system(f"ls {result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            
            td_id = sra_sample_dict[sra_id]
            if td_id[:2] == "CD":
                continue
            # print (sample_individual_dict[td_id])
            # if int(sample_individual_dict[td_id]) not in [1, 2, 3, 4,5, 6, 7, 8, 9, 10]:
            #     continue
            
            my_bkps = self.read_bkp(acc_file)
            self.cohort_data[sra_id] = my_bkps
        self.filter_hgt()
        
    def read_bkp(self, bkp_file):
        my_bkps = []
        f = open(bkp_file)
        all_rows = csv.reader(f)
        sample_dict = {}
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
                sample_dict[eb.hgt_tag] = 1
                my_bkps.append(eb)
        for eb in my_bkps:
            eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
        f.close()
        for hgt_tag in sample_dict:
            if hgt_tag not in self.hgt_count:
                self.hgt_count[hgt_tag] = 0
            self.hgt_count[hgt_tag] += 1
        return my_bkps

    def filter_hgt(self):

        for hgt_tag in self.hgt_count:
            if self.hgt_count[hgt_tag] >= sample_cutoff:
                hgt_index = len(self.all_hgt)
                self.all_hgt[hgt_tag] = hgt_index

    def get_HGT_table(self):
        for sra_id in self.cohort_data:
            self.sample_array_dict[sra_id] = [0]*len(self.all_hgt)
            for hgt in self.cohort_data[sra_id]:
                if hgt.hgt_tag not in self.all_hgt:
                    continue
                hgt_index = self.all_hgt[hgt.hgt_tag]
                self.sample_array_dict[sra_id][hgt_index] = 1 #hgt.abundance  # 1  #hgt.abundance
        self.get_HGT_list()

    def remove_both_zero(self, original_1, original_2):
        new_1 = []
        new_2 = []
        for i in range(len(original_1)):
            if original_1[i] != 0 or original_2[i] != 0:
                new_1.append(original_1[i])
                new_2.append(original_2[i])
        res = stats.spearmanr(new_1, new_2)
        return res.correlation

    def get_pearson(self):
        same_individual = []
        diff_individual = []
        diff_time_points = {}
        individual_diff_time_points = {}
        sample_list = list(self.sample_array_dict.keys())

        for i in range(len(sample_list)):
            sample1 = sample_list[i]
            ind_1 = sample_individual_dict[sra_sample_dict[sample1]]
            single_dict = {}
            for j in range(0, len(sample_list)):
                if i == j:
                    continue

                sample2 = sample_list[j]
                ind_2 = sample_individual_dict[sra_sample_dict[sample2]]
                # if ind_1 != ind_2:
                #     continue
                if ind_2 in single_dict and ind_2 != ind_1:
                    continue
                if ind_2 != ind_1:
                    single_dict[ind_2] = 1
                # print (sra_sample_dict[sample1], sra_sample_dict[sample2])
                # pearson = stats.pearsonr(self.sample_array_dict[sample1], self.sample_array_dict[sample2])
                # correlation = pearson[0]
                
                res = stats.spearmanr(self.sample_array_dict[sample1], self.sample_array_dict[sample2])
                correlation = res.correlation

                # correlation = self.remove_both_zero(self.sample_array_dict[sample1], self.sample_array_dict[sample2])
                
                if ind_1 == ind_2:
                    same_individual.append(correlation)
                    time_1 = sample_time_point[sra_sample_dict[sample1]]
                    time_2 = sample_time_point[sra_sample_dict[sample2]]
                    diff_time = abs(time_1 - time_2)
                    # print (sra_sample_dict[sample1], sra_sample_dict[sample2], diff_time)
                    if diff_time not in diff_time_points:
                        diff_time_points[diff_time] = []
                    diff_time_points[diff_time].append(correlation)

                    if ind_1 not in individual_diff_time_points:
                        individual_diff_time_points[ind_1] = [[],[],[],[],[],[],[],[],[]]
                    individual_diff_time_points[ind_1][diff_time-1] += [correlation]

                else:
                    diff_individual.append(correlation)
            # print (ind_1, single_dict, len(single_dict))
                # print (ind_1, ind_2, sample1, sample2, pearson[0], sum(self.sample_array_dict[sample1]), sum(self.sample_array_dict[sample2]))
        U1, p = mannwhitneyu(same_individual, diff_individual)
        print (p, np.mean(same_individual), np.mean(diff_individual), np.median(same_individual), np.median(diff_individual), len(same_individual), len(diff_individual))
        # print (np.mean(same_individual), np.median(same_individual), len(same_individual))

        # for diff_time in diff_time_points:
        #     print (diff_time, np.mean(diff_time_points[diff_time]), np.median(diff_time_points[diff_time]), len(diff_time_points[diff_time]))
        # for individual in individual_diff_time_points:
        #     print (individual) 
        #     for i in range(1, 10):
        #         ave = round(np.mean(individual_diff_time_points[individual][i-1]), 2)
        #         print (i, ave)

    def get_HGT_list(self):
        HGT_list = [0] * len(self.all_hgt)
        for hgt_tag in self.all_hgt:
            hgt_index = self.all_hgt[hgt_tag]
            HGT_list[hgt_index] = hgt_tag
        self.HGT_list = HGT_list

    def print_table(self):
        f = open(hgt_table, 'w')
        print ("SRA,Individual,Time,", end = '', file = f)
        for hgt_tag in self.HGT_list:
            print ("%s,"%(hgt_tag), end = '', file = f)
        print ('', file = f)
        for sra_id in self.sample_array_dict:
            individual = sample_individual_dict[sra_sample_dict[sra_id]]
            time = sample_time_point[sra_sample_dict[sra_id]]
            print ("%s,%s,%s,"%(sra_id, individual,time), end = '', file = f)
            for value in self.sample_array_dict[sra_id]:
                print ("%s,"%(value), end = '', file = f)
            print ('', file = f)
        f.close()
    
    def get_distribution(self):
        result = []
        for sra_id in self.sample_array_dict:
            result.append(self.sample_array_dict[sra_id])
        result = np.array(result)
        new_result = result.T
        data = []
        for new in new_result:
            data.append(sum(new))
        sns.displot(data, bins = 99)
        # Add labels and titles
        plt.xlabel('N.O. of Sample')
        plt.ylabel('Count')
        plt.savefig('/mnt/d/HGT/time_lines/distribution_plot.pdf')
        print (sorted(data, reverse=True)[:30])

    def count(self):
        individual_hgt_num = {}
        for sra_id in self.sample_array_dict:
            individual = sample_individual_dict[sra_sample_dict[sra_id]]
            if individual not in individual_hgt_num:
                individual_hgt_num[individual] = [0]*10
            individual_hgt_num[individual][sample_time_point[sra_sample_dict[sra_id]]-1] = sum(self.sample_array_dict[sra_id])
        for individual in individual_hgt_num:
            print ("individual:", individual, "HGT count at different time:", individual_hgt_num[individual])

    def find_conserve_HGT(self):
        individual_hgt_dict = {}
        for sra_id in self.sample_array_dict:
            individual = sample_individual_dict[sra_sample_dict[sra_id]]
            if individual not in individual_hgt_dict:
                individual_hgt_dict[individual] = {}
            for i in range(len(self.sample_array_dict[sra_id])):
                if self.sample_array_dict[sra_id][i] == 0:
                    continue
                hgt_tag = self.HGT_list[i]
                if hgt_tag not in individual_hgt_dict[individual]:
                    individual_hgt_dict[individual][hgt_tag] = 0
                individual_hgt_dict[individual][hgt_tag] += 1

        for individual in individual_hgt_dict:
            for hgt_tag in individual_hgt_dict[individual]:
                if individual_hgt_dict[individual][hgt_tag] == 10:
                    print ("individual:", individual, "conserved HGT:", hgt_tag)
        
if __name__ == "__main__":
    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    design_file = "/mnt/d/HGT/time_lines/sample_design.tsv"
    result_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    hgt_table = "/mnt/d/HGT/time_lines/SRP366030.HGT.table.csv"

    bin_size = 500
    split_cutoff = 0  #10
    sample_cutoff = 8  # 8
    abun_cutoff = 1e-7  #1e-7


    # for bin_size in [100, 200, 500, 1000]:
    #     for split_cutoff in range(6, 20, 2):
    #         for sample_cutoff in range(2, 11, 2):
    sra_sample_dict = read_meta()
    sample_individual_dict, sample_time_point = read_design()
    # print (sample_individual_dict)
    tim = Time_Lines()
    tim.read_samples()
    # print (tim.cohort_data.keys())
    tim.get_HGT_table()
    tim.get_pearson()
    print (len(tim.cohort_data), len(tim.all_hgt))
    print (bin_size, split_cutoff, sample_cutoff)
    # tim.print_table()
    # tim.get_distribution()
    # tim.find_conserve_HGT()