import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mechanism_taxonomy import Taxonomy
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from collections import defaultdict
 
def read_meta():
    
    sra_sample_dict = {}

    for line in open(meta_data):
        if line.strip() == '':
            continue
        array = line.strip().split(',')
        if array[0] != 'Run':
            if array[-6] != "ILLUMINA":
                continue
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
        for hgt_tag in sample_dict:
            if hgt_tag not in self.hgt_count:
                self.hgt_count[hgt_tag] = 0
            self.hgt_count[hgt_tag] += 1
        return my_bkps

    def get_genome_taxa(self, genome):
        # g1 = get_pure_genome(genome)
        taxa = taxonomy.taxonomy_dict[genome].split(";")[5]
        return taxa

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
            for j in range(i+1, len(sample_list)):
                if i == j:
                    continue

                sample2 = sample_list[j]
                ind_2 = sample_individual_dict[sra_sample_dict[sample2]]
                # if ind_1 != ind_2:
                #     continue
                # if ind_2 in single_dict and ind_2 != ind_1:
                #     continue
                # if ind_2 != ind_1:
                #     single_dict[ind_2] = 1
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

        data = []
        for value in same_individual:
            data.append([value, "same"])
        for value in diff_individual:
            data.append([value, "different"])
        df = pd.DataFrame(data, columns = ["correlation", "group"])
        df.to_csv("/mnt/d/R_script_files/time_line_correlation.csv", sep=',')

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

        record_conseve_species = {}
        converve_hgt_count = []
        conserve_hgt_tag = defaultdict(set)
        for individual in individual_hgt_dict:
            record_conseve_species[individual] = {}
            num = 0
            for hgt_tag in individual_hgt_dict[individual]:
                if individual_hgt_dict[individual][hgt_tag] == 10:
                    
                    num += 1
                    # print ("individual:", individual, "conserved HGT:", hgt_tag)
                    hgt_array = hgt_tag.split("&")
                    g1 = get_pure_genome(hgt_array[0])
                    g2 = get_pure_genome(hgt_array[2])
                    s1 = taxonomy.taxonomy_dict[g1].split(";")[6]
                    s2 = taxonomy.taxonomy_dict[g2].split(";")[6]
                    genus1 = taxonomy.taxonomy_dict[g1].split(";")[6]
                    genus2 = taxonomy.taxonomy_dict[g2].split(";")[6]
                    # print ("individual:", individual, s1, s2)
                    species_pair = "&".join(sorted([genus1, genus2]))
                    conserve_hgt_tag[species_pair].add(individual)

                    record_conseve_species[individual][s1] = 1
                    record_conseve_species[individual][s2] = 1
            converve_hgt_count.append(num)
            # print (list(record_conseve_species[individual].keys()))
        print ("conserve HGT count distribution", np.mean(converve_hgt_count), sorted(converve_hgt_count) )
        # get the conserve species that are common cross individuals
        species_dict = {}
        for individual in record_conseve_species:
            for species in record_conseve_species[individual]:
                if species not in species_dict:
                    species_dict[species] = 0
                species_dict[species] += 1
        sorted_dict = sorted(species_dict.items(), key=lambda item: item[1], reverse = True)
        conserve_hgt_tag_count = {}
        for species_pair in conserve_hgt_tag:
            conserve_hgt_tag_count[species_pair] = len(conserve_hgt_tag[species_pair])
        sorted_conserve_hgt_tag = sorted(conserve_hgt_tag_count.items(), key=lambda item: item[1], reverse = True)
        print ("sorted species_pair", sorted_conserve_hgt_tag[:5])

        common_species = 0
        uniq_species = 0
        for e in sorted_dict:
            # 
            if e[1] > 5:
                common_species += 1
                print (e[0], e[1])
            if e[1] == 1:
                uniq_species += 1
        species_num = len(species_dict) - 1
        print (species_num, common_species/species_num , uniq_species/species_num, uniq_species)
 
    def find_diff_HGT(self):
        # find the species with significant diff HGT 
        taxonomy = Taxonomy()
        ind_time_dict = {}

        for sra_id in self.sample_array_dict:
            individual = sample_individual_dict[sra_sample_dict[sra_id]]
            time_point = sample_time_point[sra_sample_dict[sra_id]]

            if individual not in ind_time_dict:
                ind_time_dict[individual] = {}
            if not (time_point == 1 or time_point == 10):
                continue
            ind_time_dict[individual][time_point] = self.sample_array_dict[sra_id]

        diff_dict = {} # record the genus pair that show diff in at least one sample
        data = []
        for individual in ind_time_dict:
            diff = np.array(ind_time_dict[individual][1]) - np.array(ind_time_dict[individual][10])
            species_array = {}
            for i in range(len(diff)):
                hgt_tag = self.HGT_list[i]
                hgt_array = hgt_tag.split("&")
                g1 = get_pure_genome(hgt_array[0])
                g2 = get_pure_genome(hgt_array[2])
                s1 = taxonomy.taxonomy_dict[g1].split(";")[5]
                s2 = taxonomy.taxonomy_dict[g2].split(";")[5]
                species_tag = "&".join(sorted([s1, s2]))

                if len(s1) == 3 or  len(s2) == 3:
                    continue

                if ind_time_dict[individual][1][i] == 1 or ind_time_dict[individual][10][i] == 1:
                    if species_tag not in species_array:
                        species_array[species_tag] = [[], []]
                    species_array[species_tag][0].append(ind_time_dict[individual][1][i])
                    species_array[species_tag][1].append(ind_time_dict[individual][10][i])

            # p_species = {}
            for species in species_array:
                U1, p = mannwhitneyu(species_array[species][0], species_array[species][1])
                # p_species[species] = p
                data.append([species, p, individual])
                # if p < 0.01:
                #     if species not in diff_dict:
                #         diff_dict[species] = []
                #     diff_dict[species].append(p)
            # sorted_dict = sorted(p_species.items(), key=lambda item: item[1])
            # print (list(diff_species.keys()))
            # print (individual, sorted_dict[:10])

        df = pd.DataFrame(data, columns = ["species", "p_value", "individual"])
        # print (df)
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected

        for index, row in df.iterrows():
            if row["p.adj"]  < 0.05:
                species = row["species"]
                if species not in diff_dict:
                    diff_dict[species] = []
                diff_dict[species].append(row["p.adj"])
        
        for genus_pair in diff_dict:
            if len(diff_dict[genus_pair]) < 1:
                continue
            print (genus_pair, len(diff_dict[genus_pair]), np.mean(diff_dict[genus_pair]), max(diff_dict[genus_pair]))
        print (len(diff_dict))

    def intra_inter_genus_HGT(self, level):
        # find the difference between intra and inter-taxa HGTs 
        # not significant
        taxonomy = Taxonomy()
        ind_time_dict = {}

        for sra_id in self.sample_array_dict:
            individual = sample_individual_dict[sra_sample_dict[sra_id]]
            time_point = sample_time_point[sra_sample_dict[sra_id]]

            if individual not in ind_time_dict:
                ind_time_dict[individual] = {}
            if not (time_point == 1 or time_point == 10):
                continue
            ind_time_dict[individual][time_point] = self.sample_array_dict[sra_id]

        inter_list = []
        intra_list = []
        
        for individual in ind_time_dict:
            diff = np.array(ind_time_dict[individual][1]) - np.array(ind_time_dict[individual][10])
            inter, intra  = 0, 0
            all_inter, all_intra = 0, 0
            species_dict = {}
            for i in range(len(diff)):
                hgt_tag = self.HGT_list[i]
                hgt_array = hgt_tag.split("&")
                g1 = get_pure_genome(hgt_array[0])
                g2 = get_pure_genome(hgt_array[2])
                s1 = taxonomy.taxonomy_dict[g1].split(";")[level]
                s2 = taxonomy.taxonomy_dict[g2].split(";")[level]
                species_tag = "&".join(sorted([s1, s2]))

                if s1[1:] == '__' or  s2[1:] == '__':
                    continue
                if s1 not in species_dict:
                    species_dict[s1] = [0, 0, 0, 0] # all intra, all inter, intra, inter
                if s2 not in species_dict:
                    species_dict[s2] = [0, 0, 0, 0]

                if abs(ind_time_dict[individual][1][i] - ind_time_dict[individual][10][i]) == 1:
                    if s1 != s2:
                        species_dict[s1][3] += 1
                        species_dict[s2][3] += 1
                        inter += 1
                    else:
                        species_dict[s1][2] += 1
                        species_dict[s2][2] += 1
                        intra += 1

                if ind_time_dict[individual][1][i] == 1 or ind_time_dict[individual][10][i] == 1:
                    if s1 != s2:
                        species_dict[s1][1] += 1
                        species_dict[s2][1] += 1
                        all_inter += 1
                    else:
                        species_dict[s1][0] += 1
                        species_dict[s2][0] += 1
                        all_intra += 1

            for s in species_dict:
                if species_dict[s][1] > 0:
                    inter_list.append(species_dict[s][3]/species_dict[s][1])
                if species_dict[s][0] > 0:
                    intra_list.append(species_dict[s][2]/species_dict[s][0])
            # inter_list.append(inter/all_inter)
            # intra_list.append(intra/all_intra)
            # intra_list.append()


        U1, inter_p = mannwhitneyu(inter_list, intra_list)
        print (level, inter_p, np.median(inter_list), np.median(intra_list))

    def comprare_intra_genus(self):
        # check whether intra-genus HGTs showed higher variabiligy than inter-genus HGTs
        # find the species with significant diff HGT 
        taxonomy = Taxonomy()
        ind_time_dict = {}

        for sra_id in self.sample_array_dict:
            individual = sample_individual_dict[sra_sample_dict[sra_id]]
            time_point = sample_time_point[sra_sample_dict[sra_id]]

            if individual not in ind_time_dict:
                ind_time_dict[individual] = {}
            if not (time_point == 1 or time_point == 10):
                continue
            ind_time_dict[individual][time_point] = self.sample_array_dict[sra_id]

        diff_dict = {} # record the genus pair that show diff in at least one sample
        for individual in ind_time_dict:
            diff = np.array(ind_time_dict[individual][1]) - np.array(ind_time_dict[individual][10])
            species_array = {}
            genus_HGT_types = [[], []]
            for i in range(len(diff)):
                hgt_tag = self.HGT_list[i]
                hgt_array = hgt_tag.split("&")
                g1 = get_pure_genome(hgt_array[0])
                g2 = get_pure_genome(hgt_array[2])
                s1 = taxonomy.taxonomy_dict[g1].split(";")[5]
                s2 = taxonomy.taxonomy_dict[g2].split(";")[5]

                if s1 == 'g__' or  s2 == 'g__':
                    continue

                if (ind_time_dict[individual][1][i] == 1 or ind_time_dict[individual][10][i] == 1) :

                    if s1 != s2:
                        genus_HGT_types[0].append(ind_time_dict[individual][1][i])
                        genus_HGT_types[1].append(ind_time_dict[individual][10][i])


            U1, p = mannwhitneyu(genus_HGT_types[0], genus_HGT_types[1])
            print (individual, p)


class Event(object):

    def __init__(self, array):
        self.sample = array[1]
        self.ins_genome = array[2]
        self.ins_pos = int(array[3])
        self.del_genome = array[4]
        self.del_start = int(array[5])
        self.del_end = int(array[6])
        self.reverse_flag = array[7]

class Transfer():

    def __init__(self, TD_sample_list):
        self.HGT_event_dict = {}
        self.max_diff = 50
        self.TD_sample_list = TD_sample_list

    def read_events(self, identified_hgt):
        for line in open(identified_hgt):
            array = line.strip().split(",")
            if array[1] == "sample":
                continue
            event = Event(array)

            ###### classify intra and inter-genera HGT events
            # genus_ins = taxonomy.taxonomy_dict["_".join(event.ins_genome.split("_")[:-1])].split(";")[5]
            # genus_del = taxonomy.taxonomy_dict["_".join(event.del_genome.split("_")[:-1])].split(";")[5]
            # # print (genus_ins, genus_del)
            # if genus_ins == "g__" or genus_del == "g__":
            #     continue
            # if genus_ins == genus_del:
            #     continue

            sample = array[1]
            if sample not in self.HGT_event_dict:
                self.HGT_event_dict[sample] = []
            self.HGT_event_dict[sample].append(event)
    
    def get_jaccard_dist(self, sample1, sample2):
        share_num = 0
        total_num = len(self.HGT_event_dict[sample1])
        for event2 in self.HGT_event_dict[sample2]:
            share_flag = False
            for event1 in self.HGT_event_dict[sample1]:
                if event1.ins_genome == event2.ins_genome and event1.del_genome == event2.del_genome and abs(event1.ins_pos - event2.ins_pos) < self.max_diff \
                    and abs(event1.del_start - event2.del_start) < self.max_diff and abs(event1.del_end - event2.del_end) < self.max_diff \
                    and event1.reverse_flag == event2.reverse_flag:
                    share_flag = True
                    break
            if share_flag:
                share_num += 1
            else:
                total_num += 1
        return share_num/total_num ## jacard distance

    def compare(self):
        same_individual = []
        diff_individual = []
        diff_time_points = {}
        individual_diff_time_points = {}
        sample_list = self.TD_sample_list

        for i in range(len(sample_list)):
            sample1 = sample_list[i]
            ind_1 = sample_individual_dict[sra_sample_dict[sample1]]
            single_dict = {}
            for j in range(i+1, len(sample_list)):
                if i == j:
                    continue
                sample2 = sample_list[j]
                ind_2 = sample_individual_dict[sra_sample_dict[sample2]]

                # if ind_2 in single_dict and ind_2 != ind_1:
                #     continue
                # if ind_2 != ind_1:
                #     single_dict[ind_2] = 1

                if sample1 not in self.HGT_event_dict or sample2 not in self.HGT_event_dict:
                    continue 
                # print (sample1, sample2)

                jac_dist = self.get_jaccard_dist(sample1, sample2)
                if ind_1 == ind_2:
                    same_individual.append(jac_dist)
                    time_1 = sample_time_point[sra_sample_dict[sample1]]
                    time_2 = sample_time_point[sra_sample_dict[sample2]]
                    diff_time = abs(time_1 - time_2)
                    # print (sra_sample_dict[sample1], sra_sample_dict[sample2], diff_time)
                    if diff_time not in diff_time_points:
                        diff_time_points[diff_time] = []
                    diff_time_points[diff_time].append(jac_dist)
                    if ind_1 not in individual_diff_time_points:
                        individual_diff_time_points[ind_1] = [[],[],[],[],[],[],[],[],[]]
                    individual_diff_time_points[ind_1][diff_time-1] += [jac_dist]
                else:
                    diff_individual.append(jac_dist)
        U1, p = mannwhitneyu(same_individual, diff_individual)
        print (p, np.mean(same_individual), np.mean(diff_individual), np.median(same_individual), np.median(diff_individual), len(same_individual), len(diff_individual))
        data = []
        for value in same_individual:
            data.append([value, "same"])
        for value in diff_individual:
            data.append([value, "different"])
        df = pd.DataFrame(data, columns = ["jaccard", "group"])
        df.to_csv("/mnt/d/R_script_files/time_line_jaccard.csv", sep=',')

class Country(Transfer):

    def __init__(self, TD_sample_list):
        super(Transfer, self).__init__(TD_sample_list)
    
def read_phenotype():
    phenotype_dict = {}
    pheno_result = "/mnt/d/HGT/association/phenotype.csv"
    for line in open(pheno_result):
        array = line.strip().split(",")
        ID = array[1]
        if ID == "sample":
            continue
        pheno = array[2:5]
        phenotype_dict[ID] = pheno
        cohort = array[2]
    return phenotype_dict

def get_samples_for_single():
    individual_dict = {}
    for name in sample_individual_dict:
        individual  = sample_individual_dict[name]
        if individual not in individual_dict:
            individual_dict[individual] = []
        for sra in sra_sample_dict:
            if sra_sample_dict[sra] ==  name:
                individual_dict[individual].append(sra)
    TD_sample_list = []
    for individual in individual_dict:
        # print (individual, individual_dict[individual])
        if len(individual_dict[individual]) > 1:
            TD_sample_list += individual_dict[individual]
    return TD_sample_list

def get_pure_genome(del_genome):
    return "_".join(del_genome.split("_")[:-1])

if __name__ == "__main__":
    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    design_file = "/mnt/d/HGT/time_lines/sample_design.tsv"
    result_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    hgt_table = "/mnt/d/HGT/time_lines/SRP366030.HGT.table.csv"
    identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"

    all_identified_hgt = "/mnt/d/HGT/seq_ana/bk/identified_event.csv"

    bin_size = 500
    sample_cutoff = 8  # 8
    abun_cutoff = 1e-7  #1e-7

    taxonomy = Taxonomy()
    # for bin_size in [100, 200, 500, 1000]:
    #     for split_cutoff in range(6, 20, 2):
    #         for sample_cutoff in range(2, 11, 2):
    sra_sample_dict = read_meta()
    sample_individual_dict, sample_time_point = read_design()




    #################  HGT event jaccard between samples
    # TD_sample_list = get_samples_for_single()
    # transfer = Transfer(TD_sample_list)
    # transfer.read_events(identified_hgt)
    # transfer.compare()





    ############### breakpoint correlation between samples
    tim = Time_Lines()
    tim.read_samples()
    tim.get_HGT_table()
    # tim.get_pearson()
    # print (len(tim.cohort_data), len(tim.all_hgt))
    # print (bin_size, sample_cutoff)





    # tim.print_table()
    # tim.get_distribution()
    # tim.find_conserve_HGT()
    tim.find_diff_HGT()
    # tim.comprare_intra_genus()
    # for level in range(1, 6):
    #     tim.intra_inter_genus_HGT(level)