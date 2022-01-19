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
import pandas as pd
import os
import pickle
from scipy.stats import mannwhitneyu
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums
import networkx as nx
import math

class Taxonomy():
    def __init__(self):
        self.taxonomy_dict = {}
        self.get()
        

    def read_UHGG(self):
        df = pd.read_table(UHGG_meta) 
        for index, row in df.iterrows():
            # self.read_UHGG[row["Genome"]] = row["Species_rep"]
            genome = row["Genome"]
            lineage = row["Lineage"]
            self.taxonomy_dict[genome] = lineage
        # for line in open(UHGG_meta):
        #     array = line.split()

    def get(self):
        save_file = "taxonomy_dict.pkl"
        if not os.path.isfile(save_file):
            self.read_UHGG()
            with open(save_file, 'wb') as f:
                pickle.dump(self.taxonomy_dict, f)
        else:
            with open(save_file, 'rb') as f:
                self.taxonomy_dict = pickle.load(f)

class Acc_Bkp(object):
    def __init__(self, list):
        self.from_ref = list[0]
        self.to_ref = list[2]
        self.from_bkp = list[1]
        self.to_bkp = list[3]
        self.if_reverse = list[6]
        self.from_side = list[4]
        self.to_side = list[5]
        self.from_ref_genome = "_".join(self.from_ref.split("_")[:-1])
        self.from_ref_lineage = taxonomy.taxonomy_dict[self.from_ref_genome]
        self.to_ref_genome = "_".join(self.to_ref.split("_")[:-1])
        self.to_ref_lineage = taxonomy.taxonomy_dict[self.to_ref_genome]
        # print (self.from_ref_genome, self.from_ref_lineage, self.to_ref_genome, self.to_ref_lineage)

    def print_out(self):
        print (self.from_ref, self.from_bkp, self.to_ref, self.to_bkp, self.from_side,\
         self.to_side, self.if_reverse, self.similarity)

    def write_out(self, writer):
        writer.writerow ([self.from_ref, self.from_bkp, self.to_ref, self.to_bkp, \
        self.from_side, self.to_side, self.if_reverse, self.read_str, self.ref_str, self.similarity])

class Phenotype():
    def __init__(self):
        self.name_disease = {}
        self.name_cohort = {}
        self.ID_disease = {}
        self.ID_cohort = {}
        self.ID_bases = {}
        self.read_pheno()
        self.read_sra_meta()
        

    def read_sra_meta(self):
        f = open(sra_meta)
        all_rows = csv.reader(f)
        for row in all_rows:
            if row[0] == "Run":
                continue   
            sra_ID = row[0]
            bases = int(row[3])
            sample_name = row[-2]
            # self.sra2name[sra_ID] = sample_name
            self.ID_disease[sra_ID] = self.name_disease[sample_name]
            self.ID_cohort[sra_ID] = self.name_cohort[sample_name]
            self.ID_bases[sra_ID] = bases
        f.close()

    def read_pheno(self):     
        df = pd.read_excel(pheno_file, header=None) 
        # print (df)
        for index, row in df.iterrows():
            sample_name = row[3]
            condition = row[6]
            full_disease = row[7]
            cohort = row[1]
            # print (cohort)
            self.name_disease[sample_name] = condition
            self.name_cohort[sample_name] = cohort
            # print (sample_name,condition) 

class Sample():

    def __init__(self, bkp_file, ID):
        self.bkp_file = bkp_file
        self.bkps = []
        self.filter_bkps = [] #recheck, not the same species
        self.ID = ID
        self.disease = phenotype.ID_disease[ID]
        self.cohort = phenotype.ID_cohort[ID]    
        self.bases = phenotype.ID_bases[ID]    
        self.single = {}
        self.pair = {}
        self.level = 3
        

        self.read_bkp()
        self.bkp_num = len(self.filter_bkps)
        self.average_bkp_per_species = round(self.bkp_num/len(self.single), 2)
        self.matrix = []
        self.nodes = list(self.single.keys()) #node num significant 4
        self.build_matrix()
        self.graph = nx.from_numpy_matrix(self.matrix)
        self.graph_density = nx.density(self.graph) #significant 4 5 6
        self.graph_average_clustering = nx.average_clustering(self.graph, nodes=None, weight=None, count_zeros=True) #not
        self.transtivity = nx.transitivity(self.graph) #significant 6
        self.weighted_degree = self.matrix.sum(axis=1, dtype='float') #list
        # print (self.weighted_degree)

    def read_bkp(self):
        f = open(self.bkp_file)
        all_rows = csv.reader(f)
        for row in all_rows:
            if row[0] == "from_ref":
                continue
            # print (row)
            eb = Acc_Bkp(row)
            self.bkps.append(eb)
            if eb.from_ref_genome != eb.to_ref_genome and eb.from_ref_lineage.split(";")[6] != eb.to_ref_lineage.split(";")[6]:
                self.filter_bkps.append(eb)

                from_tax = eb.from_ref_lineage.split(";")[self.level]
                to_tax = eb.to_ref_lineage.split(";")[self.level]

                if from_tax not in self.single:
                    self.single[from_tax] = 1
                else:
                    self.single[from_tax] += 1
                if to_tax not in self.single:
                    self.single[to_tax] = 1
                else:
                    self.single[to_tax] += 1

                tax = sorted([from_tax, to_tax])
                tag = "&".join(tax)
                if tag not in self.pair:
                    self.pair[tag] = 1
                else:
                    self.pair[tag] += 1
        # print ('bkp num is', len(self.bkps))
        f.close()

    def build_matrix(self):
        self.matrix = np.zeros((len(self.nodes), len(self.nodes)))
        nodes_index = {}
        for i in range(len(self.nodes)):
            nodes_index[self.nodes[i]] = i
        for pairs in self.pair.keys():
            pair_list = pairs.split("&")
            if pair_list[0] == pair_list[1]:
                continue
            a = nodes_index[pair_list[0]]
            b = nodes_index[pair_list[1]]
            self.matrix[a][b] = 1
            self.matrix[b][a] = 1
        self.matrix = np.matrix(self.matrix)

        

class Analyze():

    def __init__(self):
        self.data = []
        self.disease = ["CRC", "control", "adenoma"]
        self.disease_index = {"CRC":0, "control":1, "adenoma":2}
        self.disease_sample_num = {"CRC":0, "control":0, "adenoma":0}

        self.all_samples()

    def all_samples(self):
        all_acc_file = "acc.list"
        os.system(f"ls new_result/SRR*>{all_acc_file}")
        for line in open(all_acc_file):
            acc_file = line.strip()
            ID = acc_file.split("/")[1].split(".")[0]
            new_file = acc_file.split("/")[1]
            # os.system(f"python3 /mnt/d/breakpoints/script/remove_repeat.py {acc_file} new_result/{new_file}")
            acc_file = f"new_result/{new_file}"
            
            sample = Sample(acc_file, ID)
            # print (sample.cohort)
            if sample.cohort != "ThomasAM_2018b":
                continue
            self.data.append(sample)
            self.disease_sample_num[sample.disease] += 1
            # print (sample.cohort)
            # print (sample.ID, sample.disease, len(sample.filter_bkps))
            # break
        print ("sample num count", self.disease_sample_num)

    def taxonomy_count(self):
        
        level = 6
        count_dict = {}
        groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
        for sample in self.data:
            # if sample.disease != "control":   #different groups
            #     continue
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[level]
                to_tax = bkp.to_ref_lineage.split(";")[level]
                for tax in [from_tax, to_tax]:
                    if tax not in count_dict:
                        count_dict[tax] = 1
                    else:
                        count_dict[tax] += 1
                    if tax not in groups_cout_dict[sample.disease]:
                        groups_cout_dict[sample.disease][tax] = 1
                    else:
                        groups_cout_dict[sample.disease][tax] += 1
        for disease in groups_cout_dict.keys():
            print (disease, "Single Top")
            self.show_sort_dict(groups_cout_dict[disease])
        self.show_sort_dict(count_dict)

    def taxonomy_pair_count(self):
        level_dict = {"phylum":1, "genus":5, "species":6}
        level = 5
        count_dict = {}
        groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
        for sample in self.data:
            # if sample.disease != "control":   #different groups
            #     continue
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[level]
                to_tax = bkp.to_ref_lineage.split(";")[level]
                tax = sorted([from_tax, to_tax])
                tag = "&".join(tax)
                if tag not in count_dict:
                    count_dict[tag] = 1
                else:
                    count_dict[tag] += 1
                if tag not in groups_cout_dict[sample.disease]:
                    groups_cout_dict[sample.disease][tag] = 1
                else:
                    groups_cout_dict[sample.disease][tag] += 1
        for disease in groups_cout_dict.keys():
            print (disease, "Paired Top")
            self.show_sort_dict(groups_cout_dict[disease])
        self.show_sort_dict(count_dict)

    def taxonomy_pair_samples(self):
        level_dict = {"phylum":1, "genus":5, "species":6}
        level = 5
        count_dict = {}
        groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
        i = 0
        for sample in self.data:
            # if sample.disease != "control":   #different groups
            #     continue
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[level]
                to_tax = bkp.to_ref_lineage.split(";")[level]
                # print (bkp.from_ref_lineage, bkp.to_ref_lineage)
                tax = sorted([from_tax, to_tax])
                tag = "&".join(tax)
                if tag not in count_dict:
                    count_dict[tag] = [0] * len(self.data)
                count_dict[tag][i] = 1
                if tag not in groups_cout_dict[sample.disease]:
                    groups_cout_dict[sample.disease][tag] = [0] * len(self.data)
                groups_cout_dict[sample.disease][tag][i] = 1
            i += 1
        for key in count_dict:
            # print (count_dict[key])
            count_dict[key] = np.sum(count_dict[key])
        for key in groups_cout_dict:
            for k in groups_cout_dict[key].keys():
                groups_cout_dict[key][k] = np.sum(groups_cout_dict[key][k])
        for disease in groups_cout_dict.keys():
            print (disease, "Paired samples Top")
            self.show_sort_dict(groups_cout_dict[disease])
        self.show_sort_dict(count_dict)

    def show_sort_dict(self, count_dict):
        if len(count_dict) == 0:
            return 0
        sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
        for i in range(10):
            print (sorted_count[i])
        print ("**************")
        return sorted_count
    
    def check_CRC_related_species(self):
        CRC_related_species = ["Fusobacterium nucleatum", "Solobacterium moorei",\
         "Porphyromonas asaccharolytica","Parvimonas micra","Peptostreptococcus stomatis","Parvimonas spp"]
        for sample in self.data:
            for bkp in sample.filter_bkps:
                from_species = bkp.from_ref_lineage.split(";")[6]
                to_species = bkp.to_ref_lineage.split(";")[6]
                for CRC_species in CRC_related_species:
                    if re.search(CRC_species, from_species) or re.search(CRC_species, to_species):
                        print (sample.ID, sample.disease, from_species, to_species)

    def compare_HGT_complexity(self):
        crc = []
        control = []
        adenoma = []
        for sample in self.data:
            if sample.disease == "control":
                control.append(len(sample.filter_bkps))
            elif sample.disease == "CRC":
                crc.append(len(sample.filter_bkps))
            elif sample.disease == "adenoma":
                adenoma.append(len(sample.filter_bkps))
        p1 = self.compare_diff(crc, control)
        p2 = self.compare_diff(crc, adenoma)
        p3 = self.compare_diff(adenoma, control)
        print ("Mean: CRC:%s control:%s adenoma:%s"%(np.mean(crc), np.mean(control), np.mean(adenoma)), p1, p2, p3)
        print ("Median: CRC:%s control:%s adenoma:%s"%(np.median(crc), np.median(control), np.median(adenoma)))

        ################# level = 4 is best
        crc = []
        control = []
        adenoma = []
        crc_bases = []
        control_bases = []
        for sample in self.data:
            if sample.disease == "control":
                control.append(len(sample.single))
                control_bases.append(sample.bases)
            elif sample.disease == "CRC":
                crc.append(len(sample.single))
                crc_bases.append(sample.bases)
            elif sample.disease == "adenoma":
                adenoma.append(len(sample.single))
        p1 = self.compare_diff(crc, control)
        p2 = self.compare_diff(crc, adenoma)
        p3 = self.compare_diff(adenoma, control)
        p4 = self.compare_diff(crc_bases, control_bases)
        print ("Mean: CRC:%s control:%s adenoma:%s"%(np.mean(crc), np.mean(control), np.mean(adenoma)), p1, p2, p3, p4)
        print ("Median: CRC:%s control:%s adenoma:%s"%(np.median(crc), np.median(control), np.median(adenoma)))

        #################
        crc = []
        control = []
        adenoma = []
        for sample in self.data:
            if sample.disease == "control":
                control.append(len(sample.pair))
            elif sample.disease == "CRC":
                crc.append(len(sample.pair))
            elif sample.disease == "adenoma":
                adenoma.append(len(sample.pair))
        p1 = self.compare_diff(crc, control)
        p2 = self.compare_diff(crc, adenoma)
        p3 = self.compare_diff(adenoma, control)
        print ("Mean: CRC:%s control:%s adenoma:%s"%(np.mean(crc), np.mean(control), np.mean(adenoma)), p1, p2, p3)
        print ("Median: CRC:%s control:%s adenoma:%s"%(np.median(crc), np.median(control), np.median(adenoma)))

    def compare_diff(self, crc, control):
        if len(crc) == 0 or len(control) == 0:
            return 0
        U1, p = mannwhitneyu(crc, control, method="auto")
        print (stats.ttest_ind(crc, control))
        # print (ranksums(crc, control))
        return p

    def cohort_HGT(self):

        ##########level=4 best 
        groups_cout_dict = {}
        record_dict = {}
        for sample in self.data:
            for HGT in sample.pair:  
                if HGT not in groups_cout_dict:
                    groups_cout_dict[HGT] = [0, 0, 0]
                index = self.disease_index[sample.disease]
                groups_cout_dict[HGT][index] += 1
        for HGT in groups_cout_dict:
            record_dict[HGT] = abs(groups_cout_dict[HGT][0] - groups_cout_dict[HGT][1])
        sorted_count = self.show_sort_dict(record_dict)
        
        array = np.zeros((len(self.data), len(sorted_count)))
        for i in range(len(sorted_count)):
            # print (sorted_count[i], groups_cout_dict[sorted_count[i][0]])
            tag = sorted_count[i][0]
            j = 0
            for sample in self.data:
                if sample.disease == "CRC":
                    if tag in sample.pair:
                        array[j][i] = np.log(sample.pair[tag])
                    j += 1
            for sample in self.data:
                if sample.disease == "control":
                    if tag in sample.pair:
                        array[j][i] = np.log(sample.pair[tag])
                    j += 1
            for sample in self.data:
                if sample.disease == "adenoma":
                    if tag in sample.pair:
                        array[j][i] = np.log(sample.pair[tag])
                    j += 1
        ax = sns.heatmap(array)
        plt.savefig('test.pdf')

            
        print ("##################")

        groups_cout_dict = {}
        record_dict = {}
        for sample in self.data:
            for HGT in sample.single:  
                if HGT not in groups_cout_dict:
                    groups_cout_dict[HGT] = [0, 0, 0]
                index = self.disease_index[sample.disease]
                groups_cout_dict[HGT][index] += 1
        for HGT in groups_cout_dict:
            record_dict[HGT] = abs(groups_cout_dict[HGT][0] - groups_cout_dict[HGT][1])
        sorted_count = self.show_sort_dict(record_dict)
        for i in range(50):
            print (sorted_count[i], groups_cout_dict[sorted_count[i][0]])
        print ("##################")

        # array = np.zeros((len(self.data), len(sorted_count)))
        # for i in range(len(sorted_count)):
            
        #     HGT = sorted_count[i][0]
        #     print (sorted_count[i][0], groups_cout_dict[HGT])
        #     j = 0
        #     x = 0
        #     for sample in self.data:
        #         if sample.disease == "CRC":
        #             if HGT in sample.single:
        #                 array[j][i] = 1
        #                 x += 1
        #             j += 1
        #     for sample in self.data:
        #         if sample.disease == "control":
        #             if HGT in sample.single:
        #                 array[j][i] = 1
        #                 x += 1
        #             j += 1
        #     for sample in self.data:
        #         if sample.disease == "adenoma":
        #             if HGT in sample.single:
        #                 array[j][i] = 1
        #                 x += 1
        #             j += 1
        # # print (array)
        # ax = sns.heatmap(array)
        # plt.savefig('test.pdf')

    def sort_dict_by_value(self, count_dict):
        sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
        return sorted_count

    def per_species_HGT_comparison(self):
        control = []
        disease = []
        for sample in self.data: 
            if sample.disease == "control":
                control.append(sample.transtivity)
            if sample.disease == "CRC":
                disease.append(sample.transtivity)
        p1 = self.compare_diff(disease, control)
        print (len(disease), len(control), p1)

    def cohort_specific_HGT(self):

        ##########level=4 best 
        groups_cout_dict = {}
        record_dict = {}
        for sample in self.data:
            for HGT in sample.pair:  
                if HGT not in groups_cout_dict:
                    groups_cout_dict[HGT] = [0, 0, 0]
                index = self.disease_index[sample.disease]
                groups_cout_dict[HGT][index] += 1
        for HGT in groups_cout_dict:
            if groups_cout_dict[HGT][0] == 0:
                record_dict[HGT] = groups_cout_dict[HGT][1]
        sorted_count = self.show_sort_dict(record_dict)     

        groups_cout_dict = {}
        record_dict = {}
        for sample in self.data:
            for HGT in sample.pair:  
                if HGT not in groups_cout_dict:
                    groups_cout_dict[HGT] = [0, 0, 0]
                index = self.disease_index[sample.disease]
                groups_cout_dict[HGT][index] += 1
        for HGT in groups_cout_dict:
            if groups_cout_dict[HGT][1] == 0:
                record_dict[HGT] = groups_cout_dict[HGT][0]
        sorted_count = self.show_sort_dict(record_dict)   

    def taxonomy_barplot(self):
        data = []
        count_dict = {}
        for sample in self.data:
            # if sample.disease != "control":   #different groups
            #     continue
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[sample.level]
                to_tax = bkp.to_ref_lineage.split(";")[sample.level]
                for tax in [from_tax, to_tax]:
                    if tax not in count_dict:
                        count_dict[tax] = 1
                    else:
                        count_dict[tax] += 1
        sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
        class_list = [x[0] for x in sorted_count][:5] + ['other']     #list(count_dict.keys())
        index_dict = {}
        for i in range(len(class_list)):
            index_dict[class_list[i]] = i
        for sample in self.data:
            if sample.disease != "CRC":
                continue
            sample_list = [0]*len(class_list)
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[sample.level]
                to_tax = bkp.to_ref_lineage.split(";")[sample.level]
                for tax in [from_tax, to_tax]:
                    if tax not in index_dict:
                        sample_list[-1] += 1
                    else:
                        sample_list[index_dict[tax]] += 1
            sample_array = np.array(sample_list)/sum(sample_list)
            data.append(sample_array)
            print (sample_list)

        for sample in self.data:
            if sample.disease != "control":
                continue
            sample_list = [0]*len(class_list)
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[sample.level]
                to_tax = bkp.to_ref_lineage.split(";")[sample.level]
                for tax in [from_tax, to_tax]:
                    
                    if tax not in index_dict:
                        sample_list[-1] += 1
                    else:
                        sample_list[index_dict[tax]] += 1
            sample_array = np.array(sample_list)/sum(sample_list)
            data.append(sample_array)
            print (sample_list)
        df = pd.DataFrame(data, columns = class_list)
        print (df)
        df.plot(kind='bar', stacked=True,
        title='Stacked Bar Graph by dataframe')
        plt.savefig('bar.pdf')

    def taxonomy_circos(self):
        
        count_dict = {}
        for sample in self.data:
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref_lineage.split(";")[sample.level]
                to_tax = bkp.to_ref_lineage.split(";")[sample.level]
                if from_tax == to_tax:
                    continue
                for tax in [from_tax, to_tax]:
                    if tax not in count_dict:
                        count_dict[tax] = 1
                    else:
                        count_dict[tax] += 1
        sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
        class_list = [x[0] for x in sorted_count][:7] + ['other']     #list(count_dict.keys())
        index_dict = {}
        for i in range(len(class_list)):
            index_dict[class_list[i]] = i
        for co in sorted_count:
            if co[0] not in index_dict:
                index_dict[co[0]] = len(class_list)-1
        print (index_dict)

        data = []
        for i in range(len(class_list)):
            one_class = [0]*len(class_list)
            for sample in self.data:
                # if sample.disease != "CRC":
                #     continue
                
                for bkp in sample.filter_bkps:
                    from_tax = bkp.from_ref_lineage.split(";")[sample.level]
                    to_tax = bkp.to_ref_lineage.split(";")[sample.level]
                    if index_dict[from_tax] == i:
                        j = index_dict[to_tax]
                        if j < i:
                            one_class[j] += 1
                    if index_dict[to_tax] == i:
                        j = index_dict[from_tax]
                        if j < i:
                            one_class[j] += 1
            data.append(one_class)
        df = pd.DataFrame(data, columns = class_list)
        df.index = class_list
        df.to_csv('for_circos.csv', sep='\t')
        print (df)
        os.system("Rscript circos.R")

      

level_dict = {"phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
sra_meta = "italy.csv"
pheno_file = "CRC.xlsx"
UHGG_meta = "/mnt/d/breakpoints/HGT/UHGG/genomes-all_metadata.tsv"
phenotype = Phenotype()
taxonomy = Taxonomy()


analyze = Analyze()
# analyze.taxonomy_count()
# analyze.taxonomy_pair_count()
# analyze.check_CRC_related_species()
# analyze.compare_HGT_complexity()
# analyze.taxonomy_pair_samples()
# analyze.cohort_HGT()
# analyze.each_species()
# analyze.cohort_specific_HGT()
# analyze.per_species_HGT_comparison()
# analyze.taxonomy_barplot()
analyze.taxonomy_circos()









# print (taxonomy.taxonomy_dict["GUT_GENOME000002"])

# file = "run_result/SRR6915098.acc.csv"
# ID = "SRR6915098"
# sample = Sample(file, ID)
# sample.read_bkp()
# print (sample.disease)

# all_samples()
