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
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import roc_auc_score
from random import shuffle
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit

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
        self.from_bkp = int(list[1])
        self.to_bkp = int(list[3])
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
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/yu_2015.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/germany.csv")
        # self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/france/france.csv")
        self.read_sra_meta("/mnt/d/breakpoints/script/analysis/italy.csv")
        
        
    def read_sra_meta(self, sra_meta):
        # f = open(sra_meta)
        df = pd.read_csv(sra_meta)
        
        for i in range(len(df.index)):
            sra_ID = df["Run"][i]
            bases = df["Bases"][i]
            if "sample_name" in df.columns and df["sample_name"][i] != "Illumina":
                sample_name = df["sample_name"][i]
            else:
                sample_name = df["Sample Name"][i]
            # for col in df.columns:
            #     print (col, df[col][i])
                
            # print (sra_ID, bases, sample_name) 
            if sample_name not in self.name_disease:
                continue
            self.ID_disease[sra_ID] = self.name_disease[sample_name]
            self.ID_cohort[sra_ID] = self.name_cohort[sample_name]
            self.ID_bases[sra_ID] = bases
        # print (len(self.ID_bases))

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
        self.level = 5
        

        self.read_bkp()
        self.bkp_num = len(self.filter_bkps)
        self.average_bkp_per_species = round(self.bkp_num/len(self.single), 2)
        self.matrix = []
        self.select_feature_matrix = []
        self.select_feature_graph = []
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

    def given_nodes_make_matrix(self, common_nodes_dict, window):
        nodes = list(common_nodes_dict.keys())
        self.select_feature_matrix = np.zeros((len(common_nodes_dict), len(common_nodes_dict)))
        nodes_index = {}
        for i in range(len(nodes)):
            nodes_index[nodes[i]] = i
        for bkp in self.filter_bkps:
            from_tax = bkp.from_ref
            to_tax = bkp.to_ref
            tax = sorted([from_tax, to_tax])
            tag = "&".join(tax)
            if tax[0] == from_tax:
                from_b = bkp.from_bkp
                to_b = bkp.to_bkp
            else:
                to_b = bkp.from_bkp
                from_b = bkp.to_bkp 
            from_b = int(from_b/window)     
            to_b = int(to_b/window)
            array = tax + [from_b, to_b]
            node1 = array[0] + "&" + str(array[2])
            node2 = array[1] + "&" + str(array[3])
            if node1 in common_nodes_dict and node2 in common_nodes_dict:
                pass
            else:
                continue
            a = nodes_index[node1]
            b = nodes_index[node2]
            self.select_feature_matrix[a][b] = 1
            self.select_feature_matrix[b][a] = 1
        self.select_feature_matrix = np.matrix(self.select_feature_matrix)
        self.select_feature_graph = nx.from_numpy_matrix(self.select_feature_matrix)

    def build_matrix(self):
        self.matrix = np.zeros((len(self.nodes), len(self.nodes)))
        nodes_index = {}
        for i in range(len(self.nodes)):
            nodes_index[self.nodes[i]] = i
        for pairs in self.pair.keys():
            pair_list = pairs.split("&")
            if pair_list[0] == pair_list[1]:
                continue
            if self.check_valid_name(pair_list[0]) == False or self.check_valid_name(pair_list[1]) == False:
                continue  #remove s__
            # print (pair_list[0], self.check_valid_name(pair_list[0]))
            a = nodes_index[pair_list[0]]
            b = nodes_index[pair_list[1]]
            self.matrix[a][b] = 1
            self.matrix[b][a] = 1
        self.matrix = np.matrix(self.matrix)

    def check_valid_name(self, name):
        if name.split("__")[1] != '':
            return True
        else:
            return False
        
class Analyze():

    def __init__(self):
        self.data = []
        self.disease = ["CRC", "control", "adenoma"]
        self.disease_index = {"CRC":0, "control":1, "adenoma":2}
        self.disease_sample_num = {"CRC":0, "control":0, "adenoma":0}
        self.disease_sample_num_cohort = {}

        self.all_samples()
        self.window = 100

    def all_samples(self):
        all_acc_file = "acc.list"
        os.system(f"ls new_result/*acc.csv>{all_acc_file}")
        for line in open(all_acc_file):
            acc_file = line.strip()
            ID = acc_file.split("/")[1].split(".")[0]
            new_file = acc_file.split("/")[1]
            # os.system(f"python3 /mnt/d/breakpoints/script/remove_repeat.py {acc_file} new_result/{new_file}")
            acc_file = f"new_result/{new_file}"
            
            sample = Sample(acc_file, ID)
            # print (sample.cohort)
            if sample.disease == "adenoma":
                continue
            self.data.append(sample)
            self.disease_sample_num[sample.disease] += 1
            if sample.cohort not in self.disease_sample_num_cohort:
                self.disease_sample_num_cohort[sample.cohort] = {"CRC":0, "control":0, "adenoma":0}
            self.disease_sample_num_cohort[sample.cohort][sample.disease] += 1
            # print (sample.cohort)
            # print (sample.ID, sample.disease, len(sample.filter_bkps))
            # break
        for key in self.disease_sample_num_cohort:
            print ("Cohort", key, self.disease_sample_num_cohort[key])
        print ("Total", self.disease_sample_num)

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
        data = []
        control_transtivity = []
        disease_transtivity = []

        control_graph_density = []
        disease_graph_density = []

        control_nodes = []
        disease_nodes = []

        control_pair = []
        disease_pair = []

        for sample in self.data: 
            data.append([sample.disease, sample.transtivity, sample.graph_density, len(sample.nodes), len(sample.pair)])

            if sample.disease == "control":
                control_transtivity.append(sample.transtivity)
                control_graph_density.append(sample.graph_density)
                control_nodes.append(len(sample.nodes))
                control_pair.append(len(sample.pair))
            if sample.disease == "CRC":
                disease_transtivity.append(sample.transtivity)
                disease_graph_density.append(sample.graph_density)
                disease_nodes.append( len(sample.nodes))
                disease_pair.append(len(sample.pair))
        p1 = self.compare_diff(disease_transtivity, control_transtivity)
        p2 = self.compare_diff(disease_graph_density, control_graph_density)
        p3 = self.compare_diff(disease_nodes, control_nodes)
        p4 = self.compare_diff(disease_pair, control_pair)
        print (len(control_transtivity), len(disease_transtivity), p1, p2, p3, p4)
        # print (np.mean(control_graph_density), np.std(control_graph_density))
        # print (np.mean(disease_graph_density), np.std(disease_graph_density))
        df = pd.DataFrame(data, columns = ["Condition", "Transtivity", "Density", "Nodes_Number", "Pair_Number"])
        df.to_csv('for_box.csv', sep='\t')

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

    def plot_most_common_pair(self):
        count_dict = {}
        record_pos = {}

        condition = ["CRC"]
        for sample in self.data:
            sample_dict = {}
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref
                to_tax = bkp.to_ref
                tax = sorted([from_tax, to_tax])
                tag = "&".join(tax)
                if tax[0] == bkp.from_ref:
                    from_b = bkp.from_bkp
                    to_b = bkp.to_bkp
                else:
                    to_b = bkp.from_bkp
                    from_b = bkp.to_bkp 
                from_b = int(from_b/self.window)
                to_b = int(to_b/self.window)
                if tag not in count_dict:
                    if tag not in sample_dict:
                        count_dict[tag] = 1
                        sample_dict[tag] = 1
                    record_pos[tag] = [[from_b, to_b, sample.disease, 1]]
                else:
                    if tag not in sample_dict:
                        count_dict[tag] += 1
                        sample_dict[tag] = 1


                    flag = False
                    for point in record_pos[tag]:
                        if point[0]==  from_b and point[1] == to_b and sample.disease == point[2]:
                            point[3] += 1
                            flag = True
                            break
                    if flag == False:
                        record_pos[tag].append([from_b, to_b, sample.disease, 1])


        sorted_count_dict = self.show_sort_dict(count_dict)  
        # print (record_pos[sorted_count_dict[0][0]])
        # top_index = 14
        for i in range(len(sorted_count_dict)):
            if sorted_count_dict[i][0] == "GUT_GENOME144544_1&GUT_GENOME145378_7":
                top_index = i
                print (top_index)

        # for top_index in range(100):
        df = pd.DataFrame(record_pos[sorted_count_dict[top_index][0]], columns =sorted_count_dict[top_index][0].split("&")+["Condition", "Sample_Num"] )
        # df = pd.DataFrame(record_pos[sorted_count_dict[top_index][0]], columns =["G1", "G2", "Condition", "Sample_Num"] )
        # print (df)
        df.to_csv('for_para.csv', sep='\t')
        os.system("Rscript parallel.R")
        os.system("mv para.pdf para_%s.pdf"%(top_index))
        #     # break
        return sorted_count_dict, record_pos

    def bkp_pair_count(self):
        level_dict = {"phylum":1, "genus":5, "species":6}
        genome_lineag = {}
        level = 5
        count_dict = {}
        window = 500
        groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
        for sample in self.data:
            # if sample.disease != "control":   #different groups
            #     continue
            for bkp in sample.filter_bkps:
                from_tax = bkp.from_ref
                to_tax = bkp.to_ref
                genome_lineag[from_tax] = bkp.from_ref_lineage.split(";")[level]
                genome_lineag[to_tax] = bkp.to_ref_lineage.split(";")[level]
                tax = sorted([from_tax, to_tax])

                tag = "&".join(tax)
                if tax[0] == bkp.from_ref:
                    from_b = bkp.from_bkp
                    to_b = bkp.to_bkp
                else:
                    to_b = bkp.from_bkp
                    from_b = bkp.to_bkp 

                from_b = int(from_b/window)     
                to_b = int(to_b/window)

                tag = "&".join(tax) + "&" +str(from_b) + "&"+ str(to_b)
                if tag not in count_dict:
                    count_dict[tag] = 1
                else:
                    count_dict[tag] += 1
                if tag not in groups_cout_dict[sample.disease]:
                    groups_cout_dict[sample.disease][tag] = 1
                else:
                    groups_cout_dict[sample.disease][tag] += 1
        # for disease in groups_cout_dict.keys():
        #     print (disease, "Paired Top")
        #     self.show_sort_dict(groups_cout_dict[disease])
        # self.show_sort_dict(count_dict)   
        

        CRC_specific_HGT = {} 
        for locus in groups_cout_dict["CRC"].keys():
            if locus not in groups_cout_dict["control"]:
                CRC_specific_HGT[locus] = groups_cout_dict["CRC"][locus]
                # print (locus, groups_cout_dict["CRC"][locus])
        sort_CRC_specific_HGT = self.show_sort_dict(CRC_specific_HGT)  

        # sorted_count_dict, record_pos = self.plot_most_common_pair()
        for pair in sort_CRC_specific_HGT[:20]:
            fir = pair[0].split("&")[0]
            sec = pair[0].split("&")[1]
            print (genome_lineag[fir], genome_lineag[sec])
            """
            tag = "&".join(pair[0].split("&")[:2])
            for i in range(len(sorted_count_dict)):
                if sorted_count_dict[i][0] == tag:
                    top_index = i            
            df = pd.DataFrame(record_pos[sorted_count_dict[top_index][0]], columns =sorted_count_dict[top_index][0].split("&")+["Condition", "Sample_Num"] )
            # df = pd.DataFrame(record_pos[sorted_count_dict[top_index][0]], columns =["G1", "G2", "Condition", "Sample_Num"] )
            # print (df)
            df.to_csv('for_para.csv', sep='\t')
            os.system("Rscript parallel.R")
            os.system("mv para.pdf para_%s.pdf"%(top_index))
            """

        control_specific_HGT = {} 
        for locus in groups_cout_dict["control"].keys():
            if locus not in groups_cout_dict["CRC"]:
                control_specific_HGT[locus] = groups_cout_dict["control"][locus]
                # print (locus, groups_cout_dict["CRC"][locus])
        sort_control_specific_HGT = self.show_sort_dict(control_specific_HGT) 
        for pair in sort_control_specific_HGT[:20]:
            fir = pair[0].split("&")[0]
            sec = pair[0].split("&")[1]
            print (genome_lineag[fir], genome_lineag[sec]) 

class RF():
    def __init__(self):
        analyze = Analyze()
        self.all_data = analyze.data
        # shuffle(self.all_data)
        self.diff_cohorts = self.classify_cohorts()        
        self.feature_num = 50
        self.window = 500
        self.level = 5
        self.remove_adenoma = True

    def classify_cohorts(self):
        diff_cohorts = {}
        for sample in self.all_data:
            # if sample.cohort == "ThomasAM_2018a":
            #     continue
            if sample.cohort not in diff_cohorts:
                diff_cohorts[sample.cohort] = [sample]
            else:
                diff_cohorts[sample.cohort].append(sample)
        return diff_cohorts

    def select_tag(self, bkp):
        # from_tax = bkp.from_ref_lineage.split(";")[self.level] #  bkp.from_ref
        # to_tax = bkp.to_ref_lineage.split(";")[self.level]  #bkp.to_ref
        from_tax = bkp.from_ref
        to_tax = bkp.to_ref
        tax = sorted([from_tax, to_tax])
        tag = "&".join(tax)
        if tax[0] == from_tax:
            from_b = bkp.from_bkp
            to_b = bkp.to_bkp
        else:
            to_b = bkp.from_bkp
            from_b = bkp.to_bkp 
        from_b = int(from_b/self.window)     
        to_b = int(to_b/self.window)
        return "&".join(tax) + "&" +str(from_b) + "&"+ str(to_b)

    def select_feature(self, cohort_data):    
        specific_HGT = {} 
        train_num = len(cohort_data)

        crc_num = 0
        control_num = 0
        for sample in cohort_data:
            if sample.disease == "CRC":
                crc_num += 1
            if sample.disease == "control":
                control_num += 1
        print (train_num, "CRC:",crc_num, "control:",control_num)
        i, j = 0, 0
        for sample in cohort_data:
            for bkp in sample.filter_bkps:
                tag = self.select_tag(bkp)
                if tag not in specific_HGT:
                    specific_HGT[tag] = [[0]*crc_num, [0]*control_num]
                if sample.disease == "CRC":
                    specific_HGT[tag][0][i] += 1
                if sample.disease == "control":
                    specific_HGT[tag][1][j] += 1

            if sample.disease == "CRC":
                i += 1
            if sample.disease == "control":
                j += 1
        print (len(specific_HGT))
        p_specific_HGT = {}
        for tag in specific_HGT:
            crc = specific_HGT[tag][0]
            control = specific_HGT[tag][1]
            if float(sum(crc) + sum(control))/train_num > 0.1:
                U1, p = mannwhitneyu(crc, control, method="auto")
                p_specific_HGT[tag] = p
            # else:
            #     print (crc, control, float(sum(crc) + sum(control))/train_num)
        print ("t-test done", len(p_specific_HGT))
        self.feature_num = len(p_specific_HGT) # use all features

        sort_specific_HGT = sorted(p_specific_HGT.items(), key=lambda item: item[1], reverse = False)[:self.feature_num]
        print ("sort done")
        sort_specific_HGT_dict = {}
        for i in range(len(sort_specific_HGT)):
            # print (sort_specific_HGT[i])
            sort_specific_HGT_dict[sort_specific_HGT[i][0]] = i
        with open('saved_dictionary.pkl', 'wb') as f:
            pickle.dump(sort_specific_HGT_dict, f)
        return sort_specific_HGT_dict

    def load_dict(self):
        with open('saved_dictionary.pkl', 'rb') as f:
            loaded_dict = pickle.load(f)
            return loaded_dict

    def generate_data(self, cohort_data, flag):
        if flag == "train":
            sort_specific_HGT_dict = self.select_feature(cohort_data)
        elif flag == "test":
            sort_specific_HGT_dict = self.load_dict() #use same features
        else:
            print ("no data")
        data = []
        label = []
        for sample in cohort_data:
            if sample.disease == "adenoma" and self.remove_adenoma:   
                continue
            sample_array = [0]*self.feature_num
            
            for bkp in sample.filter_bkps:
                tag = self.select_tag(bkp)
                if tag in sort_specific_HGT_dict:
                    sample_array[sort_specific_HGT_dict[tag]] += 1
            data.append(sample_array)
            # print (sample_array, sample.disease)
            label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

    def feature_matrix(self):
        data, label = self.generate_data(self.all_data)
        ax = sns.heatmap(data)
        plt.savefig('heatmap.pdf')

    def LODO(self):
        for lack in range(len(self.diff_cohorts)):
            train_sam, test_sam = self.split_test_train(lack)

            train_data, train_label =  self.generate_data(train_sam, "train")
            test_data, test_label = self.generate_data(test_sam, "test") 
            clf = RandomForestClassifier() #max_depth=2, random_state=0
            clf.fit(train_data, train_label)     
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            print (datasets[lack] , "AUC", roc_auc)  
            print ("*************************")
                
    def random_forest(self):

        train_data, train_label = self.generate_data(self.diff_cohorts["WirbelJ_2018"] + self.diff_cohorts["ThomasAM_2018a"]+self.diff_cohorts["ThomasAM_2018b"], "train")
        clf = RandomForestClassifier() #max_depth=2, random_state=0
        clf.fit(train_data, train_label)
        print ("training is done")


        test_data, test_label = self.generate_data(self.diff_cohorts["YuJ_2015"], "test")       
        roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
        print ("AUC", roc_auc)

        # cv = ShuffleSplit(n_splits=10, test_size=0.2, random_state=0)
        # scores = cross_val_score(clf, fir_data, fir_data, cv=4, scoring='roc_auc')
        # print (scores, np.mean(scores))

        # clf.fit(sec_data, sec_label)
        # roc_auc = roc_auc_score(fir_label, clf.predict_proba(fir_data)[:,1])
        # print ("AUC", roc_auc)

        # correct = 0
        # for i in range(len(test_data)):
        #     predict_label = clf.predict([test_data[i]])[0]
        #     if predict_label == test_label[i]:
        #         correct += 1
        #     print ("Predict label is %s; True label is %s."%(predict_label, test_label[i]))
        # roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
        # print ("AUC", roc_auc, correct/len(test_data), len(test_data))

    def split_test_train(self, lack):
        datasets = list(self.diff_cohorts.keys())
        train_sam = []
        test_sam = []
        for i in range(len(datasets)):
            if lack == i:
                test_sam = self.diff_cohorts[datasets[i]]
            else:
                train_sam += self.diff_cohorts[datasets[i]]  
        return train_sam, test_sam       

    def select_common_nodes(self, cohort_data):
        specific_HGT = {} 
        train_num = len(cohort_data)

        crc_num = 0
        control_num = 0
        for sample in cohort_data:
            if sample.disease == "CRC":
                crc_num += 1
            if sample.disease == "control":
                control_num += 1
        print (train_num, "CRC:",crc_num, "control:",control_num)
        i, j = 0, 0
        for sample in cohort_data:
            for bkp in sample.filter_bkps:
                tag = self.select_tag(bkp)
                if tag not in specific_HGT:
                    specific_HGT[tag] = [[0]*crc_num, [0]*control_num]
                if sample.disease == "CRC":
                    specific_HGT[tag][0][i] += 1
                if sample.disease == "control":
                    specific_HGT[tag][1][j] += 1

            if sample.disease == "CRC":
                i += 1
            if sample.disease == "control":
                j += 1

        p_specific_HGT = {}
        for tag in specific_HGT:
            crc = specific_HGT[tag][0]
            control = specific_HGT[tag][1]
            # if float(sum(crc) + sum(control))/train_num > 0.05:
            if max([float(sum(crc))/len(crc), float(sum(control))/len(control)]) > 0.2:

                U1, p = mannwhitneyu(crc, control, method="auto")
                p_specific_HGT[tag] =  p
        print ("edges num:", len(p_specific_HGT))
        sort_specific_HGT = sorted(p_specific_HGT.items(), key=lambda item: item[1], reverse = False)[:self.feature_num]
        # sort_specific_HGT = p_specific_HGT.items()

        sort_specific_HGT_dict = {}
        select_nodes = {}
        for i in range(len(sort_specific_HGT)):
            # print (sort_specific_HGT[i])
            sort_specific_HGT_dict[sort_specific_HGT[i][0]] = i
            edge = sort_specific_HGT[i][0]
            array = edge.split("&")
            node1 = array[0] + "&" + array[2]
            node2 = array[1] + "&" + array[3]
            select_nodes[node1] = 1
            select_nodes[node2] = 1
        return select_nodes

    def complex_feature(self):
        datasets = list(self.diff_cohorts.keys())
        for lack in range(len(self.diff_cohorts)):
            train_sam, test_sam = self.split_test_train(lack)
            features_dict = self.select_common_nodes(train_sam)
            for sample in self.all_data:
                sample.given_nodes_make_matrix(features_dict, self.window)

            train_data, train_label =  self.complex_data(train_sam)
            test_data, test_label = self.complex_data(test_sam) 
            clf = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5, ) #max_depth=2, random_state=0
            clf.fit(train_data, train_label)     
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            print (datasets[lack] , len(self.diff_cohorts[datasets[lack]]), "AUC", roc_auc)  
            print ("*************************")

    def complex_data(self, cohort_data):
        data = []
        label = []
        for sample in cohort_data:
            
            density = nx.density(sample.select_feature_graph)
            transitivity = nx.transitivity(sample.select_feature_graph)
            sample_array = np.array(sample.select_feature_matrix).flatten()   #[density, transitivity]
            # sample_array = [density, transitivity]
            # print ("xx", sample_array)
            data.append(sample_array)
            # print (sample_array, sample.disease)
            label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label



      
if __name__ == "__main__":
    level_dict = {"phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
    # sra_meta = "italy.csv"
    pheno_file = "allmetadata.xlsx"#"CRC.xlsx"
    UHGG_meta = "/mnt/d/breakpoints/HGT/UHGG/genomes-all_metadata.tsv"
    phenotype = Phenotype()
    taxonomy = Taxonomy()

    # cohort = "ThomasAM_2018b"
    # analyze = Analyze()
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
    # analyze.taxonomy_circos()
    # analyze.plot_most_common_pair()
    # analyze.bkp_pair_count()

    rf = RF()
    # rf.random_forest()
    # rf.feature_matrix()
    # rf.LODO()
    rf.complex_feature()
