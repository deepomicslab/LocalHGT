#!/usr/bin/env python3

import csv
import pysam
import random
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
from scipy import linalg
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums
import networkx as nx
import math
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA
from random import shuffle
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.feature_selection import SelectFromModel
from scipy.sparse import csgraph
from graph_denoise import svd_denoise
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
np.set_printoptions(threshold=sys.maxsize)

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
        self.score = float(list[9])
        # print (self.from_ref_genome, self.from_ref_lineage, self.to_ref_genome, self.to_ref_lineage)

    def print_out(self):
        print (self.from_ref, self.from_bkp, self.to_ref, self.to_bkp, self.from_side,\
         self.to_side, self.if_reverse, self.similarity)

    def write_out(self, writer):
        writer.writerow ([self.from_ref, self.from_bkp, self.to_ref, self.to_bkp, \
        self.from_side, self.to_side, self.if_reverse, self.read_str, self.ref_str, self.similarity])

def get_abd():
    marker_species_dict = {}
    i = 0
    for species in marker_species:
        if species.split()[-1] == "spp.":
            form_name = "g__" + species.split()[0] 
        else:
            form_name = "s__" + "_".join(species.split()) 
        # print (form_name)
        marker_species_dict[form_name] = i 
        i += 1
    print (marker_species_dict)
    return marker_species_dict

def get_samples_abd():
    sample_abd = {}
    marker_species_dict = get_abd()
    marker_species_num = len(marker_species)
    for cohort in cohort_abd:
        abd_file = cohort_abd[cohort]
        f = open("use/" + abd_file, 'r')
        i = 0
        for line in f:
            array = line.strip().split()
            if i == 0:
                sample_list = array
                for sample in sample_list:
                    sample_abd[sample] = [0] * marker_species_num
            else:
                species_name = array[0].split("|")[-1]
                
                if species_name in marker_species_dict:
                    species_index = marker_species_dict[species_name]
                    for j in range(len(sample_list)):
                        sample = sample_list[j]
                        abundance = float(array[j+1])
                        sample_abd[sample][species_index] = abundance
                else:
                    genus_name = species_name = array[0].split("|")[-2]
                    if genus_name in marker_species_dict:
                        species_index = marker_species_dict[genus_name]
                        for j in range(len(sample_list)):
                            sample = sample_list[j]
                            abundance = float(array[j+1])
                            sample_abd[sample][species_index] += abundance
                            # print (sample_abd[sample])
                        
                # print (len(sample_list), len(array))
            i += 1
    # print (sample_abd)
    return sample_abd

def get_abd_file_name():
    ID_abd_file = {} 
    file = "last_gutmeta_sample.tsv"
    df = pd.read_csv(file, sep = "\t")
    for i in range(len(df.index)):
        sra_ID = df["run_name"][i]
        sample_name = df["sample_name"][i]
        project_name  = df["project_name"][i]
        if pd.isnull(sra_ID):
            continue
        #     print ("******")
        # print (sra_ID, sample_name, project_name)
        array = sra_ID.split(";")
        for ID in array:
            ID = ID.strip()
            ID_abd_file[ID] = project_name + "_" + sample_name
            # print (ID, project_name + "_" + sample_name)
        # print (sra_ID, sample_name, project_name)
    return ID_abd_file
      
class Phenotype():
    def __init__(self):
        self.name_disease = {}
        self.name_cohort = {}
        self.name_bases = {}
        self.name_basics = {}
        self.ID_disease = {}
        self.ID_cohort = {}
        self.ID_bases = {}
        self.ID_basics = {}
        self.name_sample_id = {}
        self.ID_marker_abundance = {}
        self.name_marker_abundance = {}
        self.read_pheno()
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/USA/usa.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/japan.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/yu_2015.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/germany.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/france/france.csv")
        self.read_sra_meta("/mnt/d/breakpoints/script/analysis/italy.csv")
        self.read_sra_meta("/mnt/d/breakpoints/script/analysis/new_result/usa_canada.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/austria/austria.csv")     
        
    def read_sra_meta(self, sra_meta, for_name = "no"):
        # f = open(sra_meta)
        df = pd.read_csv(sra_meta)
        for i in range(len(df.index)):
            sra_ID = df["Run"][i]
            # bases = float(df["Bases"][i])
            if "sample_name" in df.columns and df["sample_name"][i] != "Illumina" and df["BioProject"][i] != "PRJDB4176" :
                sample_name = df["sample_name"][i]
            else:
                sample_name = df["Sample Name"][i]
            if sra_meta == "/mnt/d/breakpoints/HGT/CRC/austria/austria.csv":
                sample_name = "SID" + str(df["Submitter_Id"][i])

            study = df["SRA Study"][i]
            # for col in df.columns:
            #     print (col, df[col][i])
            
            if sample_name not in self.name_disease:
                continue

            self.ID_disease[sra_ID] = self.name_disease[sample_name]
            self.ID_cohort[sra_ID] = self.name_cohort[sample_name]
            self.ID_bases[sra_ID] = self.name_bases[sample_name]
            self.ID_basics[sra_ID] = self.name_basics[sample_name]
            if sample_name not in sample_abd:
                sample_id = self.name_sample_id[sample_name]
                abundance = sample_abd[sample_id]

            else:
                abundance = sample_abd[sample_name]
            self.ID_marker_abundance[sra_ID] = abundance
            self.name_marker_abundance[sample_name] = abundance


    def read_pheno(self):     
        df = pd.read_excel(pheno_file, header=None) 
        # print (df)
        # crc = 0
        # control = 0 
        for index, row in df.iterrows():
            if index == 0:
                continue
            sample_name = row[3]
            sample_id = row[2]
            condition = row[6]
            full_disease = row[7]
            cohort = row[1]
            bases = float(row[18])
            age = str(row[8])
            BMI = str(row[25])
            gender = gender_dict[str(row[11])]
            # print (index, age, gender, BMI)
            if age == "nan":
                age = 0
            else:
                age = int(age)
            

            # print (BMI)
            if BMI == "nan":
                BMI = 0
            else:
                # print (BMI)
                BMI = round(float(BMI))
            if cohort == "ZellerG_2014" or cohort == "YachidaS_2019" or cohort == "HanniganGD_2017":
                sample_name = row[2]

            # if cohort == "VogtmannE_2016":
            #     print (cohort, sample_name, condition, row[7])

            self.name_disease[sample_name] = condition
            self.name_cohort[sample_name] = cohort
            self.name_basics[sample_name] = [age, gender, BMI]
            self.name_bases[sample_name] = bases
            self.name_sample_id[sample_name] = sample_id
            # print (sample_name,condition) 
        # print (crc, control)

class Sample():

    def __init__(self, bkp_file, ID, level):
        self.bkp_file = bkp_file
        self.bkps = []
        self.filter_bkps = [] #recheck, not the same species
        self.ID = ID
        self.tag = "normal"
        
        if ID in phenotype.ID_disease:
            self.disease = phenotype.ID_disease[ID]
            self.cohort = phenotype.ID_cohort[ID]    
            self.bases = phenotype.ID_bases[ID]   
            self.basic_features =  phenotype.ID_basics[ID] 
            self.marker_abundance = phenotype.ID_marker_abundance[ID]
        elif ID in phenotype.name_disease:
            self.disease = phenotype.name_disease[ID]
            self.cohort = phenotype.name_cohort[ID]    
            self.bases = phenotype.name_bases[ID]    
            self.basic_features = phenotype.name_basics[ID] 
            self.marker_abundance = phenotype.name_marker_abundance[ID] 
        else:
            print ("no pheno", ID) 
            self.tag = "no pheno"    
        self.single = {}
        self.pair = {}
        # self.level = 5
        self.level = level
        

        self.read_bkp()
        self.bkp_num = len(self.filter_bkps)
        if self.bkp_num == 0:
            self.tag = "no bkp"
            print ("no bkp", ID)
        if self.tag != "no bkp":
            pass
            # """
            # self.average_bkp_per_species = round(self.bkp_num/len(self.single), 2)
            # self.matrix = []
            # self.select_feature_matrix = []
            # self.select_feature_graph = []
            self.select_feature_array = []
            # self.nodes = list(self.single.keys()) #node num significant 4
            # self.build_matrix()
            # self.graph = nx.from_numpy_matrix(self.matrix)
            # self.graph_density = nx.density(self.graph) #significant 4 5 6
            # self.graph_average_clustering = nx.average_clustering(self.graph, nodes=None, weight=None, count_zeros=True) #not
            # self.transtivity = nx.transitivity(self.graph) #significant 6
            # self.weighted_degree = self.matrix.sum(axis=1, dtype='float') #list
            # print (self.weighted_degree)
            # """
        # print ("one sample")

    def read_bkp(self):
        window = 500
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
    
                if self.level == 7:
                    from_tax = eb.from_ref
                    to_tax = eb.to_ref
                    tax = sorted([from_tax, to_tax])
                    tag = "&".join(tax)
                    if tax[0] == from_tax:
                        from_b = eb.from_bkp
                        to_b = eb.to_bkp
                    else:
                        to_b = eb.from_bkp
                        from_b = eb.to_bkp 
                    from_b = int(from_b/window)     
                    to_b = int(to_b/window)
                    new_tag = "&".join(tax) + "&" +str(from_b) + "&"+ str(to_b)
                    array = tax + [from_b, to_b]
                    from_tax = array[0] + "&" + str(array[2])
                    to_tax = array[1] + "&" + str(array[3])
                else:
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

    def given_nodes_make_matrix(self, common_nodes_dict, window, select_edges):
        nodes = list(common_nodes_dict.keys())
        n = len(common_nodes_dict)
        # self.select_feature_matrix = np.zeros((n, n))
        self.select_feature_array = np.zeros(len(select_edges))
        nodes_index = {}
        for i in range(len(nodes)):
            nodes_index[nodes[i]] = i
        for bkp in self.filter_bkps:
            
            if self.level == 7:
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
                new_tag = "&".join(tax) + "&" +str(from_b) + "&"+ str(to_b)
                array = tax + [from_b, to_b]
                node1 = array[0] + "&" + str(array[2])
                node2 = array[1] + "&" + str(array[3])
            else:
                from_tax = bkp.from_ref_lineage.split(";")[self.level] #  bkp.from_ref
                to_tax = bkp.to_ref_lineage.split(";")[self.level]  #bkp.to_ref
                tax = sorted([from_tax, to_tax])
                new_tag = "&".join(tax)       
                node1 = from_tax  
                node2 = to_tax   
            if node1 in common_nodes_dict and node2 in common_nodes_dict:
                pass
            else:
                continue
            
            if new_tag in select_edges:
                self.select_feature_array[select_edges[new_tag]] = 1

    def get_species_node_feature(self, nodes_dict, window, select_edges):
        common_nodes_dict = nodes_dict.copy()
        node_index = {}
        i = 0
        for node in common_nodes_dict:
            node_index[node] = i
            i += 1

        
        for node in common_nodes_dict:
            common_nodes_dict[node] = []

        for bkp in self.filter_bkps:    
            from_tax = bkp.from_ref_lineage.split(";")[self.level] #  bkp.from_ref
            to_tax = bkp.to_ref_lineage.split(";")[self.level]  #bkp.to_ref      
            node1 = from_tax  
            node2 = to_tax   
            if node1 in common_nodes_dict and node2 in common_nodes_dict:
                pass
            else:
                continue
            if node2 not in common_nodes_dict[node1]:
                common_nodes_dict[node1].append(node2)
            if node1 not in common_nodes_dict[node2]:
                common_nodes_dict[node2].append(node1)

        new_common_nodes_dict = {}
        for node in common_nodes_dict:
            new_common_nodes_dict[node] = np.zeros(len(common_nodes_dict))
            for a in common_nodes_dict[node]:
                new_common_nodes_dict[node][node_index[a]] = 1
                for b in common_nodes_dict[a]:
                    new_common_nodes_dict[node][node_index[b]] = 1



        for node in common_nodes_dict:
            common_nodes_dict[node] = np.sum(new_common_nodes_dict[node])
        self.select_feature_array = list(common_nodes_dict.values())
        # print (self.select_feature_array)
        # print (len(common_nodes_dict), self.select_feature_array)
        
    def check_graph(self):
        select_edges = {}
        for bkp in self.filter_bkps:
            from_tax = bkp.from_ref_lineage.split(";")[self.level] #  bkp.from_ref
            to_tax = bkp.to_ref_lineage.split(";")[self.level]  #bkp.to_ref
            tax = sorted([from_tax, to_tax])
            new_tag = "&".join(tax)    
            if new_tag not in select_edges:
                select_edges[new_tag] = 1

        event_array = np.zeros((len(select_edges), len(select_edges)))
        events = list(select_edges.keys())
        for i in range(len(events)):
            for j in range(len(events)):
                event_1 = events[i]
                event_2 = events[j]
                tag_array_1 = event_1.split("&")
                tag_array_2 = event_2.split("&")
                for tag in tag_array_1:
                    if tag in tag_array_2:
                        # print (event_1, event_2)
                        event_array[i][j] = 1
        e_num = 0
        v_num = len(event_array)
        
        for array in event_array:
            e_num += np.sum(array)
        e_num = (e_num-v_num) / 2
        mean_degree = e_num/v_num
        print (self.ID, len(select_edges), v_num, e_num, mean_degree)

    def build_matrix(self):
        self.matrix = np.zeros((len(self.nodes), len(self.nodes)))
        nodes_index = {}
        for i in range(len(self.nodes)):
            nodes_index[self.nodes[i]] = i
        for pairs in self.pair.keys():
            pair_list = pairs.split("&")
            if self.level == 7:
                node1 = pair_list[0] + "&" + pair_list[1]
                node2 = pair_list[2] + "&" + pair_list[3]
            else:
                node1 = pair_list[0]
                node2 = pair_list[1]
            if node1 == node2:
                continue
            if self.check_valid_name(node1) == False or self.check_valid_name(node2) == False:
                continue  #remove s__
            # print (pair_list[0], self.check_valid_name(pair_list[0]))
            a = nodes_index[node1]
            b = nodes_index[node2]
            self.matrix[a][b] = 1
            self.matrix[b][a] = 1
        self.matrix = np.matrix(self.matrix)

    def check_valid_name(self, name):
        if self.level < 7:
            if name.split("__")[1] != '':
                return True
            else:
                return False
        else:
            return True
        
class Analyze():

    def __init__(self, level):
        self.data = []
        self.level = level
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
            
            sample = Sample(acc_file, ID, self.level)
            if sample.tag == "no pheno" or sample.tag == "no bkp":
                continue            
            if (sample.disease != pheno_group_1) and (sample.disease != pheno_group_2):
                continue
            if str(sample.disease) == "nan":
                continue
            if sample.cohort == "HanniganGD_2017":
                continue
            # if sample.bases < 2000000000:
            #     continue
            # sample.check_graph()
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

    def __del__(self):
        del self.data
        del self.disease_sample_num_cohort
        print ("object, deleted")

class RF1():

    def __init__(self, level, common_ratio):
        analyze = Analyze(level)
        self.all_data = analyze.data
        del analyze
      
        self.window = 500
        self.common_ratio = common_ratio #0.01
        self.level = self.all_data[0].level
        self.remove_adenoma = True
        self.event_array = ''
        print ("RF init done")

    def prepare_data(self):
        event_dict = {}
        i = 0
        for sample in self.all_data:
            for bkp in sample.filter_bkps:
                edge = self.select_tag(bkp)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge not in event_dict:
                    event_dict[edge] = i
                    i += 1
        event_num = len(event_dict)
        feature_dict = {}
        label_dict = {}
        cohort_dict = {}
        random.shuffle(self.all_data)
        for sample in self.all_data:
            if sample.disease == "CRC":
                label_dict[sample.ID] = 1
            if sample.disease == "control":
                label_dict[sample.ID] = 0
            cohort_dict[sample.ID] = sample.cohort
            feature_dict[sample.ID] = np.zeros(event_num)
            for bkp in sample.filter_bkps:
                edge = self.select_tag(bkp)
                if edge in event_dict:
                    feature_dict[sample.ID][event_dict[edge]] = 1

        feature_data = list(feature_dict.values())
        train_label = list(label_dict.values())
        sel = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5,random_state = np.random.seed(2021))
        sel.fit(feature_data, train_label)
        model = SelectFromModel(sel, prefit=True, threshold=-np.inf, max_features = 200)
        new_feature_data = model.transform(feature_data)
        print (len(new_feature_data[0]))

        clf = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5,random_state = np.random.seed(2021)) 
        # clf.fit(train_data, train_label)  
        results = cross_val_score(clf, new_feature_data, train_label, cv=5)
        print (np.mean(results), results)

        # return feature_dict, label_dict, cohort_dict


    def select_tag(self, bkp):
        if self.level == 7:
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
        else:
            from_tax = bkp.from_ref_lineage.split(";")[self.level] #  bkp.from_ref
            to_tax = bkp.to_ref_lineage.split(";")[self.level]  #bkp.to_ref
            tax = sorted([from_tax, to_tax])
            tag = "&".join(tax)
            return tag

    def select_ranking(self, cohort_data):
        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0

        for sample in cohort_data:
            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == "CRC":
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == "control":
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
            sample_dict = {}
            for bkp in sample.filter_bkps:
                edge = self.select_tag(bkp)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                if edge not in specific_HGT:
                    specific_HGT[edge] = [0, 0]
                if edge not in cohort_sam_HGT[sample.cohort]:
                    cohort_sam_HGT[sample.cohort][edge] = [0, 0]
                if sample.disease == "CRC":
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == "control":
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 


        i = 0
        select_nodes = {}
        select_edges = {}
        for tag in specific_HGT:
            if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) <= self.common_ratio:
                continue
            if tag not in select_edges:
                select_edges[tag] = i
                i += 1
            array = tag.split("&")
            if self.level == 7:
                node1 = array[0] + "&" + array[2]
                node2 = array[1] + "&" + array[3]
            else:
                node1 = array[0]
                node2 = array[1]
            select_nodes[node1] = 1
            select_nodes[node2] = 1
        print ( "from %s choose %s edges."%(len(specific_HGT), len(select_edges)))

        for sample in cohort_data:
            sample.given_nodes_make_matrix(select_nodes, self.window, select_edges)


        train_data, train_label =  [], []
        for sample in cohort_data:
            sample_array = sample.select_feature_array
            train_data.append(sample_array)
            train_label.append(sample.disease)
        # self.complex_data(cohort_data)
        print ("start ranking")
        clf = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5,random_state = np.random.seed(2021)) 
        clf.fit(train_data, train_label)  
        np.save('importances.npy',clf.feature_importances_) 
        with open('features.pkl', 'wb') as f:
            pickle.dump(select_edges, f)
        print ("feature ranking score saved")

    def get_top_feature(self, select_feature_num):
        feature_importances = np.load('importances.npy', allow_pickle=True) 
        with open('features.pkl', 'rb') as f:
            select_edges = pickle.load(f)

        i = 0
        feature_importances_dict = {}
        for edge in select_edges.keys():
            feature_importances_dict[edge] = feature_importances[i]
            i += 1
        sort_select_edges = sorted(feature_importances_dict.items(), key=lambda item: item[1], reverse = True)
        print (sort_select_edges[:5])

        select_nodes = {}
        select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            select_edges[edge] = i
            array = edge.split("&")
            node1 = array[0]
            node2 = array[1]
            # if node1 in marker_genus or node2 in marker_genus:
            #     print ("same feature", i, node1, node2)
            select_nodes[node1] = 1
            select_nodes[node2] = 1
        # print (select_nodes, select_edges)
        print ("features readed.")
        return select_nodes, select_edges

    def complex_feature(self, select_feature_num):
        auc_list = []
        datasets = list(self.diff_cohorts.keys())

        # train_sam, test_sam = self.split_test_train(0)
        # all_sam = train_sam + test_sam   
        # self.select_ranking(all_sam)  

        for lack in range(len(self.diff_cohorts)):
            train_sam, test_sam = self.split_test_train(lack)
            print ("samples splitted.")
            self.select_ranking(train_sam)
            features_dict, select_edges = self.get_top_feature(select_feature_num)
            print ("%s features selected."%(len(select_edges)))
            self.get_event_matrix(select_edges)
            print (self.event_array)
            check_features = {}
            for sample in self.all_data:
                sample.given_nodes_make_matrix(features_dict, self.window, select_edges)
            print ("prepare data")
            train_data, train_label =  self.complex_data(train_sam)
            test_data, test_label = self.complex_data(test_sam) 


            all_data = list(train_data) + list(test_data)
            pca = PCA(n_components = 100)
            new_data = pca.fit_transform(all_data)
            print ("explained_variance_ratio_:", pca.explained_variance_ratio_)
            train_data = new_data[:len(train_data)]
            test_data = new_data[len(train_data):]



            if len(test_data) < 10:
                continue
            print ("start training, used feature:", len(train_data[0]))
            clf = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5, random_state = np.random.seed(2021)) #max_depth=2, random_state=0 entropy gini
            clf.fit(train_data, train_label)     
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            print (datasets[lack] , len(self.diff_cohorts[datasets[lack]]), "AUC", roc_auc) 
            print ("*************************")
            auc_list.append(roc_auc)
        return auc_list

    def get_event_matrix(self, select_edges):
        event_array = np.zeros((len(select_edges), len(select_edges)))
        events = list(select_edges.keys())
        for i in range(len(events)):
            for j in range(len(events)):
                event_1 = events[i]
                event_2 = events[j]
                tag_array_1 = event_1.split("&")
                tag_array_2 = event_2.split("&")
                for tag in tag_array_1:
                    if tag in tag_array_2:
                        event_array[i][j] = 1
        # print (event_array)
        # print (list(event_array))
        self.event_array = event_array
        # return event_array

    def complex_data(self, cohort_data):
        data = []
        label = []
        for sample in cohort_data:
            # density = nx.density(sample.select_feature_graph)
            # transitivity = nx.transitivity(sample.select_feature_graph)
            # sample_array = np.array(sample.select_feature_matrix).flatten()   #[density, transitivity]
            # sample_array = [density, transitivity]
            # sample_array = sample.marker_abundance + sample.basic_features + [sample.bases] + list(sample.select_feature_array)
            sample_array = sample.select_feature_array
            # sample_array = sample.marker_abundance
            # sample_array = self.event_array * sample.select_feature_array
            # sample_array = sample_array.flatten()
            # sample_array = linalg.eigvals(sample_array)
            # sample_array, b = np.linalg.eigh(sample_array)
            # print (sample_array)
            data.append(sample_array)
            # print (sample_array, sample.disease)
            label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

    def __del__(self):
        del self.all_data
        # del self.diff_cohorts
        print ("object, deleted")

class RF():

    def __init__(self, level, common_ratio):

        # analyze = Analyze(level)
        # self.all_data = analyze.data
        # del analyze
        # with open("sample_data", "wb") as fp:
        #     pickle.dump(self.all_data, fp)

        with open("sample_data", "rb") as fp:
            self.all_data = pickle.load(fp)
        self.diff_cohorts = self.classify_cohorts()        
        self.window = 500
        self.common_ratio = common_ratio #0.01
        self.level = self.all_data[0].level
        self.remove_adenoma = True
        self.event_array = ''
        print ("RF init done")

    def classify_cohorts(self):
        diff_cohorts = {}
        for sample in self.all_data:
            if sample.cohort not in diff_cohorts:
                diff_cohorts[sample.cohort] = [sample]
            else:
                diff_cohorts[sample.cohort].append(sample)
        return diff_cohorts

    def select_tag(self, bkp):
        if self.level == 7:
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
        else:
            from_tax = bkp.from_ref_lineage.split(";")[self.level] #  bkp.from_ref
            to_tax = bkp.to_ref_lineage.split(";")[self.level]  #bkp.to_ref
            tax = sorted([from_tax, to_tax])
            tag = "&".join(tax)
            return tag

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

    def remove_rare(self, specific_HGT):
        new_specific_HGT = {}
        for tag in specific_HGT:
            flag = True
            num = 0
            for cohort in cohort_sam_num.keys():
                if tag not in cohort_sam_HGT[cohort]:
                    continue
                crc_ratio = float(cohort_sam_HGT[cohort][tag][0]) / cohort_sam_num[cohort][0]
                control_ratio = float(cohort_sam_HGT[cohort][tag][1]) / cohort_sam_num[cohort][1]
                if max([crc_ratio, control_ratio]) < self.common_ratio:
                    pass
                else:
                    num += 1
            if num >= 2:
                new_specific_HGT[tag] = 1
        specific_HGT = new_specific_HGT
        return specific_HGT

    def select_ranking(self, cohort_data):
        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0

        for sample in cohort_data:
            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == "CRC":
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == "control":
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
            sample_dict = {}
            for bkp in sample.filter_bkps:
                edge = self.select_tag(bkp)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                if edge not in specific_HGT:
                    specific_HGT[edge] = [0, 0]
                if edge not in cohort_sam_HGT[sample.cohort]:
                    cohort_sam_HGT[sample.cohort][edge] = [0, 0]
                if sample.disease == "CRC":
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == "control":
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 


        i = 0
        select_nodes = {}
        select_edges = {}
        for tag in specific_HGT:
            if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) <= self.common_ratio:
                continue
            if tag not in select_edges:
                select_edges[tag] = i
                i += 1
            array = tag.split("&")
            if self.level == 7:
                node1 = array[0] + "&" + array[2]
                node2 = array[1] + "&" + array[3]
            else:
                node1 = array[0]
                node2 = array[1]
            select_nodes[node1] = 1
            select_nodes[node2] = 1
        print ( "from %s choose %s edges."%(len(specific_HGT), len(select_edges)))

        for sample in cohort_data:
            sample.given_nodes_make_matrix(select_nodes, self.window, select_edges)


        train_data, train_label =  [], []
        for sample in cohort_data:
            sample_array = sample.select_feature_array
            train_data.append(sample_array)
            train_label.append(sample.disease)
        # self.complex_data(cohort_data)
        print ("start ranking")
        clf = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5,random_state = np.random.seed(2021)) 
        clf.fit(train_data, train_label)  
        np.save('importances.npy',clf.feature_importances_) 
        with open('features.pkl', 'wb') as f:
            pickle.dump(select_edges, f)
        print ("feature ranking score saved")

    def select_ranking_ave(self, test_cohort):
        specific_HGT = {} 
        crc_num = 0
        control_num = 0

        for sample in self.all_data:
            if sample.disease == "CRC":
                crc_num += 1
            if sample.disease == "control":
                control_num += 1
            sample_dict = {}
            for bkp in sample.filter_bkps:
                edge = self.select_tag(bkp)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                if edge not in specific_HGT:
                    specific_HGT[edge] = [0, 0]
                if sample.disease == "CRC":
                    specific_HGT[edge][0] += 1
                if sample.disease == "control":
                    specific_HGT[edge][1] += 1

        i = 0
        select_nodes = {}
        select_edges = {}
        for tag in specific_HGT:
            if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) <= self.common_ratio:
                continue
            if tag not in select_edges:
                select_edges[tag] = i
                i += 1
            array = tag.split("&")
            node1 = array[0]
            node2 = array[1]

            select_nodes[node1] = 1
            select_nodes[node2] = 1
        print ( "from %s choose %s edges."%(len(specific_HGT), len(select_edges)))

        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_nodes, self.window, select_edges)

        scores_list = []
        feature_importances_scores = np.zeros(len(select_edges))

        for cohort in self.diff_cohorts:
            if cohort == test_cohort:
                continue
            train_data, train_label =  [], []
            for sample in self.diff_cohorts[cohort]:
                sample_array = sample.select_feature_array
                train_data.append(sample_array)
                # train_label.append(sample.disease)
                if sample.disease == "CRC":
                    train_label.append(1)
                if sample.disease == "control":
                    train_label.append(0)   

            print ("start ranking for %s"%(cohort))
            clf = RandomForestClassifier(n_estimators=1000, criterion="entropy", min_samples_leaf=5,random_state = np.random.seed(2021)) 
            clf.fit(train_data, train_label)  
            # feature_importances_scores += clf.feature_importances_
            scores_list.append(clf.feature_importances_)
        scores_list = np.array(scores_list).T
        for i in range(len(select_edges)):
            feature_importances_scores[i] = np.median(scores_list[i])
            # feature_importances_scores[i] = min(scores_list[i])
        np.save('importances.npy',feature_importances_scores) 
        with open('features.pkl', 'wb') as f:
            pickle.dump(select_edges, f)
        print ("feature ranking score saved")

    def get_top_feature(self, select_feature_num):
        feature_importances = np.load('importances.npy', allow_pickle=True) 
        with open('features.pkl', 'rb') as f:
            select_edges = pickle.load(f)

        i = 0
        feature_importances_dict = {}
        for edge in select_edges.keys():
            feature_importances_dict[edge] = feature_importances[i]
            i += 1
        sort_select_edges = sorted(feature_importances_dict.items(), key=lambda item: item[1], reverse = True)
        print (sort_select_edges[:5])

        select_nodes = {}
        select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            select_edges[edge] = i
            array = edge.split("&")
            node1 = array[0]
            node2 = array[1]
            # if node1 in marker_genus or node2 in marker_genus:
            #     print ("same feature", i, node1, node2)
            select_nodes[node1] = 1
            select_nodes[node2] = 1
        # print (select_nodes, select_edges)
        print ("features readed.")
        return select_nodes, select_edges

    def LODO(self, select_feature_num):
        auc_list = []
        datasets = list(self.diff_cohorts.keys())

        train_sam, test_sam = self.split_test_train(0)
        # all_sam = train_sam + test_sam   
        # self.select_ranking(all_sam)  

        for lack in range(len(self.diff_cohorts)):
            train_sam, test_sam = self.split_test_train(lack)
            self.select_ranking(train_sam)
            # self.select_ranking_ave(datasets[lack])
            features_dict, select_edges = self.get_top_feature(select_feature_num)
            print ("%s HGT events were selected."%(len(select_edges)))
            self.get_event_matrix(select_edges)
            # self.get_species_matrix(features_dict, select_edges)

            check_features = {}
            for sample in self.all_data:
                sample.given_nodes_make_matrix(features_dict, self.window, select_edges)
                # sample.get_species_node_feature(features_dict, self.window, select_edges)
            print ("prepare data")
            train_data, train_label =  self.complex_data(train_sam)
            test_data, test_label = self.complex_data(test_sam) 
            print ("start training, used feature:", len(train_data[0]))

            clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy",\
             min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
            #max_depth=2, random_state=0 entropy gini #min_samples_leaf=5
            clf.fit(train_data, train_label)     
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])

            # clf = Lasso(alpha=0.1)
            # clf.fit(train_data, train_label)     
            # roc_auc = roc_auc_score(test_label, clf.predict(test_data))

            print (datasets[lack] , len(self.diff_cohorts[datasets[lack]]), "AUC", roc_auc) 
            print ("*************************")
            auc_list.append(roc_auc)
        return auc_list

    def cross_validation(self, select_feature_num):
        auc_list = []
        datasets = list(self.diff_cohorts.keys())
        random.Random(2022).shuffle(self.all_data)
        fold_num = 5
        each_fold = int(len(self.all_data)/fold_num)


        # train_sam, test_sam = self.split_test_train(0)
        # all_sam = train_sam + test_sam   
        # self.select_ranking(all_sam)  

        for f in range(fold_num):
            train_sam = self.all_data[:each_fold*f] + self.all_data[each_fold*f+each_fold:]
            test_sam = self.all_data[each_fold*f:each_fold*f+each_fold]
            self.select_ranking(train_sam)
            # self.select_ranking_ave(datasets[lack])
            features_dict, select_edges = self.get_top_feature(select_feature_num)
            print ("%s HGT events were selected."%(len(select_edges)))
            self.get_event_matrix(select_edges)
            for sample in self.all_data:
                sample.given_nodes_make_matrix(features_dict, self.window, select_edges)
                # sample.get_species_node_feature(features_dict, self.window, select_edges)
            print ("prepare data")
            train_data, train_label =  self.complex_data(train_sam)
            test_data, test_label = self.complex_data(test_sam) 
            print ("start training, used feature:", len(train_data[0]))
            # clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy",\
            #  min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
            # #max_depth=2, random_state=0 entropy gini #min_samples_leaf=5
            # clf.fit(train_data, train_label)     
            # roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])


            clf = Ridge(alpha=1.0)
            clf.fit(train_data, train_label)     
            roc_auc = roc_auc_score(test_label, clf.predict(test_data))
            print (f, each_fold*f, each_fold*f+each_fold , "AUC", roc_auc) 
            print ("*************************")
            auc_list.append(roc_auc)
        return auc_list

    def get_event_matrix(self, select_edges):
        event_array = np.zeros((len(select_edges), len(select_edges)))
        events = list(select_edges.keys())
        # print ("events:", events[:5])
        for i in range(len(events)):
            for j in range(len(events)):
                event_1 = events[i]
                event_2 = events[j]
                tag_array_1 = event_1.split("&")
                tag_array_2 = event_2.split("&")
                if i == j:
                    continue
                for tag in tag_array_1:
                    if tag in tag_array_2:
                        event_array[i][j] = 1
                        
        laplacian = unnormalized_laplacian(event_array)
        a, b = np.linalg.eig(laplacian)
        self.event_array = b.T.real

    def get_species_matrix(self, features_dict, select_edges):
        species_array = np.zeros((len(features_dict), len(features_dict)))
        species = list(features_dict.keys())

        for i in range(len(species)):
            for j in range(len(species)):
                event_1 = species[i]
                event_2 = species[j]
                tag = event_1 + "&" + event_2
                if i == j:
                    continue
                for tag in select_edges:
                    species_array[i][j] = 1
                        
        laplacian = unnormalized_laplacian(species_array)
        a, b = np.linalg.eig(laplacian)
        eig_vector = b.T.real
        self.event_array = eig_vector

    def convert_GCN_data(self, sample_array):
        edge_list = []
        for i in range(len(sample_array)):
            for j in range(i, len(sample_array)):
                if sample_array[i][j] == 1:
                    edge_list.append([i, j])
                    if i != j:
                        edge_list.append([j, i])
        return edge_list

    def get_GCN_data(self, cohort_data):
        data = []
        label = []
        node_feature = []
        for sample in cohort_data:
            node_feature.append(list(sample.select_feature_array))
            sample_array = self.event_array * sample.select_feature_array
            # sample_array = svd_denoise(sample_array, percent_cut)
            edge_list = self.convert_GCN_data(sample_array)
            data.append(edge_list)
            if sample.disease == "CRC":
                label.append(1)
            else:
                label.append(0)
        # data = np.array(data)
        # label = np.array(label)
        return data, label, node_feature

    def prepare_GCN(self, select_feature_num):
        auc_list = []
        datasets = list(self.diff_cohorts.keys())

        train_sam, test_sam = self.split_test_train(0)
        all_sam = train_sam + test_sam   
        self.select_ranking(all_sam)  
        features_dict, select_edges = self.get_top_feature(select_feature_num)
        self.get_event_matrix(select_edges)

        for sample in self.all_data:
            sample.given_nodes_make_matrix(features_dict, self.window, select_edges) 


        train_data, train_label, node_feature =  self.get_GCN_data(train_sam)    
        print (train_data[0])
        print (train_label[0])   
        np.save('training_graph.npy', train_data) 
        np.save('training_label.npy', train_label) 
        np.save('training_nodes.npy', node_feature)

    def prepare_NN(self, select_feature_num):
        datasets = list(self.diff_cohorts.keys())

        train_sam, test_sam = self.split_test_train(0)
        all_sam = train_sam + test_sam   
        self.select_ranking(all_sam)  

        features_dict, select_edges = self.get_top_feature(select_feature_num)
        print ("%s features selected."%(len(select_edges)))
        self.get_event_matrix(select_edges)

        for sample in self.all_data:
            sample.given_nodes_make_matrix(features_dict, self.window, select_edges)

        feature_dict = {}
        lable_dict = {}
        for cohort in self.diff_cohorts.keys():
            corhort_sam = self.diff_cohorts[cohort]
            train_data, train_label =  self.complex_data(corhort_sam)
            feature_dict[cohort] = train_data
            lable_dict[cohort] = train_label

        with open("feature_dict.pkl", 'wb') as f:
            pickle.dump(feature_dict, f)
        with open("lable_dict.pkl", 'wb') as f:
            pickle.dump(lable_dict, f)

    def complex_data(self, cohort_data):
        data = []
        label = []
        for sample in cohort_data:
            # sample_array = sample.marker_abundance

            sample_array = np.dot(self.event_array, sample.select_feature_array)
            # sample_array = list(sample_array) + sample.marker_abundance

            data.append(sample_array)
            if sample.disease == "CRC":
                label.append(1)
            if sample.disease == "control":
                label.append(0)            
            # label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

    def __del__(self):
        del self.all_data
        del self.diff_cohorts
        print ("object, deleted")

def normalized_laplacian(adj_matrix):
    R = np.sum(adj_matrix, axis=1)
    R_sqrt = 1/np.sqrt(R)
    D_sqrt = np.diag(R_sqrt)
    I = np.eye(adj_matrix.shape[0])
    return I - D_sqrt * adj_matrix * D_sqrt

def unnormalized_laplacian(adj_matrix):
    R = np.sum(adj_matrix, axis=1)
    degreeMatrix = np.diag(R)
    return degreeMatrix - adj_matrix

def select_genus():
    marker_genus = {}
    for species in marker_species:
        genus = "g__" + species.split()[0]
        marker_genus[genus] = 1
    return marker_genus


if __name__ == "__main__":
    level_dict = {"phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
    gender_dict = {"male":0, "female":1, "nan": 2}
    # sra_meta = "italy.csv"
    pheno_file = "allmetadata.xlsx"#"CRC.xlsx"
    UHGG_meta = "/mnt/d/breakpoints/HGT/UHGG/genomes-all_metadata.tsv"
    cohort_abd = {"YuJ_2015":"2021-03-31.YuJ_2015.relative_abundance.xls",
    "WirbelJ_2018":"2021-03-31.WirbelJ_2018.relative_abundance.xls",
    "HanniganGD_2017":"2021-03-31.HanniganGD_2017.relative_abundance.xls",
    "YachidaS_2019":"2021-10-14.YachidaS_2019.relative_abundance.xls",
    "ThomasAM_2018a":"2021-03-31.ThomasAM_2018a.relative_abundance.xls",
    "ThomasAM_2018b":"2021-03-31.ThomasAM_2018b.relative_abundance.xls",
    "ZellerG_2014":"2021-03-31.ZellerG_2014.relative_abundance.xls",
    "FengQ_2015":"2021-03-31.FengQ_2015.relative_abundance.xls",
    "VogtmannE_2016":"2021-03-31.VogtmannE_2016.relative_abundance.xls"
    }
    marker_species = ["Peptostreptococcus stomatis", "Fusobacterium nucleatum", "Parvimonas spp.", "Porphyromonas asaccharolytica", "Gemella morbillorum",
    "Clostridium symbiosum", "Parvimonas micra", "Escherichia coli", "Streptococcus parasanguinis", "Clostridium leptum", "Clostridium hathewayi",
    "Anaerotruncus colihominis", "Prevotella copri", "Lachnospiraceae 3 1 57FAA CT1", "Actinomyces graevenitzii", "Alistipes spp."]
    marker_genus = select_genus()
    # print (len(marker_genus), marker_genus)
    # marker_species = ["Fusobacterium nucleatum", "Parvimonas micra", "Parvimonas spp.", "Gemella morbillorum", "Peptostreptococcus stomatis",\
    # "Solobacterium moorei", "Clostridium symbiosum", "Anaerococcus vaginalis", "Porphyromonas asaccharolytica", "Prevotella intermedia",\
    #  "Bacteroides fragilis", "Porphyromonas somerae", "Anaerococcus obesiensis", "Porphyromonas uenonis", "Peptostreptococcus anaerobius",\
    #  "Streptococcus constellatus", "Granulicatella adiacens"]

    sample_abd = get_samples_abd()
    phenotype = Phenotype()
    taxonomy = Taxonomy()
    
    result_file = open("random_forest.log", 'w')
    # percent_cut = 0.9
    TREE_NUM = 1000
    LEAF_S = 5
    pheno_group_1 = "control" 
    pheno_group_2 = "CRC" 
    level, common_ratio, feature_num = 5, 0.1, 1000   #0.06, 150


    rf = RF(level, common_ratio)
    # rf.prepare_NN(feature_num)
    # auc_list = rf.cross_validation(feature_num)
    auc_list = rf.LODO(feature_num)
    print (common_ratio, feature_num, auc_list, np.mean(auc_list), np.median(auc_list))

    # for common_ratio in [0.02, 0.06, 0.1, 0.14]:
    #     # common_ratio *= 0.01
    #     rf = RF(level, common_ratio)
    #     for feature_num in [150, 500, 1000]:
    #         # auc_list = rf.cross_validation(feature_num)
    #         auc_list = rf.LODO(feature_num)
    #         print (common_ratio, feature_num, auc_list, np.mean(auc_list), np.median(auc_list))
    #         print (common_ratio, feature_num, auc_list, np.mean(auc_list), np.median(auc_list), file = result_file)
    #     del rf
    # result_file.close()


    # rf = RF(level, common_ratio)
    # for TREE_NUM in [500, 1000]:
    #     for LEAF_S in [5, 10, 20]:
    #         auc_list = rf.complex_feature(feature_num)
    #         print (TREE_NUM, LEAF_S, auc_list, np.mean(auc_list), np.median(auc_list))
    #         print (TREE_NUM, LEAF_S, auc_list, np.mean(auc_list), np.median(auc_list), file = result_file)
    #     # del rf
    # result_file.close()


