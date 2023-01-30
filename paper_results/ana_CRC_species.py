#!/usr/bin/env python3

"""
1. get differential HGTs between CRC and control.
2. establish the CRC classifier.

"""

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
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.feature_selection import SelectFromModel
from scipy.sparse import csgraph
from graph_denoise import svd_denoise
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.inspection import permutation_importance
from scipy.stats import wilcoxon
from sklearn.preprocessing import normalize
import scipy
from scipy.linalg import svd
import powerlaw
np.set_printoptions(threshold=sys.maxsize)
import scipy.special as sc
# from math import comb #The comb function is new in Python 3.8
from scipy.special import comb
# from deep_learning import Classifier
from KR_norm_juicer import KRnorm_sym

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
"VogtmannE_2016":"2021-03-31.VogtmannE_2016.relative_abundance.xls",
"KarlssonFH_2013":"2021-10-14.KarlssonFH_2013.relative_abundance.xls"
}
marker_species = ["Peptostreptococcus stomatis", "Fusobacterium nucleatum", "Parvimonas spp.", "Porphyromonas asaccharolytica", "Gemella morbillorum",
"Clostridium symbiosum", "Parvimonas micra", "Escherichia coli", "Streptococcus parasanguinis", "Clostridium leptum", "Clostridium hathewayi",
"Anaerotruncus colihominis", "Prevotella copri", "Eisenbergiella tayi", "Actinomyces graevenitzii", "Alistipes spp."]
# marker_species = ["Peptostreptococcus stomatis", "Fusobacterium nucleatum", "Parvimonas spp.", "Porphyromonas asaccharolytica", "Gemella morbillorum",
# "Clostridium symbiosum", "Parvimonas micra", "Escherichia coli", "Streptococcus parasanguinis", "Clostridium leptum", "Clostridium hathewayi",
# "Anaerotruncus colihominis", "Prevotella copri", "Actinomyces graevenitzii", "Alistipes spp."]   
# f = open("previous_16_abun_markers.csv", 'w')
# for sp in  marker_species:
#     print (sp, file = f)
# f.close()

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
        self.from_ref_lineage = taxonomy.taxonomy_dict[self.from_ref_genome]
        self.to_ref_genome = "_".join(self.to_ref.split("_")[:-1])
        self.to_ref_lineage = taxonomy.taxonomy_dict[self.to_ref_genome]
        self.score = float(list[11])
        # print (self.from_ref_genome, self.from_ref_lineage, self.to_ref_genome, self.to_ref_lineage)

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
    # print (marker_species_dict)
    return marker_species_dict

def get_genus_abd(marker_genus):
    sample_abd = {}
    marker_species_num = len(marker_genus)
    i = 0
    for genus in marker_genus:
        marker_genus[genus] = i
        i += 1
    # print (marker_species_num, marker_genus)
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
                genus_name = array[0].split("|")[-2]
                if species_name in marker_genus:
                    # print ("name", species_name, marker_genus)
                    species_index = marker_genus[species_name]
                    for j in range(len(sample_list)):
                        sample = sample_list[j]
                        abundance = float(array[j+1])
                        sample_abd[sample][species_index] += abundance
                elif genus_name in marker_genus:
                    # print ("name", species_name, marker_genus)
                    species_index = marker_genus[genus_name]
                    for j in range(len(sample_list)):
                        sample = sample_list[j]
                        abundance = float(array[j+1])
                        sample_abd[sample][species_index] += abundance
            i += 1
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
        self.ID_name = {}
        self.ID_marker_abundance = {}
        self.name_marker_abundance = {}
        self.read_pheno()
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/USA/usa.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/japan.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/yu_2015.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/germany.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/france/france.csv")
        self.read_sra_meta("/mnt/d/breakpoints/script/analysis/italy.csv")
        # self.read_sra_meta("/mnt/d/breakpoints/script/analysis/new_result/usa_canada.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/austria/austria.csv")   
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/t2d.csv")  
        
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
            if sra_meta == "/mnt/d/breakpoints/HGT/CRC/t2d.csv":
                # sample_name = df["Sample Name"][i]
                sample_name = sra_ID
            # print (sra_ID, sample_name)
            study = df["SRA Study"][i]
            # for col in df.columns:
            #     print (col, df[col][i])
            
            if sample_name not in self.name_disease:
                continue
            self.ID_name[sra_ID] = sample_name
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
            if cohort == "KarlssonFH_2013":
                # for z in range(len(row)):
                #     print (z, row[z])
                sample_name = row[21]
                if len(sample_name.split(";")) > 1:
                    sample_name = sample_name.split(";")[0]
                # print (cohort, sample_name)
            # if cohort == "VogtmannE_2016":
            #     print (cohort, sample_name, condition, row[7])

            self.name_disease[sample_name] = condition
            self.name_cohort[sample_name] = cohort
            self.name_basics[sample_name] = [age, gender, BMI]
            self.name_bases[sample_name] = bases
            self.name_sample_id[sample_name] = sample_id
            # print (sample_name,condition) 
        # print (crc, control)

def get_tag(bkp, level):
    if level == 7:
        from_tax = bkp.from_ref_genome
        to_tax = bkp.to_ref_genome
    else:
        from_tax = bkp.from_ref_lineage.split(";")[level] #  bkp.from_ref
        to_tax = bkp.to_ref_lineage.split(";")[level]  #bkp.to_ref
    from_tax = "_".join(from_tax.split())
    to_tax = "_".join(to_tax.split())
    tax = sorted([from_tax, to_tax])
    new_tag = "&".join(tax)       
    node1 = from_tax  
    node2 = to_tax  
    return new_tag

class Sample():

    def __init__(self, bkp_file, ID):
        self.bkp_file = bkp_file
        self.bkps = []
        self.ID = ID
        self.tag = "normal"
        self.genus_abundance = ''
        self.reads_num = 0
        self.level = 5
        # self.name = phenotype.ID_name[ID]
        
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
            # print ("## no pheno", ID) 
            self.tag = "no pheno"    
        self.read_bkp()
        if len(self.bkps) == 0:
            self.tag = "no bkp"
            # print ("no bkp", ID)
        if self.tag != "no bkp":
            self.select_feature_array = []

    def read_bkp(self):
        f = open(self.bkp_file)
        all_rows = csv.reader(f)
        for row in all_rows:
            # print (row)
            if row[0][0] == "#":
                self.reads_num = int(row[0].split(";")[0].split(":")[1])
                # print (row, self.reads_num)
                pass
            elif row[0] == "from_ref":
                pass
            else:
                if self.reads_num == 0:
                    print ("old bkp", self.bkp_file)
                    break
                eb = Acc_Bkp(row)
                if eb.from_ref_genome != eb.to_ref_genome:
                    self.bkps.append(eb)
        f.close()

    def given_nodes_make_matrix(self, select_edges):
        common_nodes_dict = {}
        i = 0
        for edge in select_edges:
            array = edge.split("&")
            if array[0] not in common_nodes_dict:
                common_nodes_dict[array[0]] = i
                i += 1
            elif array[1] not in common_nodes_dict:
                common_nodes_dict[array[1]] = i
                i += 1
        nodes = list(common_nodes_dict.keys())
        n = len(common_nodes_dict)
        self.select_feature_array = np.zeros(len(select_edges))
        nodes_index = {}
        for i in range(len(nodes)):
            nodes_index[nodes[i]] = i
        for bkp in self.bkps:
            new_tag = get_tag(bkp, self.level)
            if new_tag in select_edges:
                self.select_feature_array[select_edges[new_tag]] = 1 #1#bkp.cross_split_reads/self.reads_num #1 

    def get_HGT_matrix(self, level, edge_num):
        nodes_index = {}
        bkp_score = {}
        choose_edge = {}
        for bkp in self.bkps:
            edge = get_tag(bkp, level)
            support_ratio = bkp.cross_split_reads/self.reads_num
            if edge not in bkp_score:
                bkp_score[edge] = support_ratio
            if support_ratio > bkp_score[edge]:
                bkp_score[edge] = support_ratio
        sort_bkp_score = sorted(bkp_score.items(), key=lambda item: item[1], reverse = True)
        choose_edge = {}
        # min_score = 1
        # if len(sort_bkp_score) > 50:
        #     min_score = sort_bkp_score[49][1]
        total_edge_num = len(sort_bkp_score)
        # edge_num = 50
        if total_edge_num < edge_num:
            edge_num = len(sort_bkp_score)
        for i in range(edge_num):
            choose_edge[sort_bkp_score[i][0]] = 1


        i = 0
        for bkp in self.bkps:
            # print ()
            edge = get_tag(bkp, level)
            if edge not in choose_edge:
                continue
            # support_ratio = bkp.cross_split_reads/self.reads_num
            # if support_ratio < min_score:
            #     continue
            array = edge.split("&")
            node1 = array[0]
            node2 = array[1]
            if node1 not in nodes_index:
                nodes_index[node1] = i
                i += 1
            if node2 not in nodes_index:
                nodes_index[node2] = i
                i += 1

        HGT_matrix = np.zeros((len(nodes_index), len(nodes_index)))       
        for bkp in self.bkps:
            edge = get_tag(bkp, level)
            # support_ratio = bkp.cross_split_reads/self.reads_num
            # if support_ratio < min_score:
            #     continue
            if edge not in choose_edge:
                continue
            array = edge.split("&")
            node1_ind = nodes_index[array[0]]
            node2_ind = nodes_index[array[1]]

            HGT_matrix[node1_ind][node2_ind] = 1
            HGT_matrix[node2_ind][node1_ind] = 1
        # origin_matrix = HGT_matrix.copy()
        return HGT_matrix, total_edge_num

    def get_network_properties(self, level, edge_num):
        HGT_matrix, total_edge_num =  self.get_HGT_matrix(level, edge_num)
        HGT_matrix = nx.from_numpy_matrix(HGT_matrix)

        density = nx.density(HGT_matrix)
        transitivity = nx.transitivity(HGT_matrix)
        algebraic_connectivity = nx.algebraic_connectivity(HGT_matrix)
        # algebraic_connectivity = 0
        assortativity = nx.degree_assortativity_coefficient(HGT_matrix)
        # assortativity = 0
        node_num = HGT_matrix.number_of_nodes()
        edge_num = HGT_matrix.number_of_edges()

        return [round(density,3), round(transitivity,3), round(algebraic_connectivity,3),\
         round(assortativity,3), int(node_num), int(edge_num)], total_edge_num

    def get_network_features(self, level):
        HGT_matrix =  self.get_HGT_matrix(level)
        HGT_matrix = nx.from_numpy_matrix(HGT_matrix)

        # density = nx.density(HGT_matrix)
        # transitivity = nx.transitivity(HGT_matrix)
        # algebraic_connectivity = nx.algebraic_connectivity(HGT_matrix)
        assortativity = nx.degree_assortativity_coefficient(HGT_matrix)
        node_num = HGT_matrix.number_of_nodes()
        # edge_num = HGT_matrix.number_of_edges()

        return [assortativity/self.bases, node_num/self.bases]

    def judge_scale_free(self, level, edge_num):
        HGT_matrix, total_edge_num =  self.get_HGT_matrix(level, edge_num)
        HGT_matrix = nx.from_numpy_matrix(HGT_matrix)
        p1, p2, p3 = infer_scale_free(HGT_matrix)
        return p1, p2, p3, total_edge_num

    def prepare_network_vector(self):
        
        """
        level = 2
        geno_num = 100
        # # find class with maximum degree
        nodes_index = {}
        i = 0
        for bkp in self.bkps:
            edge = get_tag(bkp, level)
            array = edge.split("&")
            node1 = array[0]
            node2 = array[1]
            if node1 not in nodes_index:
                nodes_index[node1] = i
                i += 1
            if node2 not in nodes_index:
                nodes_index[node2] = i
                i += 1
        HGT_matrix = np.zeros((len(nodes_index), len(nodes_index)))
        for bkp in self.bkps:
            edge = get_tag(bkp, level)
            edge = get_tag(bkp, level)
            array = edge.split("&")
            node1_ind = nodes_index[array[0]]
            node2_ind = nodes_index[array[1]]
            HGT_matrix[node1_ind][node2_ind] = 1
            HGT_matrix[node2_ind][node1_ind] = 1
        R = np.sum(HGT_matrix, axis=1)
        degree_dict = {}
        i = 0
        for key in nodes_index.keys():
            degree_dict[key] = R[i]
            i += 1
        sort_degree_dict = sorted(degree_dict.items(), key=lambda item: item[1], reverse = True)
        # print(len(sort_degree_dict))
        class_dict = {}
        if geno_num > len(sort_degree_dict):
            num = len(sort_degree_dict)
        else:
            num = geno_num
        for i in range(num):
            class_dict[sort_degree_dict[i][0]] = i
        """


        HGT_matrix = np.zeros((len(class_dict), len(class_dict)))
        for bkp in self.bkps:
            edge = get_tag(bkp, level)
            array = edge.split("&")
            node1 = array[0]
            node2 = array[1]
            if node1 == node2:
                continue
            if node1 not in class_dict or node2 not in class_dict:
                continue
            HGT_matrix[class_dict[node1]][class_dict[node2]] = 1
            HGT_matrix[class_dict[node2]][class_dict[node1]] = 1

        # C = np.ones((len(class_dict), len(class_dict)))
        # HGT_matrix = HGT_matrix + pseudo * C
        # HGT_matrix = svd_denoise(HGT_matrix, svd_num)

        if select_m == "density":
            HGT_matrix = nx.from_numpy_matrix(HGT_matrix)
            density = nx.density(HGT_matrix)        
            return [density]  

        elif select_m == "laplacian":
            HGT_matrix = nx.from_numpy_matrix(HGT_matrix)
            # laplacian = nx.normalized_laplacian_matrix(HGT_matrix)
            # laplacian = scipy.sparse.csr_matrix.toarray(laplacian)
            # a, b = np.linalg.eig(laplacian)
            # f = sorted(a.real)[1]
            f = nx.algebraic_connectivity(HGT_matrix)
            return [f]

        elif select_m == "eigenvalue":
            a, b = np.linalg.eig(HGT_matrix)
            a = a.real
            a = sorted(a, reverse = True)
            return a
        
        elif select_m == "degree":
            R = np.sum(HGT_matrix, axis=1)
            return R

def svd_denoise(A, num):
    U, s, Vh = svd(A)
    new_A = np.dot(U[:,0:num], np.dot(np.diag(s[0:num]), Vh[0:num,:]))
    return new_A

def egen_denoise(A, num):
    vals, vecs = np.linalg.eig(A)
    new_A = np.dot(vecs[:,0:num], np.dot(np.diag(vals[0:num]), vecs.T[0:num,:]))
    return new_A

def infer_scale_free(g):
    d= g.degree()
    d = [x[1] for x in list(d)]
    # print (d)
    result = powerlaw.Fit(d)
    R1, p1 = result.distribution_compare('power_law', 'lognormal_positive')

    # result = powerlaw.Fit(d)
    # p2=result.distribution_compare('power_law', 'lognormal')

    result = powerlaw.Fit(d)
    R2, p2 = result.distribution_compare('power_law', 'exponential', normalized_ratio=True)

    result = powerlaw.Fit(d)
    R3, p3 = result.distribution_compare('power_law', 'stretched_exponential',normalized_ratio=True)  #Weibull
    return R1, R2, R3

class RF():

    def __init__(self):
        self.all_data = []
        self.disease_sample_num_cohort = {}
        # self.all_samples()
        # with open("sample_data", "wb") as fp:
        #     pickle.dump(self.all_data, fp)
        with open("sample_data", "rb") as fp:
            self.all_data = pickle.load(fp)
        self.cohorts_names = {}        
        self.HGT_network = ''
        self.species_network = ''
        self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b",\
            "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        # self.cohort_names = ["VogtmannE_2016","WirbelJ_2018","ZellerG_2014","FengQ_2015", "ThomasAM_2018a","ThomasAM_2018b"]
        print ("RF init done")

    def combine_markers(self, select_feature_num):
        # with open('selected_26_edges.pkl', 'rb') as f:
        #     select_HGTs = pickle.load(f)
        #     print ("origin:", select_HGTs)
        out_marker = open("CRC_prediction_markers.csv", 'w')
        select_HGTs, abun_related_HGTs = self.select_top_HGT(select_feature_num)
        # print (select_HGTs, len(select_HGTs))
        # with open('filtered_HGTs.pkl', 'rb') as f:
        #     abun_related_HGTs = pickle.load(f)
        # select_edges = {**select_HGTs, **abun_related_HGTs} # combine abundance-related HGT events
        select_edges = select_HGTs
        with open('selected_diff_edges.pkl', 'wb') as f:
            pickle.dump(select_edges, f) 
        # select_edges = select_HGTs
        i = 0
        for edge in select_edges:
            select_edges[edge] = i
            print (edge.split("&")[0], edge.split("&")[1], sep = ",", file = out_marker)
            i += 1
        # print ("final HGT num:", len(select_edges))
        self.get_HGT_network(select_edges)

        for edge in select_HGTs:  # add the species in HGT to abundance markers
            array = edge.split("&")
            if array[0] not in marker_genus:
                marker_genus[array[0]] = 1
            if array[1] not in marker_genus:
                marker_genus[array[1]] = 1
        i = 0
        for genus in marker_genus:
            marker_genus[genus] = i
            i += 1
        # print (len(marker_genus), "species for abundance", marker_genus)
        self.get_abundance_marker_value(marker_genus, self.all_data)
        self.get_genus_network(marker_genus)
        for genus in marker_genus:
            print (genus, '', sep = ",", file = out_marker)

        self.store_markers(marker_genus, select_edges)
        print ("HGT num", len(select_edges), "abun num", len(marker_genus))
        return select_edges, marker_genus

    def get_abundance_marker_value(self, marker_genus, sample_list):

        sample_abd = get_genus_abd(marker_genus)
        # print (sample_abd)
        for sample in sample_list:
            if sample.ID not in phenotype.ID_name:
                sample_name = sample.ID
            else:
                sample_name = phenotype.ID_name[sample.ID]
            if sample_name not in sample_abd:
                sample_id = phenotype.name_sample_id[sample_name]
                genus_abundance = sample_abd[sample_id]
            else:
                genus_abundance = sample_abd[sample_name]
            sample.genus_abundance = genus_abundance

    def store_markers(self, marker_genus, select_edges):
        # store the markers
        with open('marker_genus.pkl', 'wb') as f:
            pickle.dump(marker_genus, f) 
        with open('select_edges.pkl', 'wb') as f:
            pickle.dump(select_edges, f) 
    
    def load_markers(self):
        print ("load saved markers.")
        with open('marker_genus.pkl', 'rb') as f:
            marker_genus = pickle.load(f)
        with open('select_edges.pkl', 'rb') as f:
            select_edges = pickle.load(f)
        self.get_HGT_network(select_edges)
        self.get_abundance_marker_value(marker_genus, self.all_data)
        self.get_genus_network(marker_genus)
        return marker_genus, select_edges

    def select_top_HGT_BK(self, select_feature_num):
        specific_HGT = {} 
        abun_related_HGTs = {}
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.all_data:
            if sample.disease == "CRC" :
                crc_num += 1
            if sample.disease == "control" :
                control_num += 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                if edge not in specific_HGT:
                    specific_HGT[edge] = [0, 0]
                if sample.disease == "CRC" :
                    specific_HGT[edge][0] += 1
                if sample.disease == "control" :
                    specific_HGT[edge][1] += 1    
        select_edges = {}
        # print (crc_num, control_num, crc_num+control_num, len(specific_HGT)/(crc_num+control_num))
        genus_level_markers = {}
        for marker in marker_genus:
            if marker[0] == "s":
                marker = "g__" + marker.split("_")[2]
            genus_level_markers[marker] = 1
        print ("abundance marker-related genus:", genus_level_markers)
        for tag in specific_HGT:
            if specific_HGT[tag][0] + specific_HGT[tag][1] < 25:
                continue
            array = tag.split("&")
            species_1 = array[0]
            species_2 = array[1]
            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
            U1, p = mannwhitneyu(crc_array, control_array)
            if p < 0.01 and species_1 in genus_level_markers and species_2 in genus_level_markers:
                abun_related_HGTs[tag] = 1
            select_edges[tag] = p 
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

        final_select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            final_select_edges[edge] = i
            if edge in abun_related_HGTs:
                del abun_related_HGTs[edge]
        print ("abun_related_HGTs", len(abun_related_HGTs))

        # tree_edges = {}
        # for i in range(100):
        #     tree_edges[sort_select_edges[i][0]] = 1
        # with open('selected_100_edges.pkl', 'wb') as f:
        #     pickle.dump(tree_edges, f) 
          
        return final_select_edges, abun_related_HGTs

    def select_top_HGT(self, select_feature_num):
        specific_HGT = {} 
        abun_related_HGTs = {}
        record_all_HGTs = {}
        crc_num = 0
        control_num = 0
        edge_distribution = []
        reads_num_list = []
        for sample in self.all_data:
            sample_dict = {}
            reads_num_list.append(sample.reads_num)
            min_split_num, p = get_split_reads_cutoff(sample.reads_num)
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                support_ratio = bkp.cross_split_reads/sample.reads_num
                # if support_ratio < min_cross:
                #     continue
                if support_ratio == 0:
                    continue
                # if bkp.cross_split_reads < min_split_num:
                #     continue
                # if bkp.pair_end/sample.reads_num < 2e-08:
                #     continue
                if edge not in record_all_HGTs:
                    record_all_HGTs[edge] = 1
                    specific_HGT[edge] = [[], []]
                if edge not in sample_dict:
                    sample_dict[edge] = 0
                # if support_ratio > sample_dict[edge]:
                #     sample_dict[edge] = support_ratio
                # sample_dict[edge] += support_ratio
                # sample_dict[edge] = 1
                sample_dict[edge] += bkp.cross_split_reads
            edge_distribution.append(len(sample_dict))
            # in each sample, choose same number of edge.       
            sample_dict = self.same_HGT_number(sample_dict, max_hgt_num, min_split_num)

            for edge in sample_dict:
                if sample.disease == "CRC" :
                    specific_HGT[edge][0].append(sample_dict[edge])
                else:
                    specific_HGT[edge][1].append(sample_dict[edge])
            if sample.disease == "CRC" :
                crc_num += 1
            if sample.disease == "control" :
                control_num += 1
        # print ("sample num :%s, Edge count: mean is %s, median is %s"%(len(edge_distribution),\
        #     np.mean(edge_distribution), np.median(edge_distribution)), np.std(edge_distribution))
        # print ("read num, median %s, mean %s"%(np.median(reads_num_list),np.mean(reads_num_list)), min(reads_num_list), max(reads_num_list))
        select_edges = {}
        # print (crc_num, control_num, crc_num+control_num, len(specific_HGT)/(crc_num+control_num))
        genus_level_markers = {}
        for marker in marker_genus:
            if marker[0] == "s":
                marker = "g__" + marker.split("_")[2]
            genus_level_markers[marker] = 1
        print ("abundance marker-related genus:", len(genus_level_markers))
        for tag in specific_HGT:
            if len(specific_HGT[tag][0]) + len(specific_HGT[tag][1]) < 25:
                continue
            array = tag.split("&")
            species_1 = array[0]
            species_2 = array[1]
            crc_array = specific_HGT[tag][0] + [0] * (crc_num - len(specific_HGT[tag][0]))
            control_array = specific_HGT[tag][1] + [0] * (control_num - len(specific_HGT[tag][1]))
            U1, p = mannwhitneyu(crc_array, control_array)
            if p < 0.01 and species_1 in genus_level_markers and species_2 in genus_level_markers:
                abun_related_HGTs[tag] = 1
            select_edges[tag] = p 
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

        final_select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            # if i < 3:
            #     print (sort_select_edges[i])
            final_select_edges[edge] = i
            if edge in abun_related_HGTs:
                del abun_related_HGTs[edge]
        print ("abun_related_HGTs", len(abun_related_HGTs))

        tree_edges = {}
        for i in range(len(sort_select_edges)): # store the all differential genome pairs
            # if sort_select_edges[i][1] > 0.05:
            #     break
            tag = sort_select_edges[i][0]
            crc_ratio = len(specific_HGT[tag][0])/crc_num
            control_ratio = len(specific_HGT[tag][1])/control_num
            tree_edges[sort_select_edges[i][0]] = [crc_ratio, control_ratio, sort_select_edges[i][1]]
        # with open('selected_diff_edges.pkl', 'wb') as f:
        #     pickle.dump(tree_edges, f) 
          
        return final_select_edges, abun_related_HGTs

    def same_HGT_number(self, sample_dict, choose_num, min_split_num): #only choose top x HGTs in each sample
        new_same_dict = {}
        sort_sample_dict = sorted(sample_dict.items(), key=lambda item: item[1], reverse = True)
        if len(sort_sample_dict) < choose_num:
            choose_num = len(sort_sample_dict)
        for z in range(choose_num):
            if sort_sample_dict[z][1] < min_split_num:
                continue
            # new_same_dict[sort_sample_dict[z][0]] = sort_sample_dict[z][1]
            new_same_dict[sort_sample_dict[z][0]] = 1
        # print (sort_sample_dict[choose_num-1])
        return new_same_dict

    def all_samples(self):
        all_acc_file = "acc.list"
        # result_dir = "new_result_more"
        result_dir = "new_result"
        os.system(f"ls {result_dir}/*acc.csv>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            ID = acc_file.split("/")[1].split(".")[0]
            new_file = acc_file.split("/")[1]
            acc_file = f"{result_dir}/{new_file}"
            
            sample = Sample(acc_file, ID)

            if sample.tag == "no pheno" or sample.tag == "no bkp":
                print (ID, sample.tag)
                continue  
            # if sample.cohort == "FengQ_2015":
            #     print (sample.cohort, ID, sample.disease)          
            if (sample.disease != "CRC") and (sample.disease != "control"):
                continue
            if str(sample.disease) == "nan":
                continue
            if sample.cohort == "HanniganGD_2017":
                continue
            # if sample.cohort == "WirbelJ_2018":
            #     continue
            #     print (sample.ID, phenotype.ID_name[ID], sample.disease)
            
            # if sample.cohort == "ThomasAM_2018a" or sample.cohort == "ThomasAM_2018b":
            #     sample.cohort = "ThomasAM_2018"
            # if sample.cohort in ["ThomasAM_2018", "VogtmannE_2016", "YachidaS_2019"]:
            #     continue
            # if sample.bases < 2000000000:
            #     continue
            # sample.check_graph()
            self.all_data.append(sample)
            if sample.cohort not in self.disease_sample_num_cohort:
                self.disease_sample_num_cohort[sample.cohort] = {"CRC":0, "control":0, "adenoma":0}
            self.disease_sample_num_cohort[sample.cohort][sample.disease] += 1

        for key in self.disease_sample_num_cohort:
            print ("Cohort", key, self.disease_sample_num_cohort[key])     

    def LODO(self, select_feature_num):
        select_edges,marker_genus = self.combine_markers(select_feature_num)
        # select_edges,marker_genus = self.load_markers()
        print ("HGT marker num:", len(select_edges))
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
        auc_list = []   
        auc_total = 0    
        for lack in range(len(self.cohort_names)):
            # print ("prepare data")
            train_data, train_label =  self.complex_data("train", self.cohort_names[lack])
            test_data, test_label = self.complex_data("test", self.cohort_names[lack]) 
            # print ("start training, used feature:", len(train_data[0]))

            clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
             min_samples_leaf = LEAF_S, random_state = np.random.seed(2021), max_features=NUM_FEATURE) 
            clf.fit(train_data, train_label)  
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            auc_total += roc_auc * len(test_label)
            print ("AUC", self.cohort_names[lack], roc_auc) 
            auc_list.append(roc_auc)
            auc_data_frame.append([round(roc_auc, 3), self.cohort_names[lack], group])
        print ("weighted mean:", auc_total/len(self.all_data))
        return auc_list, auc_total/len(self.all_data)
 
    def validation(self, select_feature_num):
        select_edges,marker_genus = self.combine_markers(select_feature_num)
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)

        print ("prepare data")
        train_data, train_label =  self.complex_data("train", "use all samples for train")

        # read validation set generated by additional_validation.py
        with open('validation_data.pkl', 'rb') as f:
            validation_data = pickle.load(f) 
        test_data, test_label = [], [] 
        for one_sample in validation_data:
            disease = one_sample[1]
            genus_abundance = one_sample[2]
            select_feature_array = one_sample[3]
            sample_array_HGT = np.dot(self.HGT_network, select_feature_array)
            sample_array_genus = np.dot(self.species_network, genus_abundance)
            sample_array = list(sample_array_HGT) + list(sample_array_genus) 
            if group == "Thomas-Abun":
                # sample_array = genus_abundance
                sample_array = extract_previous_16_markers(marker_genus, genus_abundance)
            test_data.append(sample_array)
            if disease == "control" :
                test_label.append(1)
            if disease == "CRC" :
                test_label.append(0)            
        test_data = np.array(test_data)
        test_label = np.array(test_label)
        print ("validation sample num:", len(test_data))
        print ("validation control num", sum(test_label))
        print ("validation CRC num", len(test_label) - sum(test_label))
      

        clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
            min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
        clf.fit(train_data, train_label)  
        roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])

        print ("AUC", "for validation", roc_auc) 
        return [], 0

    def get_HGT_network(self, select_edges):
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
        # print (event_array)          
        # laplacian = unnormalized_laplacian(event_array)
        # a, b = np.linalg.eig(laplacian)
        # self.HGT_network = b.T.real

        I = np.diag([1]*len(event_array))
        U = np.ones((len(event_array), len(event_array)))
        # self.HGT_network = event_array*hgt_alpha + I * hgt_beta + U * hgt_gamma # matrix not normalized
        self.HGT_network = normalize(event_array, axis=1, norm='l1')*hgt_alpha + I * hgt_beta + U * hgt_gamma

    def get_genus_network(self, marker_genus):
        self.species_network = np.zeros((len(marker_genus), len(marker_genus)))
        genus_index = {}
        index = 0
        id_conver = {} # save the index:marker pairs
        events = list(marker_genus.keys())
        for index in range(len(events)):
        # for genus in marker_genus:
            genus = events[index]
            id_conver[index] = genus
            if genus[0] == "s":
                marker = "g__" + genus.split("_")[2]
            else:
                marker = genus
            if marker not in genus_index:
                genus_index[marker] = [index]
            else:
                genus_index[marker] += [index]
            # index += 1
        
        nodes_list = list(marker_genus.keys())
        edges_list = []
        

        for sample in self.all_data:
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                    if array[0] in genus_index and array[1] in genus_index:
                        
                        a_list = genus_index[array[0]]
                        b_list = genus_index[array[1]]
                        for a in a_list:
                            for b in b_list:
                                self.species_network[a][b] += 1
                                self.species_network[b][a] += 1

        for i in range(len(self.species_network)):
            for j in range(len(self.species_network)):
                if self.species_network[i][j] > 200: # the HGT should be in more than 200 samples
                    self.species_network[i][j] = 1
                    if id_conver[i] != id_conver[j]:
                        edges_list.append([id_conver[i], id_conver[j]])
                else:
                    self.species_network[i][j] = 0
        with open('selected_abun_genus_graph.pkl', 'wb') as f:
            pickle.dump([nodes_list, edges_list], f) 
        print ("abundance marker num:", len(marker_genus))
        pearson_network, pearson_edges_list = self.get_genus_matrix_pearson(marker_genus, id_conver)
        with open('selected_abun_genus_graph_pearson.pkl', 'wb') as f:
            pickle.dump([nodes_list, pearson_edges_list], f) 
        I = np.diag([1]*len(marker_genus))
        U = np.ones((len(marker_genus), len(marker_genus)))
        # self.species_network = normalize(self.species_network, axis=1, norm='l1')
        # self.species_network = self.species_network * delta + pearson_network * zeta + I * eta + U * mu
        # vals, vecs = np.linalg.eig(pearson_network)
        # Lambda = np.diag(vals*0.5)
        # new_pearson_network = np.dot(np.dot(vecs, Lambda), vecs.T)
        new_pearson_network = normalize(pearson_network, axis=1, norm='l1')
        self.species_network = normalize(self.species_network, axis=1, norm='l1') * delta +  new_pearson_network* zeta + I * eta + U * mu

    def get_genus_matrix_pearson(self, marker_genus, id_conver):
        sample_matrix = np.zeros((len(marker_genus), len(self.all_data)))
        j = 0
        for sample in self.all_data:
            for i in range(len(sample.genus_abundance)):
                sample_matrix[i][j] = sample.genus_abundance[i]
            j += 1
        event_array = np.zeros((len(marker_genus), len(marker_genus)))
        events = list(marker_genus.keys())
        pearson_edges_list = []
        for i in range(len(events)):
            for j in range(len(events)):
                if i == j:
                    continue
                x = sample_matrix[i]
                y = sample_matrix[j]
                pearson = scipy.stats.pearsonr(x, y)

                if str(pearson[0]) == "nan":
                    # print (i, j, "no pearson")
                    continue
                # if abs(pearson[0]) > 0.3:
                #     event_array[i][j] = 1
                #     pearson_edges_list.append([id_conver[i], id_conver[j]])
                if pearson[0] > corr_para:
                    event_array[i][j] = pearson[0]
                    pearson_edges_list.append([id_conver[i], id_conver[j]])
                else:
                    event_array[i][j] = 0
                # event_array[i][j] = pearson[0]
                # pearson_edges_list.append([id_conver[i], id_conver[j]])
        return event_array, pearson_edges_list

    def complex_data(self, data_usage, select_cohort):
        data = []
        label = []
        for sample in self.all_data:
            if (sample.cohort == select_cohort and data_usage == "train"):
                continue
            if (sample.cohort != select_cohort and data_usage == "test"):
                continue
            sample_array_HGT = np.dot(self.HGT_network, sample.select_feature_array)
            sample_array_genus = np.dot(self.species_network, sample.genus_abundance)
            # proper_list,s_edge_num = sample.get_network_properties(level, 50)
            # network_property = sample.get_network_features(network_level)
            # network_property = sample.prepare_network_vector()
            sample_array = list(sample_array_HGT) + list(sample_array_genus) #+ proper_list #+ [proper_list[3]]#network_property
            # sample_array = list(sample.select_feature_array) + list(np.array(sample.genus_abundance)*1)
            # sample_array = list(sample.select_feature_array) + list(sample.marker_abundance)
            # sample_array = list(sample_array_genus)
            # sample_array = list(sample_array_HGT) + list(sample.marker_abundance)
            # sample_array = list(sample.select_feature_array) + list(sample.marker_abundance)
            # sample_array = sample_array_HGT
            # sample_array = sample.marker_abundance
            # sample_array = sample.select_feature_array
            if group == "Thomas-Abun":
                sample_array = sample.marker_abundance
            elif group == "HGT":
                sample_array = list(sample.select_feature_array)
            elif group == "Abun":
                sample_array = list(sample.genus_abundance)
            data.append(sample_array)
            if sample.disease == "control" :
                label.append(1)
            if sample.disease == "CRC" :
                label.append(0)            
            # label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

    def prepare_deep_learning(self, select_feature_num):
        select_edges,marker_genus = self.combine_markers(select_feature_num)
        print (select_edges)
        train_cohort = []
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
            train_cohort.append(sample.cohort)
        train_data, train_label =  self.complex_data("train", "use all to train")

        prediction_data = [list(train_data), list(train_label), train_cohort]
        with open("prediction_data.pkl", 'wb') as f:
            pickle.dump(prediction_data, f)

    def LODO_DNN(self, select_feature_num):
        select_edges,marker_genus = self.combine_markers(select_feature_num)
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
        auc_list = []   
        auc_total = 0    
        for lack in range(len(self.cohort_names)):
            print ("prepare data")
            train_data, train_label =  self.complex_data("train", self.cohort_names[lack])
            test_data, test_label = self.complex_data("test", self.cohort_names[lack]) 
            print ("start training, used feature:", len(train_data[0]))

            cla = Classifier()

            roc_auc = cla.run_model(train_data, train_label, test_data, test_label)
            auc_total += roc_auc * len(test_label)
            print ("AUC", self.cohort_names[lack], roc_auc) 
            auc_list.append(roc_auc)
            auc_data_frame.append([round(roc_auc, 3), self.cohort_names[lack], group])
        print ("weighted mean:", auc_total/len(self.all_data))
        return auc_list, auc_total/len(self.all_data)

    def read_T2D_data(self):
        T2D_data = []
        disease_num = {}

        acc_dir = "/mnt/d/breakpoints/HGT/CRC/T2D"
        files = os.listdir(acc_dir)
        for acc_file in files:
            if not re.search("acc.csv", acc_file):
                continue
            ID = acc_file.split(".")[0]
            # print (ID)
            acc_file = os.path.join(acc_dir, acc_file)
            sample = Sample(acc_file, ID)
            if sample.tag == "no pheno" or sample.tag == "no bkp":
                print (ID, sample.tag)
                continue         
            # if (sample.disease != "CRC") and (sample.disease != "control"):
            #     continue
            # print (sample.disease)
            if sample.disease not in disease_num:
                disease_num[sample.disease] = 0
            disease_num[sample.disease] += 1
            T2D_data.append(sample)
        print (len(T2D_data))
        print (disease_num)
        return T2D_data

    def t2d_validation(self, select_feature_num):
        select_edges,marker_genus = self.combine_markers(select_feature_num)
        print (select_edges)
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
        T2D_data = self.read_T2D_data()
        for sample in T2D_data:
            sample.given_nodes_make_matrix(select_edges)
        self.get_abundance_marker_value(marker_genus, T2D_data)

        print ("prepare data")
        train_data, train_label =  self.complex_data("train", "use all to train")
        # test_data, test_label =  self.complex_data("test", "ZellerG_2014")
        # test_data = list(test_data)
        # test_label = list(test_label)
        test_data, test_label = [], [] 
        for sample in T2D_data:
            sample_array_HGT = np.dot(self.HGT_network, sample.select_feature_array)
            sample_array_genus = np.dot(self.species_network, sample.genus_abundance)
            sample_array = list(sample_array_HGT) + list(sample_array_genus) #+ proper_list #+ [proper_list[3]]#network_property
            if group == "Thomas-Abun":
                sample_array = sample.marker_abundance
            test_data.append(sample_array)
            # if sample.disease == "control" :
            #     test_label.append(1)
            if sample.disease == "CRC" :
                test_label.append(0)  
            else:
                 test_label.append(1)         
        test_data = np.array(test_data)
        test_label = np.array(test_label)
        print ("validation sample num:", len(test_data))
        print ("validation control num", sum(test_label))
        print ("validation CRC num", len(test_label) - sum(test_label))
      

        clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
            min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
        clf.fit(train_data, train_label)  
        # roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
        # print ("AUC", "for validation", roc_auc) 
        accuracy = sklearn.metrics.accuracy_score(test_label, clf.predict(test_data))
        print ("accuracy", accuracy, "FPR", 1-accuracy )
        
        return [], 0

    def __del__(self):
        del self.all_data
        print ("object, deleted")


class Fast_RF():

    def __init__(self):
        self.all_data = []
        self.disease_sample_num_cohort = {}   
        self.HGT_network = ''
        self.species_network = ''
        self.pearson_network = ''
        self.revised_HGT_network = ''
        self.revised_species_network = ''
        self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b",\
            "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        self.sample_dict = {}


    def combine_markers(self, select_feature_num):
        select_HGTs, abun_related_HGTs = self.select_top_HGT(select_feature_num)
        # select_edges = {**select_HGTs, **abun_related_HGTs} # combine abundance-related HGT events
        select_edges = select_HGTs
        i = 0
        for edge in select_edges:
            select_edges[edge] = i
            i += 1
        self.get_HGT_network(select_edges)

        for edge in select_HGTs:  # add the species in HGT to abundance markers
            array = edge.split("&")
            if array[0] not in marker_genus:
                marker_genus[array[0]] = 1
            if array[1] not in marker_genus:
                marker_genus[array[1]] = 1
        i = 0
        for genus in marker_genus:
            marker_genus[genus] = i
            i += 1
        self.get_abundance_marker_value(marker_genus, self.all_data)
        self.get_genus_network(marker_genus)
        print ("HGT num", len(select_edges), "abun num", len(marker_genus))
        return select_edges, marker_genus
    
    def basic_run(self, select_feature_num):
        with open("sample_data", "rb") as fp:
            self.all_data = pickle.load(fp) 
        print ("data loaded")
        select_edges,marker_genus = self.combine_markers(select_feature_num)  
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
        self.save_sample_info()
        del self.all_data
    
    def repeat_run(self):
        self.network_revision()
        auc_list, weighted_mean = self.LODO()
        return auc_list, weighted_mean
    
    def network_revision(self):
        I = np.diag([1]*len(marker_genus))
        U = np.ones((len(marker_genus), len(marker_genus)))
        # self.species_network = normalize(self.species_network, axis=1, norm='l1')
        # self.revised_species_network = self.species_network * delta + self.pearson_network * zeta + I * eta + U * mu # origin
        # self.revised_species_network = normalize_digraph_method1(self.pearson_network)* zeta + I * eta + U * mu 
        # self.revised_species_network = normalizeAdjacency_method2(self.pearson_network)* zeta + I * eta + U * mu 
        # vals, vecs = np.linalg.eig(self.pearson_network)
        # Lambda = np.diag(vals**2)
        # new_pearson_network = np.dot(np.dot(vecs, Lambda), vecs.T)
        new_pearson_network = normalize(self.pearson_network, axis=1, norm='l1')
        # new_pearson_network = self.pearson_network
        # new_pearson_network = egen_denoise(self.pearson_network, test_para)

        self.revised_species_network = normalize(self.species_network, axis=1, norm='l1') * delta + new_pearson_network * zeta + I * eta + U * mu 
        # self.revised_species_network = self.species_network/np.sum(self.species_network) * delta + self.pearson_network/np.sum(self.pearson_network) * zeta + I/np.sum(I) * eta + U/np.sum(U) * mu
        # self.revised_species_network = normalize(self.species_network, axis=1, norm='l1')*hgt_alpha * delta + normalize(self.pearson_network, axis=1, norm='l1') * zeta + normalize(I, axis=1, norm='l1') * eta + normalize(U, axis=1, norm='l1') * mu
        # new_pearson_network, x = KRnorm_sym(self.pearson_network)
        # self.revised_species_network = new_pearson_network * zeta + I * eta + U * mu 
        # print (normalize(self.species_network, axis=1, norm='l1')* zeta)

        I = np.diag([1]*len(self.HGT_network))
        U = np.ones((len(self.HGT_network), len(self.HGT_network)))
        # self.revised_HGT_network = self.HGT_network*hgt_alpha + I * hgt_beta + U * hgt_gamma # origin
        self.revised_HGT_network = normalize(self.HGT_network, axis=1, norm='l1')*hgt_alpha + I * hgt_beta + U * hgt_gamma
        # self.revised_HGT_network = self.HGT_network/np.sum(self.HGT_network)*hgt_alpha + I/np.sum(I) * hgt_beta + U/np.sum(U) * hgt_gamma
        # self.revised_HGT_network = normalize(self.HGT_network, axis=1, norm='l1')*hgt_alpha + normalize(I, axis=1, norm='l1') * hgt_beta + normalize(U, axis=1, norm='l1') * hgt_gamma
        # new_hgt_network, x = KRnorm_sym(self.HGT_network)
        # self.revised_HGT_network = new_hgt_network*hgt_alpha + I * hgt_beta + U * hgt_gamma

    def get_abundance_marker_value(self, marker_genus, sample_list):

        sample_abd = get_genus_abd(marker_genus)
        # print (sample_abd)
        for sample in sample_list:
            if sample.ID not in phenotype.ID_name:
                sample_name = sample.ID
            else:
                sample_name = phenotype.ID_name[sample.ID]
            if sample_name not in sample_abd:
                sample_id = phenotype.name_sample_id[sample_name]
                genus_abundance = sample_abd[sample_id]
            else:
                genus_abundance = sample_abd[sample_name]
            sample.genus_abundance = genus_abundance

    def same_HGT_number(self, sample_dict, choose_num, min_split_num): #only choose top x HGTs in each sample
        new_same_dict = {}
        sort_sample_dict = sorted(sample_dict.items(), key=lambda item: item[1], reverse = True)
        if len(sort_sample_dict) < choose_num:
            choose_num = len(sort_sample_dict)
        for z in range(choose_num):
            if sort_sample_dict[z][1] < min_split_num:
                continue
            # new_same_dict[sort_sample_dict[z][0]] = sort_sample_dict[z][1]
            new_same_dict[sort_sample_dict[z][0]] = 1
        # print (sort_sample_dict[choose_num-1])
        return new_same_dict

    def select_top_HGT(self, select_feature_num):
        specific_HGT = {} 
        abun_related_HGTs = {}
        record_all_HGTs = {}
        crc_num = 0
        control_num = 0
        edge_distribution = []
        reads_num_list = []
        for sample in self.all_data:
            sample_dict = {}
            reads_num_list.append(sample.reads_num)
            min_split_num, p = get_split_reads_cutoff(sample.reads_num)
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                support_ratio = bkp.cross_split_reads/sample.reads_num
                # if support_ratio < min_cross:
                #     continue
                if support_ratio == 0:
                    continue
                # if bkp.cross_split_reads < min_split_num:
                #     continue
                # if bkp.pair_end/sample.reads_num < 2e-08:
                #     continue
                if edge not in record_all_HGTs:
                    record_all_HGTs[edge] = 1
                    specific_HGT[edge] = [[], []]
                if edge not in sample_dict:
                    sample_dict[edge] = 0
                # if support_ratio > sample_dict[edge]:
                #     sample_dict[edge] = support_ratio
                # sample_dict[edge] += support_ratio
                # sample_dict[edge] = 1
                sample_dict[edge] += bkp.cross_split_reads
            edge_distribution.append(len(sample_dict))
            # in each sample, choose same number of edge.       
            sample_dict = self.same_HGT_number(sample_dict, max_hgt_num, min_split_num)

            for edge in sample_dict:
                if sample.disease == "CRC" :
                    specific_HGT[edge][0].append(sample_dict[edge])
                else:
                    specific_HGT[edge][1].append(sample_dict[edge])
            if sample.disease == "CRC" :
                crc_num += 1
            if sample.disease == "control" :
                control_num += 1
        # print ("sample num :%s, Edge count: mean is %s, median is %s"%(len(edge_distribution),\
        #     np.mean(edge_distribution), np.median(edge_distribution)), np.std(edge_distribution))
        # print ("read num, median %s, mean %s"%(np.median(reads_num_list),np.mean(reads_num_list)), min(reads_num_list), max(reads_num_list))
        select_edges = {}
        # print (crc_num, control_num, crc_num+control_num, len(specific_HGT)/(crc_num+control_num))
        genus_level_markers = {}
        for marker in marker_genus:
            if marker[0] == "s":
                marker = "g__" + marker.split("_")[2]
            genus_level_markers[marker] = 1
        print ("abundance marker-related genus:", len(genus_level_markers))
        for tag in specific_HGT:
            if len(specific_HGT[tag][0]) + len(specific_HGT[tag][1]) < 25:
                continue
            array = tag.split("&")
            species_1 = array[0]
            species_2 = array[1]
            crc_array = specific_HGT[tag][0] + [0] * (crc_num - len(specific_HGT[tag][0]))
            control_array = specific_HGT[tag][1] + [0] * (control_num - len(specific_HGT[tag][1]))
            U1, p = mannwhitneyu(crc_array, control_array)
            if p < 0.01 and species_1 in genus_level_markers and species_2 in genus_level_markers:
                abun_related_HGTs[tag] = 1
            select_edges[tag] = p 
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

        final_select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            # if i < 3:
            #     print (sort_select_edges[i])
            final_select_edges[edge] = i
            if edge in abun_related_HGTs:
                del abun_related_HGTs[edge]
        print ("abun_related_HGTs", len(abun_related_HGTs))

        tree_edges = {}
        for i in range(len(sort_select_edges)): # store the all differential genome pairs
            # if sort_select_edges[i][1] > 0.05:
            #     break
            tag = sort_select_edges[i][0]
            crc_ratio = len(specific_HGT[tag][0])/crc_num
            control_ratio = len(specific_HGT[tag][1])/control_num
            tree_edges[sort_select_edges[i][0]] = [crc_ratio, control_ratio, sort_select_edges[i][1]]
        # with open('selected_diff_edges.pkl', 'wb') as f:
        #     pickle.dump(tree_edges, f) 
          
        return final_select_edges, abun_related_HGTs

    def all_samples(self):
        all_acc_file = "acc.list"
        # result_dir = "new_result_more"
        result_dir = "new_result"
        os.system(f"ls {result_dir}/*acc.csv>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            ID = acc_file.split("/")[1].split(".")[0]
            new_file = acc_file.split("/")[1]
            acc_file = f"{result_dir}/{new_file}"
            
            sample = Sample(acc_file, ID)

            if sample.tag == "no pheno" or sample.tag == "no bkp":
                print (ID, sample.tag)
                continue  
            # if sample.cohort == "FengQ_2015":
            #     print (sample.cohort, ID, sample.disease)          
            if (sample.disease != "CRC") and (sample.disease != "control"):
                continue
            if str(sample.disease) == "nan":
                continue
            if sample.cohort == "HanniganGD_2017":
                continue
            self.all_data.append(sample)
            if sample.cohort not in self.disease_sample_num_cohort:
                self.disease_sample_num_cohort[sample.cohort] = {"CRC":0, "control":0, "adenoma":0}
            self.disease_sample_num_cohort[sample.cohort][sample.disease] += 1

        for key in self.disease_sample_num_cohort:
            print ("Cohort", key, self.disease_sample_num_cohort[key])     

    def LODO(self):
        auc_list = []   
        auc_total = 0    
        for lack in range(len(self.cohort_names)):
            # print ("prepare data")
            train_data, train_label =  self.complex_data("train", self.cohort_names[lack])
            test_data, test_label = self.complex_data("test", self.cohort_names[lack]) 
            # print ("start training, used feature:", len(train_data[0]))

            clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
             min_samples_leaf = LEAF_S, random_state = np.random.seed(2021), max_features=NUM_FEATURE) 
            clf.fit(train_data, train_label)  
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            auc_total += roc_auc * len(test_label)
            # print ("AUC", self.cohort_names[lack], roc_auc) 
            auc_list.append(roc_auc)
            auc_data_frame.append([round(roc_auc, 3), self.cohort_names[lack], group])
        print ("unweight mean", np.mean(auc_list))
        return auc_list, auc_total/len(self.sample_dict)
 
    def get_HGT_network(self, select_edges):
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
        self.HGT_network = event_array

    def get_genus_network(self, marker_genus):
        self.species_network = np.zeros((len(marker_genus), len(marker_genus)))
        genus_index = {}
        index = 0
        id_conver = {} # save the index:marker pairs
        events = list(marker_genus.keys())
        for index in range(len(events)):
        # for genus in marker_genus:
            genus = events[index]
            id_conver[index] = genus
            if genus[0] == "s":
                marker = "g__" + genus.split("_")[2]
            else:
                marker = genus
            if marker not in genus_index:
                genus_index[marker] = [index]
            else:
                genus_index[marker] += [index]
            # index += 1
        
        nodes_list = list(marker_genus.keys())
        edges_list = []
        

        for sample in self.all_data:
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                    if array[0] in genus_index and array[1] in genus_index:
                        
                        a_list = genus_index[array[0]]
                        b_list = genus_index[array[1]]
                        for a in a_list:
                            for b in b_list:
                                self.species_network[a][b] += 1
                                self.species_network[b][a] += 1

        for i in range(len(self.species_network)):
            for j in range(len(self.species_network)):
                if self.species_network[i][j] > 200: # the HGT should be in more than 200 samples
                    self.species_network[i][j] = 1
                    if id_conver[i] != id_conver[j]:
                        edges_list.append([id_conver[i], id_conver[j]])
                else:
                    self.species_network[i][j] = 0
        print ("abundance marker num:", len(marker_genus))
        self.pearson_network = self.get_genus_matrix_pearson(marker_genus, id_conver)

    def get_genus_matrix_pearson(self, marker_genus, id_conver):
        sample_matrix = np.zeros((len(marker_genus), len(self.all_data)))
        j = 0
        for sample in self.all_data:
            for i in range(len(sample.genus_abundance)):
                sample_matrix[i][j] = sample.genus_abundance[i]
            j += 1
        event_array = np.zeros((len(marker_genus), len(marker_genus)))
        events = list(marker_genus.keys())
        for i in range(len(events)):
            for j in range(len(events)):
                if i == j:
                    continue
                x = sample_matrix[i]
                y = sample_matrix[j]
                pearson = scipy.stats.pearsonr(x, y)

                if str(pearson[0]) == "nan":
                    # print (i, j, "no pearson")
                    continue
                # if abs(pearson[0]) > 0.3:
                # if pearson[0] > 0.3:
                #     event_array[i][j] = 1
                #     pearson_edges_list.append([id_conver[i], id_conver[j]])
                if pearson[0] > corr_para:
                    event_array[i][j] = pearson[0]
                else:
                    event_array[i][j] = 0
                # event_array[i][j] = pearson[0]
                # pearson_edges_list.append([id_conver[i], id_conver[j]])
                # else:
                #     event_array[i][j] = 0
        print ("###")
        return event_array

    def complex_data(self, data_usage, select_cohort):
        data = []
        label = []
        for sample in self.sample_dict:
            info = self.sample_dict[sample]
            if (info[0] == select_cohort and data_usage == "train"):
                continue
            if (info[0] != select_cohort and data_usage == "test"):
                continue
            # sample_array_HGT = np.dot(self.revised_HGT_network, self.norm_array(info[1]))
            # sample_array_genus = np.dot(self.revised_species_network, self.norm_array(info[2]))
            
            sample_array_HGT = np.dot(self.revised_HGT_network, info[1])
            sample_array_genus = np.dot(self.revised_species_network, np.array(info[2])) 
            sample_array = list(sample_array_HGT) + list(sample_array_genus) #+ proper_list #+ [proper_list[3]]#network_property
            # print (list(info[1]), list(sample_array_HGT), info[2], list(sample_array_genus))
            # print (info[2], list(sample_array_genus), list(np.array(info[2])-sample_array_genus))
            # sample_array = list(info[1]) + list(np.array(info[2]))
            # sample_array = list(sample.genus_abundance)
            # sample_array = list(info[1]) +  list(info[2])
            # sample_array = list(info[1])
            # sample_array = list(self.norm_array(info[1])) + list(self.norm_array(info[2]))
            # if group == "Thomas-Abun":
            #     sample_array = sample.marker_abundance
            # elif group == "HGT":
            #     sample_array = list(sample.select_feature_array)
            # elif group == "Abun":
            #     sample_array = list(sample.genus_abundance)
            data.append(sample_array)
            if info[3] == "control" :
                label.append(1)
            if info[3] == "CRC" :
                label.append(0)            
        data = np.array(data)
        label = np.array(label)
        return data, label
    
    def norm_array(self, x):
        x = np.array(x)
        norm  =  np.linalg.norm(x)
        if norm == 0: 
            return x
        return x/norm
    
    def save_sample_info(self):
        for sample in self.all_data:
            self.sample_dict[sample.ID] = [sample.cohort, sample.select_feature_array, sample.genus_abundance, sample.disease]


def extract_previous_16_markers(marker_genus, genus_abundance):
    marker_species_dict = get_abd() # previous 16 abundance markers
    previous_abundance = []
    for genus in marker_genus:
        index = marker_genus[genus]
        abundance = genus_abundance[index]
        if genus in marker_species_dict:
            previous_abundance.append(abundance)
    return previous_abundance

def get_split_reads_cutoff(g): # g is the number of reads in the specific sample
    # prior_a = 25 # prior probability
    prior_b = 50000000#63333330  # prior probability
    given_n = 42648185 # mean number of reads among samples
    # # n = 182534663 # max number of reads 
    # given_r = 2 # the cutoff with n reads  4

    for m in range(100):
        alpha = prior_a + m
        beta= prior_b + g - m
        p = 1
        for k in range(given_r + 1):
            choose = comb(given_n, k)
            beta1 = sc.beta(alpha + k, beta + given_n - k)
            beta2 = sc.beta(alpha, beta)
            # print (choose, beta1)
            p -= choose * beta1 / beta2
        if p > 0.9: #0.9
            break
        # print (m, p)
        # print (r, p, sc.beta(alpha+k, beta+n-k), sc.beta(alpha, beta))
    return m, p

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
    i = 0
    for species in marker_species:
        if species[-1] == ".":
            genus = "g__" + species.split()[0]
        else:
            genus = "s__" + "_".join(species.split())
        if genus not in marker_genus:
            marker_genus[genus] = i
            i += 1
    return marker_genus

def find_all_class():
    class_dict = {}
    i = 0
    for lineage in taxonomy.taxonomy_dict.values():
        my_class = lineage.split(";")[2]
        if my_class not in class_dict:
            class_dict[my_class] = i
            i += 1
    # print (class_dict.keys())
    # print (len(class_dict))
    return class_dict

def normalize_digraph_method1(A):
    Dl = np.sum(A, 0)  
    num_node = A.shape[0]
    Dn = np.zeros((num_node, num_node))
    for i in range(num_node):
        if Dl[i] > 0:
            Dn[i, i] = Dl[i]**(-1) 
    AD = np.dot(A, Dn)
    return AD


if __name__ == "__main__":
    level = 5
    TREE_NUM = 1000
    LEAF_S = 5
    NUM_FEATURE = "sqrt"
    # delta, zeta, eta, mu  =  1e-8, 2e-5, 0.4, 1e-10    # Genus network
    # hgt_alpha, hgt_beta, hgt_gamma = 1e-11, 1, 1e-7  # HGT network   before 2023/1/10

    # delta, zeta, eta, mu  =  9e-09, 7e-05, 0.4, 1e-10    # Genus network
    # hgt_alpha, hgt_beta, hgt_gamma = 2.4e-07, 1, 1e-7  # HGT network   
    # feature_num = 21   #64 for 0.870  24 for 0.864  10
    feature_num, TREE_NUM, LEAF_S, NUM_FEATURE = 20,  1000, 5, 2
    delta, hgt_alpha, hgt_beta, hgt_gamma, eta, zeta, mu = 0, 0, 1, 1e-8, 1, 0.048, 0
    corr_para = 0.4

    group = ""
    auc_data_frame = []
    network_level = 2

    max_hgt_num = 2150  #2150
    min_cross = 1e-11
    pseudo = 0
    svd_num = 100000
    my_methods = ["density", "laplacian", "eigenvalue", "degree"]
    select_m = my_methods[2]
    i=2

    
    result_file = open("random_forest.log", 'w')

    prior_a = 20 # prior probability
    given_r = 2 # the cutoff with n reads  4
    
    marker_genus = select_genus()
    sample_abd =  get_genus_abd(marker_genus)   
    phenotype = Phenotype() 
    taxonomy = Taxonomy()
    class_dict = find_all_class()

    ### test the classifier in the independent CRC cohort and T2D cohort
    group = "Hybrid"
    rf = RF()
    # rf.t2d_validation(feature_num) # independent T2D cohort
    # auc_list, weighted_mean = rf.LODO(feature_num)
    # must run python additional_validation.py 
    auc_list, weighted_mean = rf.validation(feature_num) # independent CRC cohort


    ### compare the integration of HGT and abundance biomarkers and the previously reported 16 biomarkers
    # for group in ["Hybrid", "Thomas-Abun"]:  # for main plot
    #     marker_genus = select_genus()
    #     sample_abd =  get_genus_abd(marker_genus)
    #     phenotype = Phenotype()
    #     taxonomy = Taxonomy()
    #     rf = RF()
    
    #     auc_list, weighted_mean = rf.LODO(feature_num)
    #     print (group, auc_list, weighted_mean)
    #     del rf
    # df = pd.DataFrame(auc_data_frame, columns = ["AUC", "Cohort", "Group"])
    # df.to_csv('./for_AUC_plot.csv', sep='\t')
    # os.system("Rscript plot_auc.R")


    ### compare the integration of HGT and abundance biomarkers, the previously reported 16 biomarkers, only HGT, and only our abundance. 
    # for group in ["Hybrid", "Thomas-Abun", "HGT", "Abun"]:  # for supplementary plot
    #     marker_genus = select_genus()
    #     sample_abd =  get_genus_abd(marker_genus)
    #     phenotype = Phenotype()
    #     taxonomy = Taxonomy()
    #     rf = RF()
    
    #     auc_list, weighted_mean = rf.LODO(feature_num)
    #     print (group, auc_list, weighted_mean)
    #     del rf
    # df = pd.DataFrame(auc_data_frame, columns = ["AUC", "Cohort", "Group"])
    # df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//for_AUC_plot.csv', sep='\t')



    ### for hyper-parameters assessment
    # marker_genus = select_genus()   
    # sample_abd =  get_genus_abd(marker_genus) 
    # # print ("basic run is finished")
    # # # hgt_alpha, hgt_beta, hgt_gamma = 1e-08, 1, 1e-09  # HGT network
    # # # delta, zeta, eta, mu =  1e-09, 0.0001, 0.9, 1e-08   # Genus network
    # result_file = open("random_forest.log", 'w+')
    # small_values = [ 0, 0.1, 0.05, 0.01, 0.005,  1e-3, 0.0005, 1e-4, 1e-5]
    
    # # for feature_num in range(16, 25):
    # corr_para = 0.4
    # zeta = 0.048
    # fast = Fast_RF()
    # fast.basic_run(feature_num)
    # auc_list, weighted_mean = fast.repeat_run()
    # weighted_mean = round(weighted_mean, 4)
    # print (TREE_NUM, LEAF_S, NUM_FEATURE, "hgt_alpha:%s,hgt_beta:%s,hgt_gamma:%s, delta:%s, zeta:%s, eta:%s, mu:%s"%(hgt_alpha, hgt_beta, hgt_gamma, delta, zeta, eta, mu), weighted_mean )
    # print ("hgt_alpha:%s,hgt_beta:%s,hgt_gamma:%s, delta:%s, zeta:%s, eta:%s, mu:%s"%(hgt_alpha, hgt_beta, hgt_gamma, delta, zeta, eta, mu), weighted_mean, file=result_file)
    # result_file.close()

    













    # def Lachnospiraceae(self):
    #     specific_HGT = {} 
    #     family_dict = {}
    #     crc_num = 0
    #     control_num = 0
    #     best_num = 0
    #     for sample in self.all_data:
    #         if sample.disease == "CRC" :
    #             crc_num += 1
    #         if sample.disease == "control" :
    #             control_num += 1
    #         sample_dict = {}
    #         # min_split_num, p = get_split_reads_cutoff(sample.reads_num)
    #         for bkp in sample.bkps:
    #             edge = get_tag(bkp, level)
    #             array = edge.split("&")
    #             from_family = bkp.from_ref_lineage.split(";")[4]
    #             to_family = bkp.from_ref_lineage.split(";")[4]
    #             family_dict[edge] = [from_family, to_family]
    #             if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
    #                 continue
    #             if edge in sample_dict:
    #                 continue
    #             else:
    #                 sample_dict[edge] = 1
    #             if edge not in specific_HGT:
    #                 specific_HGT[edge] = [0, 0]
    #             if sample.disease == "CRC" :
    #                 specific_HGT[edge][0] += 1
    #             if sample.disease == "control" :
    #                 specific_HGT[edge][1] += 1    
    #     select_edges = {}
    #     for tag in specific_HGT:
    #         if specific_HGT[tag][0] + specific_HGT[tag][1] < 25:
    #             continue
    #         array = tag.split("&")
    #         species_1 = array[0]
    #         species_2 = array[1]
    #         crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
    #         control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
    #         U1, p = mannwhitneyu(crc_array, control_array)
    #         select_edges[tag] = p 
    #     sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)


    #     # tree_edges = {}
    #     data = []
    #     for i in range(30):
    #         tag = sort_select_edges[i][0]
    #         # print (i, sort_select_edges[i])
    #         if family_dict[tag][0] == "f__Lachnospiraceae" or family_dict[tag][1] == "f__Lachnospiraceae":
    #             print (family_dict[tag], tag, specific_HGT[tag])
    #             data.append([tag, specific_HGT[tag][0]/crc_num, "CRC" ])
    #             data.append([tag, specific_HGT[tag][1]/control_num, "control"])
    #     df = pd.DataFrame(data, columns = ["HGT", "Proportion", "Group"])
    #     df.to_csv('for_Lachnospiraceae.csv', sep='\t')
