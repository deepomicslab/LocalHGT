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
import scipy
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
        self.read_sra_meta("/mnt/d/breakpoints/script/analysis/new_result/usa_canada.csv")
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/austria/austria.csv")     
        
    def read_sra_meta(self, sra_meta):
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
            if row[0] == "from_ref":
                continue
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
            new_tag = get_tag(bkp, level)
            if new_tag in select_edges:
                self.select_feature_array[select_edges[new_tag]] = 1 

    def get_HGT_matrix(self, level):
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
        HGT_matrix = nx.from_numpy_matrix(HGT_matrix)
        density = nx.density(HGT_matrix)
        # density = nx.transitivity(HGT_matrix)
        return density

    def get_common_matrix(self, select_edges, level):
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
        HGT_matrix = np.zeros((len(nodes), len(nodes)))
        nodes_index = {}
        for i in range(len(nodes)):
            nodes_index[nodes[i]] = i

        for bkp in self.bkps:
            new_tag = get_tag(bkp, level)
            if new_tag in select_edges:
                array = edge.split("&")
                node1 = nodes_index[array[0]]
                node2 = nodes_index[array[1]]
                HGT_matrix[node1][node2] = 1
                HGT_matrix[node2][node1] = 1
        HGT_matrix = nx.from_numpy_matrix(HGT_matrix)
        density = nx.density(HGT_matrix)
        return density

class RF():

    def __init__(self, common_ratio):
        self.all_data = []
        self.disease_sample_num_cohort = {}

        # self.all_samples()
        # with open("sample_data", "wb") as fp:
        #     pickle.dump(self.all_data, fp)

        with open("sample_data", "rb") as fp:
            self.all_data = pickle.load(fp)

        self.cohorts_names = {}        
        self.common_ratio = common_ratio #0.01
        self.event_array = ''
        self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b",\
        "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        # self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b","ZellerG_2014","YuJ_2015"]
        # self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018",\
        # "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        # self.cohort_names = ["FengQ_2015", "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        print ("RF init done")

    def all_samples(self):
        all_acc_file = "acc.list"
        os.system(f"ls new_result/*acc.csv>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            ID = acc_file.split("/")[1].split(".")[0]
            new_file = acc_file.split("/")[1]
            acc_file = f"new_result/{new_file}"
            
            sample = Sample(acc_file, ID)

            if sample.tag == "no pheno" or sample.tag == "no bkp":
                print (ID, sample.tag)
                continue  
            # if sample.cohort == "FengQ_2015":
            #     print (sample.cohort, ID, sample.disease)          
            if (sample.disease != pheno_group_1) and (sample.disease != pheno_group_2):
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

    def select_ranking(self, ignore_cohort):
        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0

        for sample in self.all_data:
            if ignore_cohort != "all":
                if sample.cohort == ignore_cohort:
                    continue
            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == pheno_group_2:
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == pheno_group_1:
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
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
                if edge not in cohort_sam_HGT[sample.cohort]:
                    cohort_sam_HGT[sample.cohort][edge] = [0, 0]
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 
        i = 0
        select_edges = {}
        for tag in specific_HGT:
            if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) <= self.common_ratio:
                continue
            if tag not in select_edges:
                select_edges[tag] = i
                i += 1
        print ( "from %s choose %s edges."%(len(specific_HGT), len(select_edges)))


        train_data, train_label =  [], []
        for sample in self.all_data:
            if ignore_cohort != "all":
                if sample.cohort == ignore_cohort:
                    continue
            sample.given_nodes_make_matrix(select_edges)
            sample_array = sample.select_feature_array
            train_data.append(sample_array)
            train_label.append(sample.disease)
        print ("start ranking")
        clf = RandomForestClassifier(n_estimators=1000, criterion="entropy",\
         min_samples_leaf=5,random_state = np.random.seed(2021), oob_score = True) 
        clf.fit(train_data, train_label)  
        np.save('importances.npy',clf.feature_importances_) 
        # result = permutation_importance(clf, train_data, train_label, n_repeats=5, random_state=np.random.seed(2021))        
        # np.save('importances.npy', result.importances_mean) 
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

        select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            select_edges[edge] = i
        return select_edges

    def select_feature(self, ignore_cohort, select_feature_num):
        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0

        for sample in self.all_data:
            if ignore_cohort != "all":
                if sample.cohort == ignore_cohort:
                    continue
            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == pheno_group_2:
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == pheno_group_1:
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
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
                if edge not in cohort_sam_HGT[sample.cohort]:
                    cohort_sam_HGT[sample.cohort][edge] = [0, 0]
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 
        i = 0
        select_edges = {}
        for tag in specific_HGT:
            if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) <= self.common_ratio:
                continue
            # """
            crc_ratio = float(specific_HGT[tag][0])/crc_num
            control_ratio = float(specific_HGT[tag][1])/control_num
            ratio_list = []
            for cohort in cohort_sam_HGT:
                if tag not in cohort_sam_HGT[cohort]:
                   ratio_list.append(0)
                else:
                    crc_ratio = float(cohort_sam_HGT[cohort][tag][0])/cohort_sam_num[cohort][0]
                    control_ratio = float(cohort_sam_HGT[cohort][tag][1])/cohort_sam_num[cohort][1]
                    max_ratio = max([crc_ratio, control_ratio])
                    # diff = abs(crc_ratio - control_ratio)
                    ratio_list.append(max_ratio)
            ratio_list = sorted(ratio_list)
            select_edges[tag] = abs(crc_ratio - control_ratio)
        print ("from %s choose %s edges."%(len(specific_HGT), len(select_edges)))

        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = True)
        print (sort_select_edges[:5])

        select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            select_edges[edge] = i
        return select_edges

    def select_feature_2(self, ignore_cohort, select_feature_num):
        with open('best_features.pkl', 'rb') as f:
            best_select_edges = pickle.load(f)
        # select_edges = {}
        # i = 0
        # for edge in list(best_select_edges.keys())[:70]:
        #     select_edges[edge] = i    
        #     i += 1


        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.all_data:
            if ignore_cohort != "all":
                if sample.cohort == ignore_cohort:
                    continue
            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == pheno_group_2:
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == pheno_group_1:
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
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
                if edge not in cohort_sam_HGT[sample.cohort]:
                    cohort_sam_HGT[sample.cohort][edge] = [0, 0]
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 
        i = 0
        select_edges = {}

        for tag in specific_HGT:
            if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) <= self.common_ratio:
                continue
            # """
            crc_ratio = float(specific_HGT[tag][0])/crc_num
            control_ratio = float(specific_HGT[tag][1])/control_num
            ratio_list = []
            p_list = []
            for cohort in cohort_sam_HGT:
                if tag not in cohort_sam_HGT[cohort]:
                   ratio_list.append(0)
                   p_list.append(1)
                else:
                    crc_ratio = float(cohort_sam_HGT[cohort][tag][0])/cohort_sam_num[cohort][0]
                    control_ratio = float(cohort_sam_HGT[cohort][tag][1])/cohort_sam_num[cohort][1]
                    # diff = abs(crc_ratio - control_ratio)
                    diff = max(crc_ratio, control_ratio)

                    crc_array = [1] * cohort_sam_HGT[cohort][tag][0] + [0] * (cohort_sam_num[cohort][0] - cohort_sam_HGT[cohort][tag][0])
                    control_array = [1] * cohort_sam_HGT[cohort][tag][1] + [0] * (cohort_sam_num[cohort][1] - cohort_sam_HGT[cohort][tag][1])
                    U1, p = mannwhitneyu(crc_array, control_array, method="auto")
                    p_list.append(p)
                    ratio_list.append(diff)
            if sorted(ratio_list)[4] < self.common_ratio:
                continue
            # ratio_list = sorted(ratio_list)
            hit_num = specific_HGT[tag][0] + specific_HGT[tag][1]
            hit_ratio = hit_num/(crc_num + control_num)
            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
            U1, p = mannwhitneyu(crc_array, control_array, method="auto")

            if sorted(ratio_list)[4] < self.common_ratio:
                continue
            if tag in best_select_edges:
                # print (tag, p, sorted(p_list))
                best_num += 1
            
            select_edges[tag] = i
            i += 1

        print ("from %s choose %s edges."%(len(specific_HGT), len(select_edges)), "filter best", best_num)

        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = True)
        print (sort_select_edges[:5])
        print (sort_select_edges[select_feature_num-5:select_feature_num])
        new_best = 0
        final_select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            if edge in best_select_edges:
                new_best += 1
            final_select_edges[edge] = i
        print ("final best", new_best)
        j = len(final_select_edges)
        for edge in best_select_edges:
            if edge not in final_select_edges:
                print ("lack", edge)
                # if select_edges[edge] < 0.01:
                # final_select_edges[edge] = j
                # j += 1
        
        return final_select_edges


        # random forest
        """
        train_data, train_label =  [], []
        for sample in self.all_data:
            if ignore_cohort != "all":
                if sample.cohort == ignore_cohort:
                    continue
            sample.given_nodes_make_matrix(select_edges)
            sample_array = sample.select_feature_array
            train_data.append(sample_array)
            train_label.append(sample.disease)
        print ("start ranking")


        clf = RandomForestClassifier(n_estimators=1000, criterion="gini", max_depth = 2, \
        min_samples_leaf=5,random_state = np.random.seed(2021), oob_score = True) 
        clf.fit(train_data, train_label) 
        np.save('importances.npy',clf.feature_importances_) 
        with open('features.pkl', 'wb') as f:
            pickle.dump(select_edges, f)
        print ("feature ranking score saved")
        select_edges = self.get_top_feature(select_feature_num)
        new_best = 0
        for edge in select_edges:
            if edge in best_select_edges:
                new_best += 1
        print ("final best", new_best)
        return select_edges
        """

    def select_feature_3(self, ignore_cohort, select_feature_num):
        with open('best_features.pkl', 'rb') as f:
            best_select_edges = pickle.load(f)
        for edge in best_select_edges:
            best_select_edges[edge] = []
        all_select_edges = {}

        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.all_data:
            if sample.disease == pheno_group_2:
                crc_num += 1
            if sample.disease == pheno_group_1:
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
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1
        
        
        select_edges = {}
        print (crc_num, control_num, crc_num+control_num, len(specific_HGT)/(crc_num+control_num))
        for tag in specific_HGT:
            if specific_HGT[tag][0] + specific_HGT[tag][1] < 25:
                continue
            # if max([float(specific_HGT[tag][0])/crc_num, float(specific_HGT[tag][1])/control_num]) < self.common_ratio:
            #     continue
            crc_ratio = float(specific_HGT[tag][0])/crc_num
            control_ratio = float(specific_HGT[tag][1])/control_num

            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])

            U1, p = mannwhitneyu(crc_array, control_array)
            # if p > 0.05:
            #     continue
            if tag in best_select_edges:
                best_num += 1

            select_edges[tag] = p #-1 * abs(crc_ratio - control_ratio) #p #

            # print (cohort, "from %s choose %s edges."%(len(specific_HGT), len(select_edges)))
        print ("&&")
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

        new_best = 0
        final_select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)

        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            print (sort_select_edges[i])
            if edge in best_select_edges:
                new_best += 1
            final_select_edges[edge] = i
        print (sort_select_edges[26])


        return final_select_edges

    def select_feature_4(self, ignore_cohort, select_feature_num,select_edges):
        train_data, train_label =  [], []
        for sample in self.all_data:
            if ignore_cohort != "all":
                if sample.cohort == ignore_cohort:
                    continue
            sample.given_nodes_make_matrix(select_edges)
            sample_array = sample.select_feature_array
            train_data.append(sample_array)
            train_label.append(sample.disease)
        print ("start ranking")
        clf = RandomForestClassifier(n_estimators=1000, criterion="entropy",\
         min_samples_leaf=5,random_state = np.random.seed(2021), oob_score = True) 
        clf.fit(train_data, train_label)  
        np.save('importances.npy',clf.feature_importances_) 
        # result = permutation_importance(clf, train_data, train_label, n_repeats=5, random_state=np.random.seed(2021))        
        # np.save('importances.npy', result.importances_mean) 
        with open('features.pkl', 'wb') as f:
            pickle.dump(select_edges, f)
        print ("feature ranking score saved")
        select_edges = self.get_top_feature(select_feature_num)
        return select_edges

    def select_feature_5(self, diff_cut):
        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0

        for sample in self.all_data:

            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == pheno_group_2:
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == pheno_group_1:
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
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
                if edge not in cohort_sam_HGT[sample.cohort]:
                    cohort_sam_HGT[sample.cohort][edge] = [0, 0]
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 
        with open('best_features.pkl', 'rb') as f:
            best_select_edges = pickle.load(f)

            filtered_HGT = {}
            HGT_index = 0
            for HGT in specific_HGT:
                existing_sample_num =  specific_HGT[HGT][0] + specific_HGT[HGT][1]

                # if existing_sample_num < 25 or abs(specific_HGT[HGT][0] - specific_HGT[HGT][1]) < 2: # 25 2
                #     continue
                # if max(specific_HGT[HGT]) < 22: # 22
                #     continue

                if existing_sample_num < 25:
                    continue

                cohort_num = []
                cohort_diff = []
                for cohort in cohort_sam_HGT:
                    if HGT in cohort_sam_HGT[cohort]:
                        cohort_num.append(sum(cohort_sam_HGT[cohort][HGT]))
                        cohort_diff.append(abs(cohort_sam_HGT[cohort][HGT][0] - cohort_sam_HGT[cohort][HGT][1]))
                    else:
                        cohort_num.append(0)
                        cohort_diff.append(0)


                crc_array = [1] * specific_HGT[HGT][0] + [0] * (crc_num - specific_HGT[HGT][0])
                control_array = [1] * specific_HGT[HGT][1] + [0] * (control_num - specific_HGT[HGT][1])

                U1, p = mannwhitneyu(crc_array, control_array, method="auto")
                if p > diff_cut:
                    continue

                # if sum(cohort_diff) < diff_cut: # 19  35  40
                #     continue
                # if abs(specific_HGT[HGT][0] - specific_HGT[HGT][1]) < diff_cut:
                #     continue
                array = HGT.split("&")
                species_1 = array[0]
                species_2 = array[1]
                filtered_HGT[HGT] = HGT_index 
                HGT_index += 1
            best_num = 0
            for HGT in best_select_edges:
                if HGT in filtered_HGT:
                    best_num += 1
            print ("filtered HGT num:", len(filtered_HGT), "best num:", best_num)
            with open("filtered_HGTs.pkl", "wb") as fp:
                pickle.dump(filtered_HGT, fp)

    def search_pattern(self):
        specific_HGT = {} 
        cohort_sam_num = {}
        cohort_sam_HGT = {}
        crc_num = 0
        control_num = 0

        with open('best_features.pkl', 'rb') as f:
            best_select_edges = pickle.load(f)

        for sample in self.all_data:
            # if ignore_cohort != "all":
            #     if sample.cohort == ignore_cohort:
            #         continue
            if sample.cohort not in cohort_sam_num:
                cohort_sam_num[sample.cohort] = [0, 0]
            if sample.cohort not in cohort_sam_HGT:
                cohort_sam_HGT[sample.cohort] = {}
            if sample.disease == pheno_group_2:
                crc_num += 1
                cohort_sam_num[sample.cohort][0] += 1
            if sample.disease == pheno_group_1:
                control_num += 1
                cohort_sam_num[sample.cohort][1] += 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                if edge not in best_select_edges:
                    continue
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
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                    cohort_sam_HGT[sample.cohort][edge][0] += 1 
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1
                    cohort_sam_HGT[sample.cohort][edge][1] += 1 
        i = 0
        # for tag in specific_HGT:
        for tag in best_select_edges:
            crc_ratio = float(specific_HGT[tag][0])/crc_num
            control_ratio = float(specific_HGT[tag][1])/control_num

            ratio_list = []
            for cohort in cohort_sam_HGT:
                if tag not in cohort_sam_HGT[cohort]:
                   ratio_list.append(0)
                else:
                    cohort_crc_ratio = float(cohort_sam_HGT[cohort][tag][0])/cohort_sam_num[cohort][0]
                    cohort_control_ratio = float(cohort_sam_HGT[cohort][tag][1])/cohort_sam_num[cohort][1]
                    # diff = cohort_crc_ratio - cohort_control_ratio
                    # diff = abs(diff)
                    # diff = round(diff, 2)
                    diff = max([cohort_crc_ratio, cohort_control_ratio])
                    ratio_list.append(diff)
            hit_ratio = (specific_HGT[tag][0] +specific_HGT[tag][1])/(crc_num + control_num)
            print (tag, max(ratio_list))

    def LODO(self, select_feature_num):
        # shuffle(self.all_data)
        auc_list = []
        # self.select_ranking("all")  
        # select_edges = self.get_top_feature(select_feature_num)

        select_edges = self.select_feature_3("all", select_feature_num)
        # with open('selected_26_edges.pkl', 'wb') as f:
        #     pickle.dump(select_edges, f)
        # with open('filtered_HGTs.pkl', 'rb') as f:
        #     select_edges = pickle.load(f)
        # select_edges = self.select_feature_4("all", select_feature_num, select_edges)

        # with open('best_features.pkl', 'wb') as f:
        #     pickle.dump(select_edges, f)
        sample_total  = 0
        auc_total = 0       
        for lack in range(len(self.cohort_names)):
            # self.select_ranking(self.cohort_names[lack])
            # select_edges = self.get_top_feature(select_feature_num)

            # select_edges = self.select_feature_3(self.cohort_names[lack], select_feature_num)

            # with open('filtered_HGTs.pkl', 'rb') as f:
            #     select_edges = pickle.load(f)
            # select_edges = self.select_feature_4(self.cohort_names[lack], select_feature_num, select_edges)
            print ("%s HGT events were selected."%(len(select_edges)))
            

            self.get_event_matrix(select_edges)
            for sample in self.all_data:
                sample.given_nodes_make_matrix(select_edges)

            print ("prepare data")
            train_data, train_label =  self.complex_data("train", self.cohort_names[lack])
            test_data, test_label = self.complex_data("test", self.cohort_names[lack]) 

            # train_data, train_label =  self.complex_data_density("train", self.cohort_names[lack])
            # test_data, test_label = self.complex_data_density("test", self.cohort_names[lack]) 
            print ("start training, used feature:", len(train_data[0]))

            clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
             min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
            #max_depth=2, random_state=0 entropy gini #min_samples_leaf=5
            clf.fit(train_data, train_label)  
            roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            print ("AUC", self.cohort_names[lack], roc_auc) 
            sample_total = len(train_label) + len(test_label)
            auc_total += roc_auc * len(test_label)
            auc_list.append(roc_auc)
            auc_data_frame.append([round(roc_auc, 3), self.cohort_names[lack], group])
        print ("weighted mean:", auc_total/sample_total)
        return auc_list

    def select_feature_3_2(self, train_cohort, select_feature_num):

        all_select_edges = {}

        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.all_data:
            if sample.cohort != train_cohort:
                continue
            if sample.disease == pheno_group_2:
                crc_num += 1
            if sample.disease == pheno_group_1:
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
                if sample.disease == pheno_group_2:
                    specific_HGT[edge][0] += 1
                if sample.disease == pheno_group_1:
                    specific_HGT[edge][1] += 1

        select_edges = {}
        for tag in specific_HGT:
            crc_ratio = float(specific_HGT[tag][0])/crc_num
            control_ratio = float(specific_HGT[tag][1])/control_num
            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
            U1, p = mannwhitneyu(crc_array, control_array, method="auto")
            select_edges[tag] = p 
        print ("from %s choose %s edges."%(len(specific_HGT), len(select_edges)), "filter best", best_num)
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)
        final_select_edges = {}
        if select_feature_num > len(sort_select_edges):
            select_feature_num = len(sort_select_edges)
        for i in range(select_feature_num):
            edge = sort_select_edges[i][0]
            final_select_edges[edge] = i
        return final_select_edges

    def save_one_out(self, select_feature_num):
        auc_list = []
        
        for lack in range(len(self.cohort_names)):
            # self.select_ranking(self.cohort_names[lack])
            # select_edges = self.get_top_feature(select_feature_num)

            select_edges = self.select_feature_3_2(self.cohort_names[lack], select_feature_num)
            print ("%s HGT events were selected."%(len(select_edges)))
            

            self.get_event_matrix(select_edges)
            for sample in self.all_data:
                sample.given_nodes_make_matrix(select_edges)

            print ("prepare data")
            train_data, train_label =  self.complex_data("test", self.cohort_names[lack])
            print ("start training, used feature:", len(train_data[0]))

            clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
             min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
            #max_depth=2, random_state=0 entropy gini #min_samples_leaf=5
            clf.fit(train_data, train_label)   

            single_list = []
            for j in range(len(self.cohort_names)):
                if j != lack:
                    test_data, test_label =  self.complex_data("test", self.cohort_names[j])
                    roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
                    single_list.append(roc_auc)
                    print (self.cohort_names[lack], self.cohort_names[j],roc_auc) 
            print ("mean", self.cohort_names[lack], np.mean(single_list))
            auc_list.append(np.mean(single_list))
            # break
        return auc_list

    def study2study(self, select_feature_num):
        all_study = []
        fold_num = 5
        for cohort in self.cohort_names:
            cohort_data = []
            auc_list = []
            for sample in self.all_data:
                if sample.cohort == cohort:
                    cohort_data.append(sample)
            random.Random(314).shuffle(cohort_data)
            each_fold = int(len(cohort_data)/fold_num)
            for f in range(fold_num):
                train_sam = cohort_data[:each_fold*f] + cohort_data[each_fold*f+each_fold:]
                test_sam = cohort_data[each_fold*f:each_fold*f+each_fold]

                specific_HGT = {} 
                crc_num = 0
                control_num = 0
                for sample in train_sam:
                    if sample.disease == pheno_group_2:
                        crc_num += 1
                    if sample.disease == pheno_group_1:
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
                        if sample.disease == pheno_group_2:
                            specific_HGT[edge][0] += 1
                        if sample.disease == pheno_group_1:
                            specific_HGT[edge][1] += 1
                select_edges = {}
                for tag in specific_HGT:
                    crc_ratio = float(specific_HGT[tag][0])/crc_num
                    control_ratio = float(specific_HGT[tag][1])/control_num
                    if max([crc_ratio, control_ratio]) <= self.common_ratio:
                        continue
                    crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
                    control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])

                    U1, p = mannwhitneyu(crc_array, control_array, method="auto")
                    select_edges[tag] = p 
                # print ("from %s choose %s edges."%(len(specific_HGT), len(select_edges)))
                sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)
                final_select_edges = {}

                if select_feature_num > len(sort_select_edges):
                    select_feature_num = len(sort_select_edges)
                for i in range(select_feature_num):
                    edge = sort_select_edges[i][0]
                    final_select_edges[edge] = i

                self.get_event_matrix(final_select_edges)
                for sample in train_sam:
                    sample.given_nodes_make_matrix(final_select_edges)
                for sample in test_sam:
                    sample.given_nodes_make_matrix(final_select_edges)
                train_data, train_label =  self.complex_data_study(train_sam)
                test_data, test_label = self.complex_data_study(test_sam) 

                # print (len(train_data), len(test_data))
                # print ("start training, used feature:", len(train_data[0]))
                clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
                    min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
                clf.fit(train_data, train_label)   
                roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
                # print ("AUC", roc_auc) 
                
                auc_list.append(roc_auc)
            print (cohort, np.mean(auc_list), auc_list)
            all_study.append(np.mean(auc_list))
            print ("*************************")
            # break
        print (select_feature_num, np.mean(all_study), all_study)
        return np.mean(all_study)

    def complex_data_study(self, cohort_data):
        data = []
        label = []
        for sample in cohort_data:
            sample_array = np.dot(self.event_array, sample.select_feature_array)
            # sample_array = list(sample_array) + sample.marker_abundance
            data.append(sample_array)
            if sample.disease == pheno_group_2:
                label.append(1)
            if sample.disease == pheno_group_1:
                label.append(0)            
            # label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

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
        graph = nx.from_numpy_matrix(event_array)   
        laplacian = nx.laplacian_matrix(graph)
        # laplacian = nx.normalized_laplacian_matrix(graph)
        laplacian = scipy.sparse.csr_matrix.toarray(laplacian)
        # print (laplacian)

        # laplacian = unnormalized_laplacian(event_array)
        # print (laplacian)

        a, b = np.linalg.eig(laplacian)
        self.event_array = b.T.real

    def complex_data_density(self, data_usage, select_cohort):
        data = []
        label = []
        for sample in self.all_data:
            if (sample.cohort == select_cohort and data_usage == "train"):
                continue
            if (sample.cohort != select_cohort and data_usage == "test"):
                continue
            # sample_array = sample.marker_abundance
            # sample_array = sample.select_feature_array
            sample_array = []
            for i in range(1, 7):
                density = sample.get_HGT_matrix(i)
                sample_array.append(density)
            # sample_array = np.dot(self.event_array, sample.select_feature_array)
            # sample_array = list(sample_array) + sample.marker_abundance
            if group == "Abundance":
                sample_array = sample.marker_abundance
            data.append(sample_array)
            if sample.disease == pheno_group_2:
                label.append(1)
            if sample.disease == pheno_group_1:
                label.append(0)            
            # label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

    def complex_data(self, data_usage, select_cohort):
        data = []
        label = []
        for sample in self.all_data:
            if (sample.cohort == select_cohort and data_usage == "train"):
                continue
            if (sample.cohort != select_cohort and data_usage == "test"):
                continue
            # sample_array = sample.marker_abundance
            # sample_array = sample.select_feature_array

            sample_array = np.dot(self.event_array, sample.select_feature_array)
            sample_array = list(sample_array) + sample.marker_abundance
            if group == "Abundance":
                sample_array = sample.marker_abundance
            data.append(sample_array)
            if sample.disease == pheno_group_2:
                label.append(1)
            if sample.disease == pheno_group_1:
                label.append(0)            
            # label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        return data, label

    def prepare_NN(self, select_feature_num):
        with open('filtered_HGTs.pkl', 'rb') as f:
        # with open('best_features.pkl', 'rb') as f:
            select_edges = pickle.load(f)
        print ("%s features selected."%(len(select_edges)))
        self.get_event_matrix(select_edges)
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
        feature_dict = {}
        lable_dict = {}
        for cohort in self.cohort_names:
            train_data, train_label =  self.complex_data("test", cohort)
            feature_dict[cohort] = train_data
            lable_dict[cohort] = train_label
            print (cohort, len(train_data), len(train_label))
        with open("feature_dict.pkl", 'wb') as f:
            pickle.dump(feature_dict, f)
        with open("lable_dict.pkl", 'wb') as f:
            pickle.dump(lable_dict, f)

    def __del__(self):
        del self.all_data
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
    TREE_NUM = 1000
    LEAF_S = 5
    pheno_group_1 = "CRC"   #"adenoma"
    pheno_group_2 = "control"  #"CRC"   "control"
    level, common_ratio, feature_num = 5, 0.06, 150   # 0.06, 150  0.2, 300 0.05, 300     0.1 700
    # group = "+ HGT"
    group = ""
    auc_data_frame = []

    # rf = RF(common_ratio)
    # for group in ["Abundance", "+HGT"]:
    #     print (group == "Abundance")
    #     auc_list = rf.LODO(26)
    # df = pd.DataFrame(auc_data_frame, columns = ["AUC", "Cohort", "Group"])
    # df.to_csv('for_AUC_plot.csv', sep='\t')
    # os.system("Rscript plot_auc.R")



    rf = RF(common_ratio)
    diff_cut = 1e-6
    feature_num = 26
    # for feature_num in range(10, 50, 2):
        # rf.select_feature_5(diff_cut)
    auc_list = rf.LODO(feature_num)
    # rf.prepare_NN(feature_num)
    print (diff_cut, feature_num, auc_list, np.mean(auc_list))
    print (diff_cut, feature_num, auc_list, np.mean(auc_list), file = result_file)





