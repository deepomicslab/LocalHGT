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
from scipy.special import comb

previous_marker_species = ["Peptostreptococcus stomatis", "Fusobacterium nucleatum", "Parvimonas spp.", "Porphyromonas asaccharolytica", "Gemella morbillorum",
"Clostridium symbiosum", "Parvimonas micra", "Escherichia coli", "Streptococcus parasanguinis", "Clostridium leptum", "Clostridium hathewayi",
"Anaerotruncus colihominis", "Prevotella copri", "Eisenbergiella tayi", "Actinomyces graevenitzii", "Alistipes spp."]

from ana_CRC_species import Sample, Taxonomy

def get_abd():
    marker_species_dict = {}
    i = 0
    for species in previous_marker_species:
        if species.split()[-1] == "spp.":
            form_name = "g__" + species.split()[0] 
        else:
            form_name = "s__" + "_".join(species.split()) 
        # print (form_name)
        marker_species_dict[form_name] = i 
        i += 1
    print (marker_species_dict)
    return marker_species_dict


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


class Vali_Sample(Sample):

    def __init__(self, bkp_file, ID):
        self.bkp_file = bkp_file
        self.bkps = []
        self.ID = ID
        self.tag = "normal"
        self.genus_abundance = ''
        self.reads_num = 0
        self.bases = 0
        self.basic_features = []
        self.level = 5
        # self.name = phenotype.ID_name[ID]
        
        # self.disease = ''
        # self.cohort = ''
        # self.marker_abundance = ''
        # self.disease = phenotype.ID_disease[ID]
        # self.cohort = phenotype.ID_cohort[ID]    
        # self.bases = phenotype.ID_bases[ID]   
        # self.basic_features =  phenotype.ID_basics[ID] 
        # self.marker_abundance = phenotype.ID_marker_abundance[ID]
  
        self.read_bkp()
        if len(self.bkps) == 0:
            self.tag = "no bkp"
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


class Validation():

    def __init__(self):
        self.status_dict = {} # SRA ID : status
        self.ID_name = {} # ID: sample_name
        self.sample_abd = {} # sample_name: [abundance list]
        pass

    def read_result_list(self, all_acc_file, cohort):
        validation_data = []
        for line in open(all_acc_file):
            acc_file = line.strip()
            ID = acc_file.split("/")[-1].split(".")[0]
            if ID not in self.status_dict:
                print (f"{ID} is not CRC.")
                continue
            if self.ID_name[ID] not in self.sample_abd:
                print (f"{ID} has no abundance file.")
                continue
            sample = Vali_Sample(acc_file, ID)

            sample.disease = self.status_dict[ID]
            sample.cohort = cohort
            sample.genus_abundance = self.sample_abd[self.ID_name[ID]]
            sample.given_nodes_make_matrix(select_edges)
            # print (sample.ID, sample.disease, sample.genus_abundance, sample.select_feature_array)
            validation_data.append([sample.ID, sample.disease, sample.genus_abundance, sample.select_feature_array])
        with open('validation_data.pkl', 'wb') as f:
            pickle.dump(validation_data, f) 

    def read_meta(self, meta_file):
        f = open(meta_file)
        for line in f:
            array = line.strip().split("\t")
            if array[0] != "SRP128485":
                continue
            sample_name = array[4]
            ID = array[5]
            status = array[8]
            if status not in status_conversion:
                continue
            status = status_conversion[status]
            # print (ID, status)
            self.status_dict[ID] = status
            self.ID_name[ID] = sample_name

    def read_abundance(self, marker_genus, abundance_dir):
        self.sample_abd = {}
        marker_species_num = len(marker_genus)
        i = 0
        for genus in marker_genus:
            marker_genus[genus] = i
            i += 1
        # print (marker_species_num, marker_genus)

        files = os.listdir(abundance_dir)
        for abd_file in files:
            sample = abd_file.split("_")[0]
            self.sample_abd[sample] = [0] * marker_species_num
            abd_file = os.path.join(abundance_dir, abd_file)
            # print (abd_file)
        # for cohort in cohort_abd:
        #     abd_file = cohort_abd[cohort]
            f = open(abd_file, 'r')
            i = 0
            for line in f:
                if line[0] == "#":
                    continue
                else:
                    array = line.strip().split()
                    # print (array)
                    abundance = float(array[2])
                    if len(array[0].split("|")) < 6:
                        continue
                    species_name = array[0].split("|")[-1]
                    genus_name = array[0].split("|")[-2]
                    if species_name in marker_genus:
                        # print ("name", species_name, marker_genus)
                        species_index = marker_genus[species_name]
                        self.sample_abd[sample][species_index] += abundance

                    elif genus_name in marker_genus:
                        # print ("name", species_name, marker_genus)
                        species_index = marker_genus[genus_name]
                        self.sample_abd[sample][species_index] += abundance
                i += 1
        # print (self.sample_abd)
        return self.sample_abd





def read_markers():

    # marker_genus = {"s__Peptostreptococcus_stomatis":1, "s__Fusobacterium_nucleatum":1}
    # select_edges = {"g__Porphyromonas&g__Porphyromonas":0}
    with open('marker_genus.pkl', 'rb') as f:
        marker_genus = pickle.load(f)
    with open('select_edges.pkl', 'rb') as f:
        select_edges = pickle.load(f)
    return marker_genus, select_edges

if __name__ == "__main__":

    status_conversion = {"Colorectal Neoplasms":"CRC", "Health":"control" }
    meta_file = "validation/last_gutmeta_sample.tsv"
    result_list = "validation/wenkui.list"
    abundance_dir = "validation/wenkui_metaphlan3/"

    marker_genus, select_edges = read_markers()
    taxonomy = Taxonomy()
    validation = Validation()
    validation.read_abundance(marker_genus, abundance_dir)
    validation.read_meta(meta_file)
    validation.read_result_list(result_list, "wenkui_yiqi")