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
# from KR_norm_juicer import KRnorm_sym
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

level_dict = {"phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
gender_dict = {"male":0, "female":1, "nan": 2}
# sra_meta = "italy.csv"
pheno_file = "/mnt/d/breakpoints/script/analysis/allmetadata.xlsx"#"CRC.xlsx"
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
"KarlssonFH_2013":"2021-10-14.KarlssonFH_2013.relative_abundance.xls",
"NielsenHB_2014":"2021-03-31.NielsenHB_2014.relative_abundance.xls",
"HallAB_2017":"2021-10-14.HallAB_2017.relative_abundance.xls",
"QinJ_2012":"2021-10-14.QinJ_2012.relative_abundance.xls",
"DavidLA_2015":"2021-03-31.DavidLA_2015.relative_abundance.xls",
"KieserS_2018":"2021-10-14.KieserS_2018.relative_abundance.xls"
}
marker_species = ["Peptostreptococcus stomatis", "Fusobacterium nucleatum", "Parvimonas spp.", "Porphyromonas asaccharolytica", "Gemella morbillorum",
"Clostridium symbiosum", "Parvimonas micra", "Escherichia coli", "Streptococcus parasanguinis", "Clostridium leptum", "Clostridium hathewayi",
"Anaerotruncus colihominis", "Prevotella copri", "Eisenbergiella tayi", "Actinomyces graevenitzii", "Alistipes spp."]

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
        save_file = "/mnt/d/breakpoints/script/analysis//taxonomy_dict.pkl"
        if not os.path.isfile(save_file):
            self.read_UHGG()
            with open(save_file, 'wb') as f:
                pickle.dump(self.taxonomy_dict, f)
        else:
            with open(save_file, 'rb') as f:
                self.taxonomy_dict = pickle.load(f)

level_list = ["phylum", "class", "order", "family", "genus", "species", "genome"]
bin_size = 100


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
        self.abundance = None
        self.split_abundance = None
        self.read = None

        self.from_ref_lineage = taxonomy.taxonomy_dict[self.from_ref_genome]
        self.to_ref_lineage = taxonomy.taxonomy_dict[self.to_ref_genome]

        taxa1 = self.from_ref_lineage.split(";")[level]
        taxa2 = self.to_ref_lineage.split(";")[level]
        if taxa1[1:] == "__" or taxa2[1:] == "__":
            self.hgt_tag = "NA"
        else:
            self.hgt_tag = "&".join(sorted([taxa1, taxa2]))
        # self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        
def get_genome_taxa(genome, level):
    # g1 = get_pure_genome(genome)
    taxa = taxonomy.taxonomy_dict[genome].split(";")[level]
    return taxa

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

    def __init__(self, bkp_list, ID, pheno, reads_num):
        self.bkps = bkp_list
        self.ID = ID  
        self.cohort = pheno[0]  
        self.disease = pheno[1]  
        self.full_disease = pheno[2].split(";")
        self.reads_num = reads_num
        self.select_feature_array = None
        self.genus_abundance = None
        self.thomas_abundance = None
        self.level = 5

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

class Data_load():

    def __init__(self):
        self.sample_obj_list = []
        self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b",\
            "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        # self.validate_cohort = ["KarlssonFH_2013", "NielsenHB_2014", "HallAB_2017", "QinJ_2012", "DavidLA_2015", "KieserS_2018"]
        self.validate_cohort = ["KarlssonFH_2013", "NielsenHB_2014", "HallAB_2017", "QinJ_2012"]
        self.validate_samples_list = []
        
    def read_samples(self):

        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        # os.system(f"ls {tgs_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        # os.system(f"ls {wenkui_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            if sra_id not in phenotype_dict:
                continue
            sample = Sample([], sra_id, phenotype_dict[sra_id], 0)
            if sample.cohort not in self.cohort_names and sample.cohort not in self.validate_cohort:
                continue
            if len(sample.full_disease) > 1:
                continue
            # if sample.disease not in ["control", "CRC"]:
            #     continue
            if sample.disease == '' or sample.disease == 'adenoma':
                continue

            my_bkps, reads_num = read_bkp(acc_file)
            sample = Sample(my_bkps, sra_id, phenotype_dict[sra_id], reads_num)
                # if sample.disease == '':
                #     sample.disease = 'control'
            if sample.cohort in self.cohort_names:
                self.sample_obj_list.append(sample)
                if len(self.sample_obj_list) % 100 == 0:
                    print ("read %s samples"%(len(self.sample_obj_list)))

            if sample.cohort in self.validate_cohort:
                self.validate_samples_list.append(sample)
                if len(self.validate_samples_list) % 100 == 0:
                    print ("read %s vali samples"%(len(self.validate_samples_list)))           
        print ("data is loaded.", "CRC cohort", len(self.sample_obj_list), "validate", len(self.validate_samples_list))
        
def read_bkp(bkp_file):
    my_bkps = []
    f = open(bkp_file)
    all_rows = csv.reader(f)
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
            my_bkps.append(eb)
    for eb in my_bkps:
        eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
    f.close()
    return my_bkps, reads_num

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
    return phenotype_dict

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
        f = open("/mnt/d/breakpoints/script/analysis/use/" + abd_file, 'r')
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
    file = "/mnt/d/breakpoints/script/analysis/last_gutmeta_sample.tsv"
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
        self.read_sra_meta("/mnt/d/HGT/CRC/IBD/ERP002061.txt")
        self.read_sra_meta("/mnt/d/HGT/CRC/IBD/HallAB_2017.txt")

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
                sample_name = sra_ID
            if sra_meta == "/mnt/d/HGT/CRC/IBD/ERP002061.txt":
                sra_ID = df["BioSample"][i] 
                sample_name = df["Sample Name"][i]   
                sample_name = sample_name.replace(".", "_")  
                sample_name = sample_name.replace("-", "_")  
            if sra_meta == "/mnt/d/HGT/CRC/IBD/HallAB_2017.txt":
                sample_name = df["Library Name"][i] 
                # print (sample_name)
                sample_name = sample_name.replace("mo", "_mo")  
                # print (sra_ID, sample_name)

            self.ID_name[sra_ID] = sample_name

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
        
            if cohort == "ZellerG_2014" or cohort == "YachidaS_2019" or cohort == "HanniganGD_2017":
                sample_name = row[2]
            if cohort == "KarlssonFH_2013":
                # for z in range(len(row)):
                #     print (z, row[z])
                sample_name = row[21]
                if len(sample_name.split(";")) > 1:
                    sample_name = sample_name.split(";")[0]
            if cohort == "QinJ_2012" or cohort == "DavidLA_2015" or cohort == "KieserS_2018":
                sra_ID = str(row[21])
                sample_name = row[3]
                # print (cohort, sra_ID, sample_name)
                self.ID_name[sra_ID] = sample_name

                # print (cohort, sample_name)
            # if cohort == "VogtmannE_2016":
            #     print (cohort, sample_name, condition, row[7])
            self.name_sample_id[sample_name] = sample_id
            # print (sample_name,condition) 
        # print (crc, control)

###############
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

class RF():

    def __init__(self, sample_obj_list):
        self.all_data = sample_obj_list
        self.disease_sample_num_cohort = {}
        self.cohorts_names = {}        
        self.HGT_network = ''
        self.species_network = ''
        self.cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b",\
            "WirbelJ_2018","ZellerG_2014","YuJ_2015"]
        print ("RF init done")

    def combine_markers(self, select_feature_num):
        out_marker = open("/mnt/d/breakpoints/script/analysis/CRC_prediction_markers.csv", 'w')
        select_HGTs, abun_related_HGTs = self.select_top_HGT(select_feature_num)
        select_edges = select_HGTs
        i = 0
        for edge in select_edges:
            select_edges[edge] = i
            print (edge.split("&")[0], edge.split("&")[1], sep = ",", file = out_marker)
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
        # print (len(marker_genus), "species for abundance", marker_genus)
        self.get_abundance_marker_value(marker_genus, self.all_data)
        self.get_thomas_abundance_value(self.all_data)
        self.get_genus_network(marker_genus)
        for genus in marker_genus:
            print (genus, '', sep = ",", file = out_marker)
        self.store_markers(marker_genus, select_edges)
        print ("HGT biomarker num", len(select_edges), "abun biomarker num", len(marker_genus))
        return select_edges, marker_genus

    def store_markers(self, marker_genus, select_edges):
        # store the markers
        with open('/mnt/d/breakpoints/script/analysis/marker_genus.pkl', 'wb') as f:
            pickle.dump(marker_genus, f) 
        with open('/mnt/d/breakpoints/script/analysis/select_edges.pkl', 'wb') as f:
            pickle.dump(select_edges, f) 

    def get_abundance_marker_value(self, marker_genus, sample_list):

        sample_abd = get_genus_abd(marker_genus)
        # print (sample_abd)
        for sample in sample_list:
            found_genus = False
            if sample.ID in sample_abd:
                genus_abundance = sample_abd[sample.ID]
                found_genus = True
            else:
                if sample.ID in phenotype.ID_name:
                    sample_name = phenotype.ID_name[sample.ID]
                    if sample_name in sample_abd:
                        genus_abundance = sample_abd[sample_name]
                        found_genus = True

                    elif sample_name in phenotype.name_sample_id:
                        sample_id = phenotype.name_sample_id[sample_name]
                        genus_abundance = sample_abd[sample_id]
                        found_genus = True

                if found_genus == False and sample.ID in phenotype.name_sample_id:
                    sample_name = phenotype.name_sample_id[sample.ID]
                    if sample_name in sample_abd:
                        genus_abundance = sample_abd[sample_name]
                        found_genus = True

            if found_genus == False:
                print ("# no abun", sample.ID, sample.cohort)
            sample.genus_abundance = genus_abundance

    def get_thomas_abundance_value(self, sample_list):

        # print (sample_abd)
        for sample in sample_list:
            found_genus = False
            if sample.ID in thomas_16_sample_abd:
                genus_abundance = thomas_16_sample_abd[sample.ID]
                found_genus = True
            else:
                if sample.ID in phenotype.ID_name:
                    sample_name = phenotype.ID_name[sample.ID]
                    if sample_name in thomas_16_sample_abd:
                        genus_abundance = thomas_16_sample_abd[sample_name]
                        found_genus = True

                    elif sample_name in phenotype.name_sample_id:
                        sample_id = phenotype.name_sample_id[sample_name]
                        genus_abundance = thomas_16_sample_abd[sample_id]
                        found_genus = True

                if found_genus == False and sample.ID in phenotype.name_sample_id:
                    sample_name = phenotype.name_sample_id[sample.ID]
                    if sample_name in thomas_16_sample_abd:
                        genus_abundance = thomas_16_sample_abd[sample_name]
                        found_genus = True

            if found_genus == False:
                print ("# no abun", sample.ID, sample.cohort)
            sample.thomas_abundance = genus_abundance

    def select_top_HGT(self, select_feature_num):
        specific_HGT = {} 
        abun_related_HGTs = {}
        record_all_HGTs = {}
        crc_num = 0
        control_num = 0
        edge_distribution = []
        for sample in self.all_data:
            if sample.disease == "CRC" :
                crc_num += 1
            elif sample.disease == "control" :
                control_num += 1
            else:
                print ( sample.disease, sample.full_disease)
                continue
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                support_ratio = bkp.cross_split_reads/sample.reads_num
                if support_ratio == 0:
                    continue
                if edge not in record_all_HGTs:
                    record_all_HGTs[edge] = 1
                    specific_HGT[edge] = [[], []]
                if edge not in sample_dict:
                    sample_dict[edge] = 0
                sample_dict[edge] += bkp.cross_split_reads
            edge_distribution.append(len(sample_dict))
            # in each sample, choose same number of edge.       
            sample_dict = self.same_HGT_number(sample_dict, max_hgt_num)

            for edge in sample_dict:
                if sample.disease == "CRC" :
                    specific_HGT[edge][0].append(sample_dict[edge])
                else:
                    specific_HGT[edge][1].append(sample_dict[edge])

        select_edges = {}
        genus_level_markers = {}
        for marker in marker_genus:
            if marker[0] == "s":
                marker = "g__" + marker.split("_")[2]
            genus_level_markers[marker] = 1
        print ("abundance marker-related genus:", len(genus_level_markers))

        diff_data = []
        for tag in specific_HGT:
            # if len(specific_HGT[tag][0]) + len(specific_HGT[tag][1]) < 25:
            #     continue
            array = tag.split("&")
            species_1 = array[0]
            species_2 = array[1]

            a = len(specific_HGT[tag][0])
            b = crc_num - len(specific_HGT[tag][0])
            c = len(specific_HGT[tag][1])
            d = control_num - len(specific_HGT[tag][1])
            # print (a, b, c, d, control_num, len(specific_HGT[tag][1]))

            oddsratio, p_value = fisher_exact([[a, b], [c, d]])
            diff_data.append([tag, p_value, oddsratio, species_1, species_2])

        df = pd.DataFrame(diff_data, columns = ["genus_pair", "p_value", "oddsratio", "species_1", "species_2"])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected
        filtered_df = df[df['p.adj'] < 0.05]
        # print (filtered_df)

        print ("diff HGT genums pair number is ", len(filtered_df))

        for index, row in filtered_df.iterrows():
            select_edges[row["genus_pair"]] = row["p.adj"]
            if row["p.adj"] < 0.05 and (row["species_1"] in genus_level_markers and row["species_2"] in genus_level_markers):
                abun_related_HGTs[row["genus_pair"]] = 1


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

    def same_HGT_number(self, sample_dict, choose_num): #only choose top x HGTs in each sample
        new_same_dict = {}
        sort_sample_dict = sorted(sample_dict.items(), key=lambda item: item[1], reverse = True)
        if len(sort_sample_dict) < choose_num:
            choose_num = len(sort_sample_dict)
        for z in range(choose_num):
            # new_same_dict[sort_sample_dict[z][0]] = sort_sample_dict[z][1]
            new_same_dict[sort_sample_dict[z][0]] = 1
        # print (sort_sample_dict[choose_num-1])
        return new_same_dict
   
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
            train_data, train_label =  self.complex_data("train", self.cohort_names[lack], self.all_data)
            test_data, test_label = self.complex_data("test", self.cohort_names[lack], self.all_data) 
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
        with open('/mnt/d/breakpoints/script/analysis/validation_data.pkl', 'rb') as f:
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
        with open('/mnt/d/breakpoints/script/analysis/selected_abun_genus_graph.pkl', 'wb') as f:
            pickle.dump([nodes_list, edges_list], f) 
        print ("abundance marker num:", len(marker_genus))
        pearson_network, pearson_edges_list = self.get_genus_matrix_pearson(marker_genus, id_conver)
        with open('/mnt/d/breakpoints/script/analysis/selected_abun_genus_graph_pearson.pkl', 'wb') as f:
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

    def complex_data(self, data_usage, select_cohort, sample_list):
        data = []
        label = []
        for sample in sample_list:
            if (sample.cohort == select_cohort and data_usage == "train"):
                continue
            if (sample.cohort != select_cohort and data_usage == "test"):
                continue
            sample_array_HGT = np.dot(self.HGT_network, sample.select_feature_array)
            sample_array_genus = np.dot(self.species_network, sample.genus_abundance)
            sample_array = list(sample.select_feature_array) + list(sample_array_genus) 

            if group == "Thomas-Abun":
                sample_array = sample.thomas_abundance
            elif group == "HGT":
                sample_array = list(sample.select_feature_array)
            elif group == "Abun":
                sample_array = list(sample.genus_abundance)
            elif group == "Abun-rectify":
                sample_array = list(sample_array_genus)
            # sample_array = list(sample_array_HGT)
            data.append(sample_array)

            if sample.disease == "CRC" :
                label.append(1)   
            else:
                 label.append(0)
        
            # label.append(sample.disease)
        data = np.array(data)
        label = np.array(label)
        print ("sample No.", len(data), len(label), "biomarker No.", len(data[0]))
        return data, label

    def t2d_validation(self, select_feature_num):
        select_edges,marker_genus = self.combine_markers(select_feature_num)
        # print (select_edges)
        
        for sample in self.all_data:
            sample.given_nodes_make_matrix(select_edges)
        train_data, train_label =  self.complex_data("train", "use all to train", self.all_data)

        validate_cohort = "KarlssonFH_2013"
        for validate_cohort in dat.validate_cohort:
            print (validate_cohort)
            T2D_data = []
            for sample in dat.validate_samples_list:
                if sample.cohort != validate_cohort:
                    continue
                sample.given_nodes_make_matrix(select_edges)
                T2D_data.append(sample)
            self.get_abundance_marker_value(marker_genus, T2D_data)
            self.get_thomas_abundance_value(T2D_data)
            test_data, test_label = self.complex_data("test", validate_cohort, dat.validate_samples_list)

            print ("validation sample num:", len(test_data))
            print ("validation CRC num", sum(test_label))
            print ("validation non-CRC num", len(test_label) - sum(test_label))
        
            clf = RandomForestClassifier(n_estimators=TREE_NUM, criterion="entropy", n_jobs = 10, \
                min_samples_leaf = LEAF_S, random_state = np.random.seed(2021)) 
            clf.fit(train_data, train_label)  
            # roc_auc = roc_auc_score(test_label, clf.predict_proba(test_data)[:,1])
            # print ("AUC", "for validation", roc_auc) 
            accuracy = sklearn.metrics.accuracy_score(test_label, clf.predict(test_data))
            print (validate_cohort, "accuracy", accuracy, "FPR", 1-accuracy )
        
        return [], 0

    def __del__(self):
        del self.all_data
        print ("object, deleted")

def extract_previous_16_markers(marker_genus, genus_abundance):
    marker_species_dict = get_abd() # previous 16 abundance markers
    previous_abundance = []
    for genus in marker_genus:
        index = marker_genus[genus]
        abundance = genus_abundance[index]
        if genus in marker_species_dict:
            previous_abundance.append(abundance)
    return previous_abundance

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


if __name__ == "__main__":
    level = 5
    TREE_NUM = 1000
    LEAF_S = 5
    NUM_FEATURE = "sqrt"
    feature_num, TREE_NUM, LEAF_S, NUM_FEATURE = 20,  1000, 5, 2
    delta, hgt_alpha, hgt_beta, hgt_gamma, eta, zeta, mu = 0, 0, 1, 1e-8, 1, 0.048, 0
    corr_para = 0.4

    feature_num = 18
    hgt_gamma = 1e-8
    zeta = 0.057

    group = ""
    auc_data_frame = []

    max_hgt_num = 2150  #2150
    i=2

    
    result_file = open("/mnt/d/breakpoints/script/analysis/random_forest.log", 'w')
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"
    abun_cutoff = 0  #1e-7
    level = 5
    phenotype_dict = read_phenotype()
    
    marker_genus = select_genus()
    thomas_16_sample_abd =  get_genus_abd(marker_genus)   
    phenotype = Phenotype() 
    taxonomy = Taxonomy()

    ### test the classifier in the independent CRC cohort and T2D cohort
    # group = "Thomas-Abun"
    group = "Hybrid"
    dat = Data_load()
    dat.read_samples()

    # rf = RF(dat.sample_obj_list)
    # rf.t2d_validation(feature_num) # independent T2D cohort
    # del rf

    # auc_list, weighted_mean = rf.LODO(feature_num)
    # must run python additional_validation.py 
    # auc_list, weighted_mean = rf.validation(feature_num) # independent CRC cohort


    ## compare the integration of HGT and abundance biomarkers and the previously reported 16 biomarkers
    # for group in ["Hybrid", "Thomas-Abun"]:  # for main plot
    for group in ["Hybrid", "Thomas-Abun", "HGT", "Abun"]:

        marker_genus = select_genus()
        phenotype = Phenotype()
        taxonomy = Taxonomy()
        rf = RF(dat.sample_obj_list)
        rf.t2d_validation(feature_num) # independent T2D cohort
        del rf

    # for feature_num in range(4, 30, 2):
        marker_genus = select_genus()
        phenotype = Phenotype()
        taxonomy = Taxonomy()
        rf = RF(dat.sample_obj_list)
    
        auc_list, weighted_mean = rf.LODO(feature_num)
        print (feature_num, group, auc_list, weighted_mean)
        del rf
    # df = pd.DataFrame(auc_data_frame, columns = ["AUC", "Cohort", "Group"])
    # df.to_csv('/mnt/d/R_script_files//for_AUC_plot.csv', sep='\t')
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
    # df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files/for_AUC_plot.csv', sep='\t')


