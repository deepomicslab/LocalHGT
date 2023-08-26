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

gender_dict = {"male":0, "female":1, "nan": 2}

class Phenotype():
    def __init__(self):
        self.name_full_disease = {}
        self.name_disease = {}
        self.name_cohort = {}
        self.name_bases = {}
        self.name_basics = {}
        self.ID_full_disease = {}
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
        self.read_sra_meta("/mnt/d/breakpoints/HGT/CRC/IBD/ERP002061.txt")  #IBD
        
    def read_sra_meta(self, sra_meta, for_name = "no"):
        df = pd.read_csv(sra_meta)
        for i in range(len(df.index)):
            sra_ID = df["Run"][i]
            if "sample_name" in df.columns and df["sample_name"][i] != "Illumina" and df["BioProject"][i] != "PRJDB4176" :
                sample_name = df["sample_name"][i]
            else:
                sample_name = df["Sample Name"][i]
            if sra_meta == "/mnt/d/breakpoints/HGT/CRC/austria/austria.csv":
                sample_name = "SID" + str(df["Submitter_Id"][i])
            if sra_meta == "/mnt/d/breakpoints/HGT/CRC/t2d.csv":
                sample_name = sra_ID
            if sra_meta == "/mnt/d/breakpoints/HGT/CRC/IBD/ERP002061.txt":
                
                sample_name = sample_name.replace(".", "_")
                sample_name = sample_name.replace("-", "_")
                sra_ID = df["BioSample"][i]
                # print (sample_name, sra_ID)
                # print (sample_name, sra_ID, self.name_disease[sample_name])
            study = df["SRA Study"][i]
            
            if sample_name not in self.name_disease:
                continue
            self.ID_name[sra_ID] = sample_name
            self.ID_disease[sra_ID] = self.name_disease[sample_name]
            self.ID_cohort[sra_ID] = self.name_cohort[sample_name]
            self.ID_bases[sra_ID] = self.name_bases[sample_name]
            self.ID_basics[sra_ID] = self.name_basics[sample_name]
            self.ID_full_disease[sra_ID] = self.name_full_disease[sample_name]
            
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
            sra_ID = str(row[21])
            if age == "nan":
                age = 0
            else:
                age = int(age)
            if BMI == "nan":
                BMI = 0
            else:
                BMI = round(float(BMI))
            if cohort == "ZellerG_2014" or cohort == "YachidaS_2019" or cohort == "HanniganGD_2017":
                sample_name = row[2]
            if cohort == "KarlssonFH_2013":
                sample_name = row[21]
                if len(sample_name.split(";")) > 1:
                    sample_name = sample_name.split(";")[0]
            if cohort == "NielsenHB_2014":
                sample_name = row[2]
                # print (sample_name)
            self.name_full_disease[sample_name] = full_disease
            self.name_disease[sample_name] = condition
            self.name_cohort[sample_name] = cohort
            self.name_basics[sample_name] = [age, gender, BMI]
            self.name_bases[sample_name] = bases
            self.name_sample_id[sample_name] = sample_id

            if sra_ID != "NA":
                self.ID_name[sra_ID] = sample_name
                self.ID_disease[sra_ID] = self.name_disease[sample_name]
                self.ID_cohort[sra_ID] = self.name_cohort[sample_name]
                self.ID_bases[sra_ID] = self.name_bases[sample_name]
                self.ID_basics[sra_ID] = self.name_basics[sample_name]
                self.ID_full_disease[sra_ID] = self.name_full_disease[sample_name]

            # print (sample_name,condition) 
        # print (crc, control)

class Sample():

    def __init__(self, ID):
        self.ID = ID
        self.tag = "normal"
        if ID in phenotype.ID_disease:
            self.disease = phenotype.ID_disease[ID]
            self.full_disease = phenotype.ID_full_disease[ID]
            self.cohort = phenotype.ID_cohort[ID]    
            self.bases = phenotype.ID_bases[ID]   
            self.basic_features =  phenotype.ID_basics[ID] 
            # self.marker_abundance = phenotype.ID_marker_abundance[ID]
        elif ID in phenotype.name_disease:
            self.disease = phenotype.name_disease[ID]
            self.full_disease = phenotype.name_full_disease[ID]
            self.cohort = phenotype.name_cohort[ID]    
            self.bases = phenotype.name_bases[ID]    
            self.basic_features = phenotype.name_basics[ID] 
            # self.marker_abundance = phenotype.name_marker_abundance[ID] 
        else:
            print ("## no pheno", ID) 
            self.tag = "no pheno"    

def find_sample_name_for_sra(ID_name, sra_meta):
    # for the sample whose SRA ID does not exist in meta table
    df = pd.read_csv(sra_meta)
    for i in range(len(df.index)):
        sra_ID = df["Run"][i]
        sample_name = sra_ID
        if sra_meta == "/mnt/d/breakpoints/HGT/CRC/austria/austria.csv":
            sample_name = "SID" + str(df["Submitter_Id"][i])
        elif "sample_name" in df.columns and df["sample_name"][i] != "Illumina" and df["BioProject"][i] != "PRJDB4176" :
            sample_name = df["sample_name"][i]
        else:
            sample_name = df["Sample Name"][i]
        # ID_name[sra_ID] = sample_name
        ID_name[sample_name] = sra_ID
    return ID_name

def get_ID_name():
    ID_name = {}
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/USA/usa.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/japan.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/yu_2015.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/germany.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/france/france.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/script/analysis/italy.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/austria/austria.csv")
    ID_name = find_sample_name_for_sra(ID_name, "/mnt/d/breakpoints/HGT/CRC/t2d.csv")
    return ID_name

def get_samples(hgt_result_dir):
    ID_dict = {}
    files = os.listdir(hgt_result_dir)
    for acc_file in files:
        if not re.search("acc.csv", acc_file):
            continue
        ID = acc_file.split(".")[0]
        ID_dict[ID] = 0
    return ID_dict

def get_pheno(ID_dict):     
    num = 0
    found_cohort = {}
    df = pd.read_excel(pheno_file, header=None) 
    for index, row in df.iterrows():
        if index == 0:
            continue
        cohort = row[1]
        find_name = False
        sample_name = "NA"
        for i in range(len(row)):
            if row[i] in ID_dict:
                find_name = True
                sample_name = row[i]
            # elif row[i] in ID_name and ID_name[row[i]] in ID_dict:
            #     find_name = True
            #     sample_name = ID_name[row[i]]

            else:
                array = str(row[i]).split(";")
                for arr in array:
                    if arr in ID_dict:
                        find_name = True
                        sample_name = arr
                        
        if find_name:
            ID_dict[sample_name] = 1
            if cohort not in found_cohort:
                found_cohort[cohort] = 1
            num += 1
    
    for ID in ID_dict:
        if ID_dict[ID] == 0:
            print ("not found", ID)
    print (len(ID_dict), num, found_cohort)


#### read tgs phenotype
def read_meta_tgs():
    
    sra_sample_dict = {}

    for line in open(tgs_meta):
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

def get_pheno_for_tgs(tgs_dir):
    sra_sample_dict = read_meta_tgs()
    add_data = []
    files = os.listdir(tgs_dir)
    for acc_file in files:
        if not re.search("acc.csv", acc_file):
            continue
        ID = acc_file.split(".")[0]
        if sra_sample_dict[ID][:2] == "TD":
            cohort = "Time-series"
        elif sra_sample_dict[ID][:2] == "CD":
            cohort = "Cross-sectional"
        else:
            print ("!!!!!!!!wrong")
        add_data.append([ID, cohort, "control", "healthy", "NA", "NA","NA","NA"])
    return add_data

### read wenkui CRC phenotype

def read_meta_wenkui(meta_file):
    status_dict = {}

    status_conversion = {"Colorectal Neoplasms":"CRC", "Health":"control" } #, "Irritable Bowel Syndrome":"IBS"
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
        status_dict[ID] = status
    return status_dict


def get_pheno_for_wenkui_CRC(wenkui_dir):
    status_dict = read_meta_wenkui(wenkui_meta_file)
    add_data = []
    files = os.listdir(wenkui_dir)
    cohort = "YangJ_2020"
    for acc_file in files:
        if not re.search("acc.csv", acc_file):
            continue
        ID = acc_file.split(".")[0]
        if ID not in status_dict:
            continue
        disease = status_dict[ID]
        if disease == "control":
            full_disease = "healthy"
        else:
            full_disease = "CRC"
        add_data.append([ID, cohort, disease, full_disease, "NA", "NA","NA","NA"])
    return add_data

if __name__ == "__main__":

    pheno_file = "/mnt/d/breakpoints/script/analysis/allmetadata.xlsx"
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    
    pheno_result = "/mnt/d/HGT/association/phenotype.csv"
    tgs_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    tgs_meta = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    wenkui_dir = "/mnt/d/breakpoints/HGT/CRC/wenkui/"
    wenkui_meta_file = "/mnt/d/breakpoints/script/analysis/validation/last_gutmeta_sample.tsv"

    data = []
    found_cohort = {}
    ID_dict = get_samples(hgt_result_dir)
    phenotype = Phenotype()
    num = 0
    for ID in ID_dict:   
        sample = Sample(ID)
        if sample.tag != "no pheno":
            num += 1
            data.append([sample.ID, sample.cohort, sample.disease, sample.full_disease, sample.bases] + sample.basic_features )
            if sample.cohort not in found_cohort:
                found_cohort[sample.cohort] = 1
    print (num, len(found_cohort), found_cohort)

    add_data = get_pheno_for_tgs(tgs_dir)
    data += add_data
    add_data_2 = get_pheno_for_wenkui_CRC(wenkui_dir)
    data += add_data_2
    df = pd.DataFrame(data, columns = ["sample", "cohort", "disease", "full_disease", "bases", "age", "gender",  "BMI"])
    df.to_csv(pheno_result, sep=',')

