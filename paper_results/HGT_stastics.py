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
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
import pydot
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
from scipy.linalg import svd
np.set_printoptions(threshold=sys.maxsize)
# from math import comb
import scipy.special as sc
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA

from ana_CRC_species import Sample, Acc_Bkp, Taxonomy, get_split_reads_cutoff

level = 5
window = 100

COG_categories = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z"]
COG_annotation = ["A: RNA processing and modification", "B: Chromatin structure and dynamics", "C: Energy production and conversion", "D: Cell cycle control, cell division, chromosome partitioning", "E: Amino acid transport and metabolism", "F: Nucleotide transport and metabolism", "G: Carbohydrate transport and metabolism", "H: Coenzyme transport and metabolism", "I: Lipid transport and metabolism", "J: Translation, ribosomal structure and biogenesis", "K: Transcription", "L: Replication, recombination and repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: Posttranslational modification, protein turnover, chaperones", "P: Inorganic ion transport and metabolism", "Q: Secondary metabolites biosynthesis, transport and catabolism", "R: General function prediction only", "S: Function unknown", "T: Signal transduction mechanisms", "U: Intracellular trafficking, secretion, and vesicular transport", "V: Defense mechanisms", "W: Extracellular structures", "Y: Nuclear structure", "Z: Cytoskeleton"]
level_dict = {"phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
level_list = ["phylum", "class", "order", "family", "genus", "species", "genome"]
cohort_names = ["YachidaS_2019","FengQ_2015","VogtmannE_2016","ThomasAM_2018a","ThomasAM_2018b",\
"WirbelJ_2018","ZellerG_2014","YuJ_2015"]

def get_tag(bkp, level):
    
    if level == 7:
        from_tax = bkp.from_ref
        to_tax = bkp.to_ref
    elif level == 8:
        from_tax = bkp.from_ref + "|" + str(round(bkp.from_bkp/window))
        to_tax = bkp.to_ref + "|" + str(round(bkp.to_bkp/window))
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

def obtain_all_feature(all_data):
    specific_HGT = {} 
    crc_num = 0
    control_num = 0
    cohort_sam_num = {}
    cohort_sam_HGT = {}
    for sample in all_data:
        if sample.cohort not in cohort_sam_num:
            cohort_sam_num[sample.cohort] = [0, 0]
            cohort_sam_HGT[sample.cohort] = {}
        if sample.disease == "control" :
            crc_num += 1
            cohort_sam_num[sample.cohort][0] += 1
        if sample.disease == "CRC" :
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
            if edge not in cohort_sam_HGT[sample.cohort]:
                cohort_sam_HGT[sample.cohort][edge] = [0, 0]
            if edge not in specific_HGT:
                specific_HGT[edge] = [0, 0]
            if sample.disease == "control" :
                specific_HGT[edge][0] += 1
                cohort_sam_HGT[sample.cohort][edge][0] += 1 
            if sample.disease == "CRC" :
                specific_HGT[edge][1] += 1
                cohort_sam_HGT[sample.cohort][edge][1] += 1 
    print ("HGT event num:", len(specific_HGT))
    return specific_HGT, cohort_sam_HGT

def get_data():
    with open("sample_data", "rb") as fp:
        all_data = pickle.load(fp)
    sample_num = len(all_data)
    print ("sample num:", sample_num)
    return all_data, sample_num

def filter_feature(specific_HGT, cohort_sam_HGT):
    filtered_HGT = {}
    HGT_index = 0
    for HGT in specific_HGT:
        existing_sample_num =  specific_HGT[HGT][0] + specific_HGT[HGT][1]
        if existing_sample_num < 25 or abs(specific_HGT[HGT][0] - specific_HGT[HGT][1]) < 2: # 25 2
            continue
        if max(specific_HGT[HGT]) < 22: # 22
            continue
        array = HGT.split("&")
        species_1 = array[0]
        species_2 = array[1]
        filtered_HGT[HGT] = HGT_index 
        HGT_index += 1
    return filtered_HGT

def tree_to_newick(g, root=None):
    # refer to https://stackoverflow.com/questions/46444454/save-networkx-tree-in-newick-format
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, g.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    for child in g[root]:
        if len(g[child]) > 0:
            subgs.append(tree_to_newick(g, root=child))
        else:
            subgs.append(child)
    return "(" + ','.join(subgs) + ")"

class Network():

    def __init__(self, all_data):
        self.data = all_data

    def compare_network(self):
        properties = ['density', 'transitivity', 'algebraic_connectivity', 'assortativity', 'Node', "Edge"]
        data = []
        # edge_num_list = [10, 10, 30, 40, 100, 120]
        edge_num_list = [17, 24, 21, 63, 88, 460]
        for level in range(1, 7):
            edge_num = edge_num_list[level-1]
            cohort_sam_num = {}
            cohort_base = {}
            properties_dict = {}
            for pro in properties:
                properties_dict[pro] = [[],[]]
            for sample in self.data:
                pro_list, total_edge_num = sample.get_network_properties(level, edge_num)
                if total_edge_num < edge_num:
                    continue
                if sample.disease == "control" :
                    index = 0
                else:
                    index = 1
                for i in range(len(properties)):
                    value = pro_list[i]#/sample.bases
                    origin = pro_list[i]
                    properties_dict[properties[i]][index].append(value)
                    data.append([properties[i], value, sample.disease, sample.cohort, level_list[level-1], origin, sample.bases])
            for i in range(len(properties)):
                # U1, p = mannwhitneyu(properties_dict[properties[i]][0], properties_dict[properties[i]][1])
                U1, p = scipy.stats.ranksums(properties_dict[properties[i]][0], properties_dict[properties[i]][1])
                print (len(properties_dict[properties[i]][0]), len(properties_dict[properties[i]][1]), level_list[level-1], properties[i], p, np.mean(properties_dict[properties[i]][0]), np.mean(properties_dict[properties[i]][1]), sep = "\t")
        df = pd.DataFrame(data, columns = ["Property", "Value", "Group", "Cohort", "Level", "Origin", "Bases"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files/network_comparison_normalized.csv', sep=',')
        # print (df)

    def compare_network_test(self):
        properties = ['density', 'transitivity', 'algebraic_connectivity', 'assortativity', 'Node', "Edge"]
        f = open("significant_net_compare.txt", 'w')
        data = []
        # 10 24 x 70
        edge_num_list = [[5, 21], [5, 31], [20, 61], [5, 101], [5, 201], [210, 500]]
        for level in range(3, 7):
            edge_num_interval = edge_num_list[level-1]
            for edge_num in range(edge_num_interval[0], edge_num_interval[1], 1):
                properties_dict = {}
                for pro in properties:
                    properties_dict[pro] = [[],[]]

                cohort_sam_num = {}
                cohort_base = {}
                for sample in self.data:
                    pro_list, total_edge_num = sample.get_network_properties(level, edge_num)
                    if total_edge_num < edge_num:
                        continue
                    if sample.disease == "control" :
                        index = 0
                    else:
                        index = 1
                    for i in range(len(properties)):
                        value = pro_list[i]#/sample.bases
                        origin = pro_list[i]
                        properties_dict[properties[i]][index].append(value)
                        data.append([properties[i], value, sample.disease, sample.cohort, level_list[level-1], origin, sample.bases])
                significant_num = 0
                for i in range(len(properties)):
                    # U1, p = mannwhitneyu(properties_dict[properties[i]][0], properties_dict[properties[i]][1])
                    U1, p = scipy.stats.ranksums(properties_dict[properties[i]][0], properties_dict[properties[i]][1])
                    if p < 0.05:
                        significant_num += 1
                    # print ("control num", len(properties_dict[properties[i]][0]), "CRC num", len(properties_dict[properties[i]][1]), level_list[level-1], properties[i], p, sep = "\t")
                    print ("control num",len(properties_dict[properties[i]][0]), "CRC num", len(properties_dict[properties[i]][1]), level_list[level-1], properties[i], p, np.mean(properties_dict[properties[i]][0]), np.mean(properties_dict[properties[i]][1]), sep = "\t")
                    print ("control num", len(properties_dict[properties[i]][0]), "CRC num", len(properties_dict[properties[i]][1]), level_list[level-1], properties[i], p, sep = "\t", file=f)
                print ("###", edge_num, level_list[level-1], significant_num)
                print ("###", edge_num, level_list[level-1], significant_num, file = f)
        f.close()
        # df = pd.DataFrame(data, columns = ["Property", "Value", "Group", "Cohort", "Level", "Origin", "Bases"])
        # df.to_csv('/mnt/c/Users/swang66/Documents/network_comparison_normalized.csv', sep=',')
        # print (df)

    def infer_sale(self):
        data = []
        # edge_num_list = [10, 10, 30, 40, 100, 120]
        f = open("scale_free_count.txt", 'w')
        edge_num_list = [17, 24, 21, 63, 88, 460]
        for level in range(1, 7):
            edge_num = edge_num_list[level-1]
            network_num = 0
            scale_free_num = 0
            for sample in self.data:
                p1, p2, p3, total_edge_num = sample.judge_scale_free(level, edge_num)
                if total_edge_num < edge_num:
                    continue
                level_name = level_list[level-1]
                data.append([p1, "lognormal_positive", sample.disease, sample.cohort, level_name])
                data.append([p2, "exponential", sample.disease, sample.cohort, level_name])
                data.append([p3, "Weibull", sample.disease, sample.cohort, level_name])
                network_num += 1
                if p1 >0 and p2>0 and p3>0:
                    scale_free_num += 1
            print (level, level_list[level-1], scale_free_num, network_num, scale_free_num/network_num, file = f)
        f.close()
        df = pd.DataFrame(data, columns = ["ratio", "Comparison", "Group", "Cohort", "Level"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//scale_free.csv', sep=',')

class Analyze():

    def __init__(self, level, all_data):
        self.data = all_data
        self.level = level
        self.disease = ["CRC", "control", "adenoma"]
        self.disease_index = {"CRC":0, "control":1, "adenoma":2}
        self.disease_sample_num = {"CRC":0, "control":0, "adenoma":0}
        self.disease_sample_num_cohort = {}
        self.window = 100

    def taxonomy_count(self):  # the taxa with most HGTs
        data = []
        for level in range(1, 7):
            count_dict = {}
            groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
            for sample in self.data:
                for bkp in sample.bkps:
                    from_tax = bkp.from_ref_lineage.split(";")[level]
                    to_tax = bkp.to_ref_lineage.split(";")[level]
                    for tax in [from_tax, to_tax]:
                        if len(tax) == 3:
                            continue
                        if tax not in count_dict:
                            count_dict[tax] = 1
                        else:
                            count_dict[tax] += 1
                        if tax not in groups_cout_dict[sample.disease]:
                            groups_cout_dict[sample.disease][tax] = 1
                        else:
                            groups_cout_dict[sample.disease][tax] += 1
            # for disease in groups_cout_dict.keys():
            #     print (disease, "Single Top")
            #     self.show_sort_dict(groups_cout_dict[disease])
            # self.show_sort_dict(count_dict)
            
            sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
            select_num = 20
            if len(sorted_count) < select_num:
                select_num = len(sorted_count)
            for i in range(select_num):
                print (sorted_count[i])
                data.append(list(sorted_count[i]) + [level] )
        df = pd.DataFrame(data, columns = ["Order", "The Number of HGT", "Level"])
        df.to_csv('order_sort.csv', sep='\t')
        print (df)
        os.system("Rscript order_sort.R")
        print ("done")
        # order_sort_count.py compare the number and ratio

    def taxonomy_count_cohort(self):  # the taxa with most HGTs
        data = []
        for cohort in cohort_names:
            for level in range(1, 7):
                print (cohort, level_list[level-1])
                count_dict = {}
                groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
                for sample in self.data:
                    if sample.cohort != cohort:
                        continue
                    for bkp in sample.bkps:
                        from_tax = bkp.from_ref_lineage.split(";")[level]
                        to_tax = bkp.to_ref_lineage.split(";")[level]
                        for tax in [from_tax, to_tax]:
                            if len(tax) == 3:
                                continue
                            if tax not in count_dict:
                                count_dict[tax] = 1
                            else:
                                count_dict[tax] += 1
                            if tax not in groups_cout_dict[sample.disease]:
                                groups_cout_dict[sample.disease][tax] = 1
                            else:
                                groups_cout_dict[sample.disease][tax] += 1
                
                sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
                select_num = 20
                if len(sorted_count) < select_num:
                    select_num = len(sorted_count)
                for i in range(select_num):
                    # print (sorted_count[i])
                    data.append(list(sorted_count[i]) + [level_list[level-1]] + [cohort] )
        df = pd.DataFrame(data, columns = ["Taxa", "Number", "Level", "Cohort"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//taxa_sort_cohort.csv', sep=',')

    def taxonomy_count_group(self):  # the taxa with most HGTs
        data = []
        for group in ["control", "CRC"]:
            for level in range(1, 7):
                print (group, level_list[level-1])
                count_dict = {}
                groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
                for sample in self.data:
                    if sample.disease != group:
                        continue
                    for bkp in sample.bkps:
                        from_tax = bkp.from_ref_lineage.split(";")[level]
                        to_tax = bkp.to_ref_lineage.split(";")[level]
                        for tax in [from_tax, to_tax]:
                            if len(tax) == 3:
                                continue
                            if tax not in count_dict:
                                count_dict[tax] = 1
                            else:
                                count_dict[tax] += 1
                            if tax not in groups_cout_dict[sample.disease]:
                                groups_cout_dict[sample.disease][tax] = 1
                            else:
                                groups_cout_dict[sample.disease][tax] += 1
                
                sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)
                select_num = 20
                if len(sorted_count) < select_num:
                    select_num = len(sorted_count)
                for i in range(select_num):
                    # print (sorted_count[i])
                    data.append(list(sorted_count[i]) + [level_list[level-1]] + [group] )
        df = pd.DataFrame(data, columns = ["Taxa", "Number", "Level", "Group"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//taxa_sort_group.csv', sep=',')

    def taxonomy_pair_count(self):
        level_dict = {"phylum":1, "genus":5, "species":6}
        level = 5
        count_dict = {}
        groups_cout_dict = {"CRC":{}, "control":{}, "adenoma":{}}
        for sample in self.data:
            # if sample.disease != "control":   #different groups
            #     continue
            for bkp in sample.bkps:
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
            for bkp in sample.bkps:
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
        num = 20
        if len(count_dict) < num:
            num = len(count_dict)
        for i in range(num):
            print (sorted_count[i])
        print ("**************")
        return sorted_count
    
    def check_CRC_related_species(self):
        CRC_related_species = ["Fusobacterium nucleatum", "Solobacterium moorei",\
         "Porphyromonas asaccharolytica","Parvimonas micra","Peptostreptococcus stomatis","Parvimonas spp"]
        for sample in self.data:
            for bkp in sample.bkps:
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
                control.append(len(sample.bkps))
            elif sample.disease == "CRC":
                crc.append(len(sample.bkps))
            elif sample.disease == "adenoma":
                adenoma.append(len(sample.bkps))
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

    def per_species_HGT_comparison(self): # HGT network density
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
            for bkp in sample.bkps:
                from_tax = bkp.from_ref_lineage.split(";")[self.level]
                to_tax = bkp.to_ref_lineage.split(";")[self.level]
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
            for bkp in sample.bkps:
                from_tax = bkp.from_ref_lineage.split(";")[self.level]
                to_tax = bkp.to_ref_lineage.split(";")[self.level]
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
            for bkp in sample.bkps:
                from_tax = bkp.from_ref_lineage.split(";")[self.level]
                to_tax = bkp.to_ref_lineage.split(";")[self.level]
                for tax in [from_tax, to_tax]:
                    
                    if tax not in index_dict:
                        sample_list[-1] += 1
                    else:
                        sample_list[index_dict[tax]] += 1
            sample_array = np.array(sample_list)/sum(sample_list)
            data.append(sample_array)
            # print (sample_list)
        df = pd.DataFrame(data, columns = class_list)
        print (df)
        df.plot(kind='bar', stacked=True,
        title='Stacked Bar Graph by dataframe', figsize=(10,3))
        plt.savefig('bar.pdf')

    def taxonomy_circos(self): #compare intra and inter proportions, and simple circos 
        data = []
        total = []
        for self.level in range(1, 7):
            count_dict = {}
            same_dict = {}
            total_diff = 0
            total_same = 0
            for sample in self.data:
                for bkp in sample.bkps:
                    from_tax = bkp.from_ref_lineage.split(";")[self.level]
                    to_tax = bkp.to_ref_lineage.split(";")[self.level]
                    if from_tax == to_tax:
                        if from_tax not in same_dict:
                            same_dict[from_tax] = 2
                        else:
                            same_dict[from_tax] += 2
                        continue
                    for tax in [from_tax, to_tax]:
                        if tax not in count_dict:
                            count_dict[tax] = 1
                        else:
                            count_dict[tax] += 1
            sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)

            num = 30 
            if num > len(sorted_count):
                num = len(sorted_count)
            for i in range(len(sorted_count)):
                if len(sorted_count[i][0]) == 3:
                    continue
                diff = sorted_count[i][1]
                if sorted_count[i][0] in same_dict:
                    same = same_dict[sorted_count[i][0]] 
                else:
                    same = 0
                array = np.array([diff, same])
                array = array/sum(array)
                total_diff += diff
                total_same += same
                if i < num:
                    print (self.level, sorted_count[i][0][3:], diff, same, array)  
                    data.append([array[0], sorted_count[i][0][3:]]+ ["Inter", self.level])
                    data.append([array[1], sorted_count[i][0][3:]] + ["Intra", self.level])
            array = np.array([total_diff, total_same])
            array = array/sum(array)
            print ("total", array)
            total.append([array[0], level_list[self.level-1], "Inter"])
            total.append([array[1], level_list[self.level-1], "Intra"])

        df = pd.DataFrame(data, columns = ["Proportion", "Order", "Group", "level"])
        df.to_csv('/mnt/c/Users/swang66/Documents/intra_order.csv', sep=',')
        df = pd.DataFrame(total, columns = ["Proportion", "Level", "Group"])
        df.to_csv('/mnt/c/Users/swang66/Documents/intra_order_total.csv', sep=',')
        # ax = df.plot(x="Order", kind='bar', stacked=True, figsize=(10,7), rot=90, fontsize=4)
        # ax.set_ylabel("Proportion")
        # ax.set_xlabel("")

        # plt.legend(loc='lower right')
        # plt.savefig('bar_inter_order.pdf')

        """
        select_num = 5
        class_list = [x[0] for x in sorted_count][:select_num] + ['others']     #list(count_dict.keys())
        index_dict = {}
        for i in range(len(class_list)):
            index_dict[class_list[i]] = i
        for co in sorted_count:
            if co[0] not in index_dict:
                index_dict[co[0]] = len(class_list)-1
        # print (index_dict)

        head_name = [x[0][3:] for x in sorted_count][:select_num] + ['others']
        data = []
        for i in range(len(class_list)):
            one_class = [0]*len(class_list)
            for sample in self.data:
                # if sample.disease != "CRC":
                #     continue
                
                for bkp in sample.bkps:
                    from_tax = bkp.from_ref_lineage.split(";")[self.level]
                    to_tax = bkp.to_ref_lineage.split(";")[self.level]
                    if index_dict[from_tax] == i:
                        j = index_dict[to_tax]
                        if j < i:
                            one_class[j] += 1
                    if index_dict[to_tax] == i:
                        j = index_dict[from_tax]
                        if j < i:
                            one_class[j] += 1
            data.append(one_class)
        df = pd.DataFrame(data, columns = head_name)
        df.index = head_name
        df.to_csv('for_circos.csv', sep='\t')
        print (df)
        os.system("Rscript circos.R")
        """

    def taxonomy_circos_CRC(self): #compare intra and inter proportions, and simple circos 
        data = []
        total = []
        for group in ["control", "CRC"]:
            for self.level in range(1, 7):
                count_dict = {}
                same_dict = {}
                total_diff = 0
                total_same = 0
                for sample in self.data:
                    if sample.disease != group:
                        continue
                    for bkp in sample.bkps:
                        from_tax = bkp.from_ref_lineage.split(";")[self.level]
                        to_tax = bkp.to_ref_lineage.split(";")[self.level]
                        if from_tax == to_tax:
                            if from_tax not in same_dict:
                                same_dict[from_tax] = 2
                            else:
                                same_dict[from_tax] += 2
                            continue
                        for tax in [from_tax, to_tax]:
                            if tax not in count_dict:
                                count_dict[tax] = 1
                            else:
                                count_dict[tax] += 1
                sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)

                num = 30 
                if num > len(sorted_count):
                    num = len(sorted_count)
                for i in range(len(sorted_count)):
                    if len(sorted_count[i][0]) == 3:
                        continue
                    diff = sorted_count[i][1]
                    if sorted_count[i][0] in same_dict:
                        same = same_dict[sorted_count[i][0]] 
                    else:
                        same = 0
                    array = np.array([diff, same])
                    array = array/sum(array)
                    total_diff += diff
                    total_same += same
                    # if i < num:
                    #     print (self.level, sorted_count[i][0][3:], diff, same, array)  
                    #     data.append([array[0], sorted_count[i][0][3:]]+ ["Inter", self.level])
                    #     data.append([array[1], sorted_count[i][0][3:]] + ["Intra", self.level])
                array = np.array([total_diff, total_same])
                array = array/sum(array)
                print ("total", array)
                total.append([array[0], level_list[self.level-1], "Inter", group])
                total.append([array[1], level_list[self.level-1], "Intra", group])

        # df = pd.DataFrame(data, columns = ["Proportion", "Order", "Group"])
        # df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files/intra_order_group.csv', sep=',')
        df = pd.DataFrame(total, columns = ["Proportion", "Level", "Type", "Group"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//intra_order_total_crc.csv', sep=',')

    def taxonomy_circos_cohort(self): #compare intra and inter proportions, and simple circos 
        data = []
        total = []
        shuffle(cohort_names)
        for cohort in cohort_names:
            for self.level in range(1, 7):
                count_dict = {}
                diff_dict = {}
                same_dict = {}
                total_diff = 0
                total_same = 0
                bkp_num = 0
                for sample in self.data:
                    if sample.cohort != cohort:
                        continue
                    for bkp in sample.bkps:
                        bkp_num += 1
                        from_tax = bkp.from_ref_lineage.split(";")[self.level]
                        to_tax = bkp.to_ref_lineage.split(";")[self.level]
                        if from_tax == to_tax:
                            total_same += 1
                        else:
                            total_diff += 1
                array = np.array([total_diff, total_same])
                array = array/sum(array)
                print ("total", cohort, level_list[self.level-1], array, bkp_num)
                total.append([array[0], level_list[self.level-1], "Inter", cohort])
                total.append([array[1], level_list[self.level-1], "Intra", cohort])

        df = pd.DataFrame(total, columns = ["Proportion", "Level", "Type", "Cohort"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//intra_order_total_cohort.csv', sep=',')

    def bkp_num_histogram_cohort(self): #compare intra and inter proportions, and simple circos 
        total = []
        all_bkp_num = []
        for cohort in cohort_names:
            for sample in self.data:
                if sample.cohort != cohort:
                    continue
                bkp_num = 0
                for bkp in sample.bkps:
                    bkp_num += 1
                all_bkp_num.append(bkp_num)
                total.append([bkp_num, cohort, sample.ID])

        df = pd.DataFrame(total, columns = ["Number", "Cohort", "ID"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//bkp_num_cohort.csv', sep=',')
        print (np.mean(all_bkp_num), np.median(all_bkp_num))
    
    def output_samples(self):
        data = []
        for sample in self.data:
            data.append([sample.ID, sample.disease, sample.cohort, sample.bases])
        df = pd.DataFrame(data, columns = ["ID", "Disease", "Cohort", "Bases"])
        df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//sample_info.csv', sep=',')

    def inter_taxa_circos(self): # see inter-taxa HGT count
        select_num = 5
        for self.level in range(1, 7):
            count_dict = {}
            same_dict = {}
            total_diff = 0
            total_same = 0
            for sample in self.data:
                for bkp in sample.bkps:
                    from_tax = bkp.from_ref_lineage.split(";")[self.level]
                    to_tax = bkp.to_ref_lineage.split(";")[self.level]
                    if from_tax == to_tax:
                        if from_tax not in same_dict:
                            same_dict[from_tax] = 2
                        else:
                            same_dict[from_tax] += 2
                        continue
                    for tax in [from_tax, to_tax]:
                        if tax not in count_dict:
                            count_dict[tax] = 1
                        else:
                            count_dict[tax] += 1
            sorted_count = sorted(count_dict.items(), key=lambda item: item[1], reverse = True)

            class_list = [x[0] for x in sorted_count][:select_num] + ['Others']     #list(count_dict.keys())
            index_dict = {}
            for i in range(len(class_list)):
                index_dict[class_list[i]] = i
            for co in sorted_count:
                if co[0] not in index_dict:
                    index_dict[co[0]] = len(class_list)-1

            head_name = [x[0][3:] for x in sorted_count][:select_num] + ['Others']
            data = []
            for i in range(len(class_list)):
                one_class = [0]*len(class_list)
                for sample in self.data:
                    for bkp in sample.bkps:
                        from_tax = bkp.from_ref_lineage.split(";")[self.level]
                        to_tax = bkp.to_ref_lineage.split(";")[self.level]
                        if index_dict[from_tax] == i:
                            j = index_dict[to_tax]
                            if j < i:
                                one_class[j] += 1
                        if index_dict[to_tax] == i:
                            j = index_dict[from_tax]
                            if j < i:
                                one_class[j] += 1
                data.append(one_class)
            df = pd.DataFrame(data, columns = head_name)
            df.index = head_name
            df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files/inter_taxa_files/inter_taxa_%s.csv'%(level_list[self.level-1]), sep='\t')
            print (level_list[self.level-1], "done")

    def cal_average_bkp_num(self):
        count_list = []
        for sample in self.data:
            count_list.append(len(sample.bkps))
        print (sum(count_list), np.mean(count_list))

    def plot_most_common_pair(self):
        count_dict = {}
        record_pos = {}

        for sample in self.data:
            sample_dict = {}
            for bkp in sample.bkps:
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

        remove_noise = {}
        for tag in record_pos:
            remove_noise[tag] = []
            for ele in record_pos[tag].copy():
                if ele[3] > 20:
                    remove_noise[tag].append(ele)
                
        for top_index in range(10):
            df = pd.DataFrame(record_pos[sorted_count_dict[top_index][0]], columns =sorted_count_dict[top_index][0].split("&")+["Condition", "Sample_Num"] )
            df2 = pd.DataFrame(remove_noise[sorted_count_dict[top_index][0]], columns =sorted_count_dict[top_index][0].split("&")+["Condition", "Sample_Num"] )
            # df = pd.DataFrame(record_pos[sorted_count_dict[top_index][0]], columns =["G1", "G2", "Condition", "Sample_Num"] )
            # print (df)
            print (top_index)
            df.to_csv('link_figure/for_para_%s.csv'%(top_index), sep='\t')
            df2.to_csv('link_figure/for_para_line_%s.csv'%(top_index), sep='\t')
            # os.system("Rscript parallel.R")
            # os.system("mv para.pdf link_figure/para_%s.pdf"%(top_index))
        #     # break
        return sorted_count_dict, record_pos

    def count_genome_density(self):
        count_dict = {}
        ref_len_dict = get_ref_len()
        self.window = 1
        for sample in self.data:
            sample_dict = {}
            for bkp in sample.bkps:
                from_ref = bkp.from_ref
                to_ref = bkp.to_ref

                from_b = int(bkp.from_bkp/self.window)
                to_b = int(bkp.to_bkp/self.window)

                if from_ref not in count_dict:
                    count_dict[from_ref] = {}
                if from_b not in count_dict[from_ref]:
                    count_dict[from_ref][from_b] = 0
                if from_ref + str(from_b) not in sample_dict:
                    count_dict[from_ref][from_b] += 1

                if to_ref not in count_dict:
                    count_dict[to_ref] = {}
                if to_b not in count_dict[to_ref]:
                    count_dict[to_ref][to_b] = 0
                if to_ref + str(to_b) not in sample_dict:
                    count_dict[to_ref][to_b] += 1

                sample_dict[to_ref + str(to_b)] = 1
                sample_dict[from_ref + str(from_b)] = 1
        origin_count_dict = count_dict.copy()

        density_list = []
        for genome in count_dict:
            num = 0
            for ele in count_dict[genome]:
                # if count_dict[genome][ele] > 5:
                num += 1   
            density =  (num/ref_len_dict[genome])
            if density > 0.8:
                print (genome, density)
            # density = round(density, 3)
            density_list.append(density)
            average_sample_count = sum(list(count_dict[genome].values()))/ref_len_dict[genome]
            count_dict[genome] = average_sample_count
        print (max(density_list))
        sorted_count_dict = self.show_sort_dict(count_dict)  
        df = pd.DataFrame(density_list, columns = ["proportion"])
        df.to_csv('/mnt/c/Users/swang66/Documents/genome_density_distribution.csv', sep = ",")
        # plt.figure(figsize=(15, 8))
        # plt.figure()
        # ax = sns.displot(df, x="Breakpoint_Num", bins=100, kde=True)
        # ax.set(xlabel='The number of breakpoint positions per bp')
        # plt.savefig("genome_density.pdf")

        """
        data = []
            # genome_index = 0
        print (len(sorted_count_dict))
        i = 0
        for genome_index in range(1000):
            genome = sorted_count_dict[genome_index][0]
            
            if int(ref_len_dict[genome]/self.window) > 450:
                i += 1
                for pos in range(450): #int(ref_len_dict[genome]/self.window)
                    if pos in origin_count_dict[genome]:
                        sam_num = origin_count_dict[genome][pos]
                    else:
                        sam_num = 0
                    data.append([pos*self.window, sam_num, genome])
                if i == 20:
                    break
        df = pd.DataFrame(data, columns = ["Position", "Sample_Number", "Genome"] )
        # print (df)
        # plt.figure()
        # sns.set_theme(style="whitegrid")
        # ax = sns.scatterplot(data=df, x="Position", y="Sample_Number", hue = "Genome", palette="ch:r=-.2,d=.3_r", sizes=(1, 8), linewidth=0)
        # ax.set(xlabel='Position')
        # plt.savefig("genome_sample.pdf")
        df.to_csv('genome_sample.csv', sep='\t')
        os.system("Rscript genome_density.R")
        """
        
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
            for bkp in sample.bkps:
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

    def differential_bkp(self, ann):
        level = 8
        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.data:
            if sample.disease == "CRC" :
                crc_num += 1
            if sample.disease == "control" :
                control_num += 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                family_edge = get_tag(bkp, 4)
                # if family_edge != "f__Lachnospiraceae&f__Lachnospiraceae":
                #     continue
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                if edge not in specific_HGT:
                    specific_HGT[edge] = [0, 0]
                if sample.disease ==  "CRC":
                    specific_HGT[edge][0] += 1
                if sample.disease ==  "control"  :
                    specific_HGT[edge][1] += 1    
        
        del self.data
        select_edges = {}
        for tag in specific_HGT:
            if specific_HGT[tag][0] + specific_HGT[tag][1] < 25:
                continue
            array = tag.split("&")
            species_1 = array[0]
            species_2 = array[1]
            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
            U1, p = mannwhitneyu(crc_array, control_array)
            select_edges[tag] = p 
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)
        with open("diff_bkps", "wb") as fp:
            pickle.dump(sort_select_edges[:30], fp)
        """
        gene_count = {}
        for i in range(30):
            print (sort_select_edges[i])
            edge = sort_select_edges[i][0]
            array = edge.split("&")
            genome_1 = array[0].split("|")[0]
            genome_1_pos = int(array[0].split("|")[1]) * window
            detail_1 = ann.given_point(genome_1, genome_1_pos)

            genome_2 = array[1].split("|")[0]
            genome_2_pos = int(array[1].split("|")[1]) * window
            detail_2 = ann.given_point(genome_2, genome_2_pos)
            
            # print (detail_1, detail_2)     
            for gene in detail_1 + detail_2:
                if gene not in gene_count:
                    gene_count[gene] = 1
                else:
                    gene_count[gene] += 1  
        print ("gene num", len(gene_count))
        sorted_gene_count = self.show_sort_dict(gene_count)
        data = []
        for i in range(len(sorted_gene_count)):
            data.append([sorted_gene_count[i][0], sorted_gene_count[i][1]])
        # df.to_csv('diff.gene_production.csv', sep='\t')
        """

    def most_HGT_bkp(self, ann): #the genome pair with most HGTs
        level = 7
        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.data:
            if sample.disease == "CRC" :
                crc_num += 1
            if sample.disease == "control"  :
                control_num += 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                # array = edge.split("&")
                # if "GUT_GENOME144544_1" not in array:
                #     continue
                if edge not in plot_genomes:
                    continue
                position = get_tag(bkp, 8)
                family_edge = get_tag(bkp, 4)
                if edge not in specific_HGT:
                    specific_HGT[edge] = {}
                if position not in specific_HGT[edge]:
                    specific_HGT[edge][position] = 1

        
        del self.data
        select_edges = {}
        for tag in specific_HGT:
            select_edges[tag] = len(specific_HGT[tag])

        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = True)
        new = {}
        for i in range(len(sort_select_edges)):
            edge = sort_select_edges[i][0]
            new[edge] = specific_HGT[edge]
        print (sort_select_edges)
        with open("diff_bkps_genomes", "wb") as fp:
            pickle.dump(new, fp)

    def compare_Lachnospiraceae_occurence(self): #the genome pair of Lachnospiraceae between CRC and control
        level = 6
        data = []
        for sample in self.data:
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                family_edge = get_tag(bkp, 4)
                if family_edge != "f__Lachnospiraceae&f__Lachnospiraceae":
                    continue
                if edge not in sample_dict:
                    sample_dict[edge] = 1
            data.append([len(sample_dict), sample.cohort, sample.disease, sample.bases])
        df = pd.DataFrame(data, columns = ["Events", "Cohort","Group", "Bases"])
        df.to_csv('/mnt/c/Users/swang66/Documents/lachnospiraceae_diff_HGTs.csv', sep=',')

    def differential_breakpoint_Lachnospiraceae(self, ann):
        level = 8
        ann.near = 5000  # the bkp near the point 
        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        best_num = 0
        for sample in self.data:
            min_split_num, p = get_split_reads_cutoff(sample.reads_num)

            if sample.disease == "CRC" :
                crc_num += 1
            if sample.disease == "control" :
                control_num += 1
            sample_dict = {}
            for bkp in sample.bkps:
                if bkp.cross_split_reads < min_split_num:
                    continue
                edge = get_tag(bkp, level)
                family_edge = get_tag(bkp, 4)
                if family_edge != "f__Lachnospiraceae&f__Lachnospiraceae":
                    continue
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                for bk in array:
                    if bk in sample_dict:
                        continue
                    else:
                        sample_dict[bk] = 1

                    if bk not in specific_HGT:
                        specific_HGT[bk] = [0, 0]
                    if sample.disease == "CRC":
                        specific_HGT[bk][0] += 1
                    if sample.disease == "control" :
                        specific_HGT[bk][1] += 1    
        
        del self.data
        ratio_data = []
        select_edges = {}
        for tag in specific_HGT:
            if specific_HGT[tag][0] + specific_HGT[tag][1] < 25:
                continue
            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
            # the HGT should occur in controls more frequently.
            crc_ratio = specific_HGT[tag][0]/crc_num
            control_ratio = specific_HGT[tag][1]/control_num

            U1, p = mannwhitneyu(crc_array, control_array)
            
            if p >= 0.05:
                continue
            ratio_data.append([crc_ratio, control_ratio, p])

            # if crc_ratio > control_ratio:
            #     continue
            select_edges[tag] = p 

        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

        gene_count = {}
        data = []
        for i in range(len(sort_select_edges)):
            if sort_select_edges[i][1] > 0.05:
                break
            if i < 5:
                print (i, sort_select_edges[i])
            bk = sort_select_edges[i][0]
            genome_1 = bk.split("|")[0]
            genome_1_pos = int(bk.split("|")[1]) * window
            genes_around = ann.given_point(genome_1, genome_1_pos)
   
            for gene in genes_around:
                if gene != "NA":
                    data.append(gene)
                # if gene in concern_genes:
                #     print (sort_select_edges[i])
        #         if gene not in gene_count:
        #             gene_count[gene] = 1
        #         else:
        #             gene_count[gene] += 1  
        # print ("gene num", len(gene_count))
        # sorted_gene_count = self.show_sort_dict(gene_count)
        
        # for i in range(len(sorted_gene_count)):
        #     print (sorted_gene_count[i])
        #     data.append(sorted_gene_count[i][0])
        # df = pd.DataFrame(ratio_data, columns = ["CRC", "control","Pvalue"])
        # df.to_csv('/mnt/c/Users/swang66/Documents/lachnospiraceae_diff_points.csv', sep=',')
        print ("diff genes number:", len(data))
        with open("diff_genes_Lach", "wb") as fp:
            pickle.dump(data, fp)

    def anno_compare(self, ann):
        level = 8
        specific_HGT = {} 
        for sample in self.data:
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                family_edge = get_tag(bkp, 4)
                # if family_edge != "f__Lachnospiraceae&f__Lachnospiraceae":
                #     continue
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                if edge in sample_dict:
                    continue
                else:
                    sample_dict[edge] = 1
                if edge not in specific_HGT:
                    specific_HGT[edge] = 0
                specific_HGT[edge] += 1
        print ("HGT num:", len(specific_HGT))
 
        gene_count = {}
        group_count = [{}, {}]
        crc_num = 0
        control_num = 0
        for sample in self.data:
            if sample.disease == "control" :
                group = 0
            if sample.disease == "CRC" :
                group = 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                family_edge = get_tag(bkp, 4)
                # if family_edge != "f__Lachnospiraceae&f__Lachnospiraceae":
                #     continue 
                # if specific_HGT[edge] < 200:
                #     continue 
        
                array = edge.split("&")
                genome_1 = array[0].split("|")[0]
                genome_1_pos = int(array[0].split("|")[1]) * window
                detail_1 = ann.given_point(genome_1, genome_1_pos)

                genome_2 = array[1].split("|")[0]
                genome_2_pos = int(array[1].split("|")[1]) * window
                detail_2 = ann.given_point(genome_2, genome_2_pos)

                for gene in detail_1 + detail_2:
                    if gene not in gene_count:
                        gene_count[gene] = 1
                    else:
                        gene_count[gene] += 1
                    if gene not in group_count[group]:
                        group_count[group][gene] = 1
                    else:
                        group_count[group][gene] += 1
        print ("gene num", len(gene_count))
        sorted_gene_count = self.show_sort_dict(gene_count)
        self.show_sort_dict(group_count[0])
        self.show_sort_dict(group_count[1])
        data = []
        for i in range(100):
            data.append([sorted_gene_count[i][0], sorted_gene_count[i][1]])
        df.to_csv('all.gene_production.csv', sep='\t')
        # print (df)

    def draw_bkp_plots(self, ann):
        level = 7
        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        for sample in self.data:
            if sample.disease == "control" :
                crc_num += 1
            if sample.disease == "CRC" :
                control_num += 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                locus_tag = get_tag(bkp, 8)
                if edge not in specific_HGT:
                    specific_HGT[edge] = {}
                if locus_tag not in specific_HGT[edge]:
                    specific_HGT[edge][locus_tag] = 0
                specific_HGT[edge][locus_tag] +=1

        print (len(specific_HGT))
        for edge in specific_HGT:
            # print (len(specific_HGT[edge]))
            for locus_tag in specific_HGT[edge].copy():
                if specific_HGT[edge][locus_tag] < 50:
                    del specific_HGT[edge][locus_tag]
            if len(specific_HGT[edge]) > 1:
                print (edge, specific_HGT[edge])


        # sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

        # gene_count = {}
        # for i in range(30):
        #     print (sort_select_edges[i])

        # df.to_csv('diff.gene_production.csv', sep='\t')

    def HGT_count_with_samples(self): #the genome pair with most HGTs      
        random.seed(0)
        random.shuffle(self.data)
        dat = []
        group_dat = []
        for level in range(1, 8):
            print (level)
            sample_index = 1
            specific_HGT = {} 

            group_index = [1, 1]
            group_HGT = [{},{}]

            for sample in self.data:
                if sample.disease == "CRC":
                    i = 0
                else:
                    i = 1
                for bkp in sample.bkps:
                    HGT = get_tag(bkp, level)
                    if HGT not in specific_HGT:
                        specific_HGT[HGT] = 1
                    if HGT not in group_HGT[i]:
                        group_HGT[i][HGT] = 1

                dat.append([sample_index, len(specific_HGT), level, sample.disease])
                group_dat.append([group_index[i], len(group_HGT[i]), level, sample.disease])
                sample_index += 1
                group_index[i] += 1
                # print (sample_index, len(specific_HGT))


        df = pd.DataFrame(dat, columns = ["sample", "HGT_count", "level", "group"])
        df.to_csv('/mnt/c/Users/swang66/Documents/HGT_count_sample.csv', sep='\t')
        df = pd.DataFrame(group_dat, columns = ["sample", "HGT_count", "level", "group"])
        df.to_csv('/mnt/c/Users/swang66/Documents/HGT_count_sample_group.csv', sep='\t')

    def HGT_individual_with_samples(self): #the genome pair with most HGTs      
        random.seed(0)
        random.shuffle(self.data)
        dat = []
        group_dat = []
        for level in range(1, 8):
            print (level)
            sample_index = 1
            specific_HGT = {} 

            group_index = [1, 1]
            group_HGT = [{},{}]

            for sample in self.data:
                if sample.disease == "CRC":
                    i = 0
                else:
                    i = 1
                sample_dict = {}
                for bkp in sample.bkps:
                    HGT = get_tag(bkp, level)
                    if HGT not in sample_dict:
                        sample_dict[HGT] = 1

                dat.append([len(sample_dict), level, sample.disease])
                group_dat.append([len(sample_dict), level, sample.disease])
                sample_index += 1
                group_index[i] += 1
                # print (sample_index, len(specific_HGT))


        df = pd.DataFrame(dat, columns = ["genome_pair", "level", "group"])
        df.to_csv('/mnt/c/Users/swang66/Documents/HGT_individual_sample.csv', sep='\t')
        df = pd.DataFrame(group_dat, columns = ["genome_pair", "level", "group"])
        df.to_csv('/mnt/c/Users/swang66/Documents/HGT_individual_sample_group.csv', sep='\t')

    def select_circos_genome(self):
        select_genome = {}
        level = 7
        specific_HGT = {} 
        crc_num = 0
        control_num = 0
        for sample in self.data:
            if sample.disease == "control" :
                crc_num += 1
            if sample.disease == "CRC" :
                control_num += 1
            sample_dict = {}
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                family_edge = get_tag(bkp, 4)
                # if family_edge != "f__Lachnospiraceae&f__Lachnospiraceae":
                #     continue

                if edge not in specific_HGT:
                    specific_HGT[edge] = [0, 0]
                if sample.disease == "control" :
                    specific_HGT[edge][0] += 1
                if sample.disease == "CRC" :
                    specific_HGT[edge][1] += 1    
        
        select_edges = {}
        for tag in specific_HGT:
            if specific_HGT[tag][0] + specific_HGT[tag][1] < 25:
                continue
            crc_array = [1] * specific_HGT[tag][0] + [0] * (crc_num - specific_HGT[tag][0])
            control_array = [1] * specific_HGT[tag][1] + [0] * (control_num - specific_HGT[tag][1])
            U1, p = mannwhitneyu(crc_array, control_array)
            select_edges[tag] = p 
        sort_select_edges = sorted(select_edges.items(), key=lambda item: item[1], reverse = False)

    def prepare_circos(self, ann):
        count_dict = {}
        ref_len_dict = get_ref_len()
        point_frequency = {}
        group_dict = [{}, {}]
        select_genome = {'GUT_GENOME143505_1':1, "GUT_GENOME144544_1":1, "GUT_GENOME147149_1":1,'GUT_GENOME143712_1':1,'GUT_GENOME143131_1':1\
        ,'GUT_GENOME096067_1':1,'GUT_GENOME096063_2':1,'GUT_GENOME147164_1':1,'GUT_GENOME095938_1':1,'GUT_GENOME147678_1':1}
        select_genome_bkps = {}
        links_dict = {}

        self.window = 100000
        for sample in self.data:
            if sample.disease == "CRC":
                group_index = 0
            else:
                group_index = 1
            
            sample_dict = {}
            for bkp in sample.bkps:
                from_ref = bkp.from_ref
                to_ref = bkp.to_ref

                if from_ref in select_genome:
                    if from_ref not in select_genome_bkps:
                        select_genome_bkps[from_ref] = {}
                    select_genome_bkps[from_ref][round(bkp.from_bkp/100)] = 1
                if to_ref in select_genome:
                    if to_ref not in select_genome_bkps:
                        select_genome_bkps[to_ref] = {}
                    select_genome_bkps[to_ref][round(bkp.to_bkp/100)] = 1
                if from_ref in select_genome and to_ref in select_genome:
                    link = from_ref + "|" + str(round(bkp.from_bkp/10000)) + "|" + to_ref + "|" + str(round(bkp.to_bkp/10000))
                    if link not in links_dict:
                        links_dict[link] = 0 
                    links_dict[link] += 1

                from_b = round(bkp.from_bkp/self.window)
                to_b = round(bkp.to_bkp/self.window)

                from_point = from_ref+ "|" + str(from_b)
                to_point = to_ref+ "|" + str(to_b)

                if from_ref not in count_dict:
                    count_dict[from_ref] = {}
                if from_point not in count_dict[from_ref]:
                    count_dict[from_ref][from_point] = 0
                if from_point not in sample_dict:
                    count_dict[from_ref][from_point] += 1

                if to_ref not in count_dict:
                    count_dict[to_ref] = {}
                if to_point not in count_dict[to_ref]:
                    count_dict[to_ref][to_point] = 0
                if to_point not in sample_dict:
                    count_dict[to_ref][to_point] += 1

                for my_point in [from_point, to_point]:
                    if my_point not in sample_dict:
                        if my_point not in group_dict[group_index]:
                            group_dict[group_index][my_point] = 0
                        group_dict[group_index][my_point] += 1

                sample_dict[to_point] = 1
                sample_dict[from_point] = 1

                for my_point in [from_point, to_point]:
                    if my_point not in point_frequency:
                        point_frequency[my_point] = 0
                    point_frequency[my_point] += 1

       
        point_dict = {}
        for genome in count_dict:
            num = 0
            for ele in count_dict[genome]:
                if count_dict[genome][ele] > 100:
                    num += 1   
            point_dict[genome] = num
        data = []
        sorted_point_dict = self.show_sort_dict(point_dict)  
        for_heatmap = []
        for i in range(len(select_genome)):
            genome = sorted_point_dict[i][0]
            for point in count_dict[genome]:
                locus = int(point.split("|")[1]) * self.window
                frequency = point_frequency[point]
                crc_num = 0
                if point in group_dict[0]:
                    crc_num = group_dict[0][point]
                control_num = 0
                if point in group_dict[1]:
                    control_num = group_dict[1][point]

                data.append([genome, locus, frequency, crc_num, control_num])
                for_heatmap.append([genome, locus, locus+self.window, crc_num/373, control_num/395])

        my_bkps = []
        for genome in select_genome_bkps:
            for key in select_genome_bkps[genome]:
                my_bkps.append([genome, int(key)*100])

        bed1 = []
        bed2 = []
        for link in links_dict:
            array = link.split("|")
            support = links_dict[link]
            if support > 0:
                bed1.append([array[0], int(array[1])*10000, int(array[1])*10000+10000, support])
                bed2.append([array[2], int(array[3])*10000, int(array[3])*10000+10000, support])

        df = pd.DataFrame(data, columns = ["Genome", "position", "freqency", "CRC", "control"] )
        df.to_csv('/mnt/c/Users/swang66/Documents/circos_data.csv', sep='\t')
        df = pd.DataFrame(for_heatmap, columns = ["Genome", "start", "end", "CRC", "control"] )
        df.to_csv('/mnt/c/Users/swang66/Documents/for_heatmap_circos_data.csv', sep='\t')
        df = pd.DataFrame(my_bkps, columns = ["Genome", "bkp"] )
        df.to_csv('/mnt/c/Users/swang66/Documents/for_hist_circos_data.csv', sep='\t')
        df = pd.DataFrame(bed1, columns = ["Genome", "start", "end", "support"] )
        df.to_csv('/mnt/c/Users/swang66/Documents/for_bed1_circos_data.csv', sep='\t')
        df = pd.DataFrame(bed2, columns = ["Genome", "start", "end", "support"] )
        df.to_csv('/mnt/c/Users/swang66/Documents/for_bed2_circos_data.csv', sep='\t')

    def compare_depth(self):

        record = [[], []]

        data = []
        for sample in self.data:

            if sample.disease == "control" :
                index = 0
            else:
                index = 1
            record[index].append(sample.bases)
            data.append([sample.bases, sample.disease, sample.cohort])
        U1, p = mannwhitneyu(record[0], record[1])         
        print (p, sep = "\t")
        df = pd.DataFrame(data, columns = ["Bases", "Group", "Cohort"])
        df.to_csv('/mnt/c/Users/swang66/Documents/depth_comparison.csv', sep=',')
        # print (df) 

    def __del__(self):
        # del self.data
        del self.disease_sample_num_cohort
        print ("object, deleted")

class Function():
    def __init__(self):
        self.sample_genes = []

    def find_genes(self, all_data, ann):
        level = 8
        sample_index = 0
        # ann.near = 5000  ###
        for sample in all_data:
            if sample.disease == "control" :
                group = 0
            if sample.disease == "CRC" :
                group = 1
            genes = {}
            genes["group"] = group
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                for bk in array:
                    genome_1 = bk.split("|")[0]
                    genome_1_pos = int(bk.split("|")[1]) * window
                    genes_around = ann.given_point(genome_1, genome_1_pos)
        
                    for gene in genes_around:
                        if gene == "NA":
                            continue
                        if gene not in genes:
                            genes[gene] = 1
            self.sample_genes.append(genes)
            sample_index += 1
            print (sample_index, len(genes))
        with open("sample_genes", "wb") as fp:
            pickle.dump(self.sample_genes, fp)

    def screen_COG(self, ann):
        annotation_dict = {}
        for cog in COG_annotation:
            name = cog.split(":")[0]
            annotation_dict[name] = cog
        # print (annotation_dict)
        with open("sample_genes", "rb") as fp:
            self.sample_genes = pickle.load(fp)
        data = []
        j = 0
        shuffle(self.sample_genes)
        for sample in self.sample_genes:   
            dict = {}
            for key in COG_categories:
                dict[key] = 0
            # if sample["group"] == 0:
            #     continue
            for ID in sample:
                if ID == "group":
                    continue
                if ID == "NA":
                    continue
                info = ann.gene_function[ID].strip()
                array = info.split(";")
                anno_dict = {}
                for ele in array:
                    arr = ele.split("=")
                    anno_dict[arr[0]] = arr[1]
                if "COG" not in anno_dict:
                    continue
                category = anno_dict["COG"]
                for i in range(len(category)):
                    dict[category[i]] += 1
            # dict["S"] = 0
            del dict["S"]
            category_count = np.array(list(dict.values()))
            COG_keys = list(dict.keys())
            if sum(category_count) > 0:
                frequency = category_count/sum(category_count)
            for i in range(len(frequency)):
                data.append([j, frequency[i], COG_keys[i], sample["group"], annotation_dict[COG_keys[i]] ])
            # data.append(frequency)
            # print (frequency)
            j += 1
        df = pd.DataFrame(data, columns = ["Sample", "Frequency", "COG", "group", "Function"])
        df.to_csv('/mnt/c/Users/swang66/Documents/COG_categories.csv', sep=',')
        # os.system("Rscript COG_heatmap.R")

    def screen_Product(self, ann):
        with open("sample_genes", "rb") as fp:
            self.sample_genes = pickle.load(fp)
        data = []
        j = 0
        shuffle(self.sample_genes)
        product_count = {}
        for sample in self.sample_genes:   
            dict = {}
            # if sample["group"] == 0:
            #     continue
            sample_dict = {}
            for ID in sample:
                if ID == "group" or ID == "NA":
                    continue
                info = ann.gene_function[ID].strip()
                array = info.split(";")
                anno_dict = {}
                for ele in array:
                    arr = ele.split("=")
                    anno_dict[arr[0]] = arr[1]
                if "product" not in anno_dict:
                    continue
                product = anno_dict["product"]
                if product not in sample_dict:
                    if product not in product_count:
                        product_count[product] = 0
                    product_count[product] += 1
                # sample_dict[product] = 1
            j += 1
        sorted_count = sorted(product_count.items(), key=lambda item: item[1], reverse = True)

        for i in range(2, 30):
            print (sorted_count[i])
            data.append(list(sorted_count[i]))
        df = pd.DataFrame(data, columns = ["Product", "Count"])
        df.to_csv('product_count.csv', sep='\t')      

    def screen_Product_Lachnospiraceae(self, ann):
        with open("diff_genes_Lach", "rb") as fp:
            my_genes = pickle.load(fp)
        print (len(my_genes), "diff_genes_Lach")
        data = []
        j = 0
        product_count = {}
        for ID in my_genes:   
            if ID == "group" or ID == "NA":
                continue
            info = ann.gene_function[ID].strip()
            array = info.split(";")
            anno_dict = {}
            for ele in array:
                arr = ele.split("=")
                anno_dict[arr[0]] = arr[1]
            if "product" not in anno_dict:
                continue
            product = anno_dict["product"]

            if product not in product_count:
                product_count[product] = 0
            product_count[product] += 1
            # sample_dict[product] = 1
            j += 1
            # if j > 1000:
            #     break
        sorted_count = sorted(product_count.items(), key=lambda item: item[1], reverse = True)

        for i in range(len(sorted_count)):
            # print (sorted_count[i])
            data.append(list(sorted_count[i]))
        print (j)
        df = pd.DataFrame(data, columns = ["Product", "Count"])
        df.to_csv('/mnt/c/Users/swang66/Documents/product_count_Lachnospiraceae.csv', sep='\t')   

    def screen_class(self, ann):
        with open("sample_genes", "rb") as fp:
        # with open("diff_genes", "rb") as fp:
            self.sample_genes = pickle.load(fp)
        # print ("self.sample_genes", len(self.sample_genes ))
        # print ("ann.gene_classification", ann.gene_classification)
        # count = [[], []]
        shuffle(self.sample_genes)
        group_count = {}
        for sample in self.sample_genes:  
            
            # if sample["group"] == 0:
            #     continue 
            for ID in sample:
        # if True:
        #     for ID in self.sample_genes:
                if ID in ann.gene_classification:
                    classify_group = ann.gene_classification[ID]
                # else:
                #     classify_group = "others"

                    if classify_group not in group_count:
                        group_count[classify_group] = 0
                    group_count[classify_group] += 1

            # x =  np.array(list(group_count.values()))
            # for key in group_count:
            #     if key == "Transposon":
            #         print (sample["group"], group_count[key]/sum(x))
            #         count[sample["group"]].append(group_count[key]/sum(x))

        # U1, p = mannwhitneyu(count[0], count[1])
        # print (p, np.mean(count[0]), np.mean(count[1]))

        # print (group_count)
        sum_count = sum(list(group_count.values()))
        plt.pie(list(group_count.values()),
                labels=list(group_count.keys())
            )  
        plt.savefig("gene_pie.pdf")  

        # print (list(group_count.keys()))
        # print (x/sum(x)) 
        sorted_count = sorted(group_count.items(), key=lambda item: item[1], reverse = True)

        for i in range(len(sorted_count)):
            print (sorted_count[i], round(sorted_count[i][1]/sum_count, 3) )
        #     data.append(list(sorted_count[i]))
        # df = pd.DataFrame(data, columns = ["Product", "Count"])
        # df.to_csv('group_count.csv', sep='\t')  

    def transfer_genes(self, ann):
        with open("transferred_record", "rb") as fp:
            transferred_record = pickle.load(fp)
        sample_genes = {}
        dict = {}
        for key in COG_categories:
            dict[key] = 0
        for sample in transferred_record:
            sample_genes[sample] = []
            trans_genes = []
            for event in transferred_record[sample]:
                insert_geno = event[2]
                insert_pos = event[3]
                trans_geno = event[0]
                trans_interval = event[1]
                genes_around = ann.given_seg(trans_geno, trans_interval)
                for gene in genes_around:
                    if gene == "NA":
                        continue
                    trans_genes.append(gene)

            for ID in trans_genes:
                info = ann.gene_function[ID].strip()
                array = info.split(";")
                anno_dict = {}
                for ele in array:
                    arr = ele.split("=")
                    anno_dict[arr[0]] = arr[1]
                if "COG" not in anno_dict:
                    continue
                category = anno_dict["COG"]
                for i in range(len(category)):
                    dict[category[i]] += 1
        print (dict)
        sorted_count = sorted(dict.items(), key=lambda item: item[1], reverse = True)

        for i in range(len(sorted_count)):
            print (sorted_count[i])

            # sample_genes[sample] = trans_genes
        #     print (sample, trans_genes)
        # with open("transferred_genes", "wb") as fp:
        #     pickle.dump(sample_genes, fp)

    def diff_COG(self, ann):
        with open("diff_genes", "rb") as fp:
            diff_genes = pickle.load(fp)
        dict = {}
        for key in COG_categories:
            dict[key] = 0
        for ID in diff_genes:
            if ID == "NA":
                continue
            info = ann.gene_function[ID].strip()
            array = info.split(";")
            anno_dict = {}
            for ele in array:
                arr = ele.split("=")
                anno_dict[arr[0]] = arr[1]
            if "COG" not in anno_dict:
                continue
            category = anno_dict["COG"]
            for i in range(len(category)):
                if category[i] == "S":
                    continue
                dict[category[i]] += 1
        print (dict)
        sorted_count = sorted(dict.items(), key=lambda item: item[1], reverse = True)

        for i in range(len(sorted_count)):
            print (sorted_count[i])

    def diff_gene_pathway(self, ann):
        # beta_Lactam = ["K02171", "K02172", "K02545", "K02547", "K15580", "K15581", "K15582", "K15583", "K17838", "K18149", "K19209", "K19213", "K21276", "K22335", "K22352"] 

        with open("diff_genes_Lach", "rb") as fp:
            diff_genes = pickle.load(fp)
        f = open("diff_genes_Lach_ko.list", 'w')
        i = 0
        for ID in diff_genes:
            if ID == "NA":
                continue
            info = ann.gene_function[ID].strip()
            array = info.split(";")
            anno_dict = {}
            for ele in array:
                arr = ele.split("=")
                anno_dict[arr[0]] = arr[1]
            if "Name" not in anno_dict:
                continue
            if "KEGG" in anno_dict:
                for k in anno_dict["KEGG"].split(','):
                    ko = k[3:]
                    # if ko in beta_Lactam:
                    #     print (ID, anno_dict["Name"])
                    print (i+1, ko, file = f)
                    # print (i+1, ko)
            i += 1
            # if i > 1000:
            #     break
            name = anno_dict["Name"]  
            # print (name) 
        print (i)       
        f.close()

class Family(): # study family f__Lachnospiraceae

    def __init__(self):
        with open('selected_diff_edges.pkl', 'rb') as f:
            self.select_edges = pickle.load(f)
        self.taxonomy = Taxonomy()
        self.cho_top_hgt_num = 30
        self.genus2family = {}

    def get_family(self, lineage):
        array = lineage.split(";")
        return array[-3]
    
    def get_marker_tree(self):
        marker_genus = {}
        for edge in list(self.select_edges.keys())[:self.cho_top_hgt_num]:
            array = edge.split("&")
            if array[0] not in marker_genus:
                marker_genus[array[0]] = 1
            if array[1] not in marker_genus:
                marker_genus[array[1]] = 1
        print ("genus number:", len(marker_genus), len(marker_genus)/self.cho_top_hgt_num)
        
        for lineage in self.taxonomy.taxonomy_dict.values():
            genus = lineage.split(";")[5]
            if genus in marker_genus:
                marker_genus[genus] = lineage
        nodes_list = []
        edges_list = []
        lac_count = 0 # f__Lachnospiraceae count
        for genus in marker_genus:
            lineage = marker_genus[genus]
            array = lineage.split(";")
            # print (array)
            for i in range(6):
                nodes_list.append(array[i])
                if i <5:
                    edges_list.append([array[i], array[i+1]]) 
            if array[-3] == "f__Lachnospiraceae":
                lac_count += 1
        print ("f__Lachnospiraceae count", lac_count, lac_count/self.cho_top_hgt_num)

        self.plot_tree(nodes_list, edges_list)
        self.compare_ratio(marker_genus)

    def compare_ratio(self, marker_genus):
        data = []
        lac_pair = 0
        for edge in list(self.select_edges.keys())[:self.cho_top_hgt_num]:
            array = edge.split("&")
            if self.get_family(marker_genus[array[0]]) == "f__Lachnospiraceae" and self.get_family(marker_genus[array[1]]) == "f__Lachnospiraceae":
                lac_pair += 1
                data.append([edge, self.select_edges[edge][0], "CRC" ])
                data.append([edge, self.select_edges[edge][1], "control"])
            # else:
            #     print (edge, select_edges[edge])
        df = pd.DataFrame(data, columns = ["HGT", "Proportion", "Group"])
        df.to_csv('for_Lachnospiraceae.csv', sep='\t')
        os.system("Rscript Lachnospiraceae.R")
        print ("f__Lachnospiraceae pair count", lac_pair, lac_pair/self.cho_top_hgt_num)

    def genus_family(self):
        for lineage in self.taxonomy.taxonomy_dict.values():
            genus = lineage.split(";")[5]  
            family = lineage.split(";")[4]  
            self.genus2family[genus] = family
    
    def plot_all_diff_hgts(self):
        self.genus_family()
        ratio_data = []
        family_num = 0
        for edge in self.select_edges:
            array = edge.split("&")
            if self.genus2family[array[0]] == "f__Lachnospiraceae" and self.genus2family[array[1]] == "f__Lachnospiraceae":
                ratio_data.append(self.select_edges[edge])
                family_num += 1
        print ("total family pair", family_num, len(self.select_edges), family_num/len(self.select_edges) )
        df = pd.DataFrame(ratio_data, columns = ["CRC", "control","Pvalue"])
        df.to_csv('/mnt/c/Users/swang66/Documents/lachnospiraceae_diff_points.csv', sep=',')

    def plot_tree(self, nodes_list, edges_list):
        new_file = open("/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files/hgt_tree.nwk", 'w')
        DG = nx.DiGraph()
        DG.add_nodes_from(nodes_list)
        DG.add_edges_from(edges_list)   
        # print (marker_genus, nx.is_tree(DG))
        color_map  = []
        for node in DG:
            if node[0] == "g":
                color_map.append("#E14D2A")
            if node[0] == "f":
                color_map.append("#FD841F")
            if node[0] == "o":
                color_map.append("#3E6D9C")
            if node[0] == "c":
                color_map.append("#001253")   
            if node[0] == "p":
                color_map.append("#66A61E")    
            if node[0] == "d":  
                color_map.append("#E6AB02") 
            # print (node)

        pos = nx.nx_agraph.graphviz_layout(DG, prog="twopi", args="")
        plt.figure(figsize=(10, 10))
        nx.draw(DG, pos, node_size=2500, alpha=0.7, node_color=color_map, \
        with_labels=True, font_size = 5, font_color="black", width=1, edge_color="grey")
        plt.axis("equal")
        plt.savefig("HGT_tree.pdf")
        nwk_tree = tree_to_newick(DG, 'd__Bacteria')
        print (nwk_tree)
        print (nwk_tree, file = new_file)
        new_file.close()

    def compare_cccur_in_each_cohort(self):
        all_data, sample_num = get_data()
        cohort_hgt = {}
        cohort_dict = {}
        for sample in all_data:
            sample_dict = {}
            min_split_num, p = get_split_reads_cutoff(sample.reads_num)
            # print (m, p)
            for bkp in sample.bkps:
                edge = get_tag(bkp, level)
                array = edge.split("&")
                # print (edge, array)
                if len(array[0].strip()) <= 3 or len(array[1].strip()) <= 3:
                    continue
                family_edge = get_tag(bkp, 4)
                array = family_edge.split("&")
                if not (array[0] == "f__Lachnospiraceae" and array[0] == "f__Lachnospiraceae"):
                    continue
                if edge not in sample_dict:
                    sample_dict[edge] = 0
                    continue
                sample_dict[edge] = bkp.cross_split_reads # count genus pair one time.

            for edge in sample_dict:  
                if sample_dict[edge] < min_split_num:
                    continue
                if edge not in cohort_hgt:
                    cohort_hgt[edge] = {}
                if sample.cohort not in cohort_hgt[edge]:
                    cohort_hgt[edge][sample.cohort] = [0, 0]
                if sample.disease == "CRC" :
                    cohort_hgt[edge][sample.cohort][0] += 1
                else:
                    cohort_hgt[edge][sample.cohort][1] += 1
                

            if sample.cohort not in cohort_dict:
                cohort_dict[sample.cohort] = [0, 0]
            if sample.disease == "CRC" :
                cohort_dict[sample.cohort][0] += 1
            else:
                cohort_dict[sample.cohort][1] += 1
        # print (cohort_dict)
        # print (cohort_hgt)
        data = []
        for cohort in cohort_dict:
            crc_ratio_list = []
            control_ratio_list = []
            # print (cohort)
            for hgt in cohort_hgt:
                # print (cohort_hgt[hgt])
                sample_num = 0
                for each in cohort_hgt[hgt]:
                    sample_num += sum(cohort_hgt[hgt][each])
                if sample_num < 25:
                    continue
                if cohort in cohort_hgt[hgt]:
                    crc_ratio = cohort_hgt[hgt][cohort][0] / cohort_dict[cohort][0]
                    crc_ratio_list.append(crc_ratio)

                    control_ratio = cohort_hgt[hgt][cohort][1] / cohort_dict[cohort][1]
                    control_ratio_list.append(control_ratio)

                    data.append([control_ratio, "control", cohort])
                    data.append([crc_ratio, "CRC", cohort])
                    
            # print (crc_ratio_list, control_ratio_list)
            U1, p = mannwhitneyu(crc_ratio_list, control_ratio_list, method="auto")
            print (cohort, p)
        df = pd.DataFrame(data, columns = ["Ratio", "Group","cohort"])
        df.to_csv('/mnt/c/Users/swang66/Documents/lachnospiraceae_each_cohort_ratio.csv', sep=',')            
        
class Annotation():

    def __init__(self, gff):
        self.gff = gff
        self.near = 100
        self.gene_annotation = {}
        self.gene_function = {}
        self.gene_classification = {}

    def read_gff(self):
        f = open(self.gff)
        for line in f:
            array = line.split("\t")
            genome = array[0]
            g_type = array[2]
            detail = array[8]
            start = int(array[3])
            end = int(array[4])
            if genome not in self.gene_annotation:
                self.gene_annotation[genome] = {}
                self.gene_annotation[genome]["intervals"] = []
            self.gene_annotation[genome]["intervals"].append([start, end])
            name = self.get_name(detail)
            self.gene_annotation[genome][str(start)+ "_"+str(end)] = name 
            self.gene_function[name] = detail
        f.close() 
        with open("gene_annotation", "wb") as fp:
            pickle.dump(self.gene_annotation, fp)
        with open("gene_function", "wb") as fp:
            pickle.dump(self.gene_function, fp)

    def read_dict(self):
        with open("gene_annotation", "rb") as fp:
            self.gene_annotation = pickle.load(fp)
        with open("gene_function", "rb") as fp:
            self.gene_function = pickle.load(fp)

    def given_point(self, genome, locus):
        if genome not in self.gene_annotation:
            return ["NA"]
        intervals = self.gene_annotation[genome]["intervals"]
        genes_around = []
        for inter in intervals:
            # print (inter)
            if locus >= inter[0] - self.near and locus <= inter[1] + self.near:
                gene_ID = self.gene_annotation[genome][str(inter[0])+ "_"+str(inter[1])]
                genes_around.append(gene_ID)
        if len(genes_around) > 0:
            return genes_around
        else:
            return ["NA"]

    def given_seg(self, genome, gene_interval):
        if genome not in self.gene_annotation:
            return ["NA"]
        intervals = self.gene_annotation[genome]["intervals"]
        genes_around = []
        for inter in intervals:
            # print (inter)
            if gene_interval[0] >= inter[0] - 5 and gene_interval[1] <= inter[1] + 5:
                gene_ID = self.gene_annotation[genome][str(inter[0])+ "_"+str(inter[1])]
                genes_around.append(gene_ID)
        if len(genes_around) > 0:
            return genes_around
        else:
            return ["NA"]

    def split_anno(self, info):
        match = re.search("product=(.*?);", info)
        # return info.split(";")[7]
        if match:
            return match.group(1)
        elif re.search("product=(.*?)$", info):
            return re.search("product=(.*?)$", info).group(1)
        else:
            return info

    def get_name(self, info):
        match = re.search("ID=(.*?);", info)
        # return info.split(";")[7]
        if match:
            return match.group(1)
        elif re.search("ID=(.*?)$", info):
            return re.search("ID=(.*?)$", info).group(1)
        # else:
        #     match = re.search("ID=(.*?);", info)
        #     return match.group(1)

    def classify(self):
        Transposon = "transpos*; insertion; resolv*; Tra[A-Z]; Tra[0-9]; IS[0-9]; conjugate transposon"
        Plasmid = "resolv*; relax*; conjug*; trb; mob*; plasmid; “type IV”; toxin; “chromosome partitioning”; “chromosome segregation”"
        Phage ="capsid; phage; tail; head; “tape measure”; antitermination"
        Other_HGT_machinery=" integrase; excision*; exonuclease; recomb; toxin; CRISPR; restrict*; resolv*; topoisomerase; reverse transcrip"
        Carbohydrate_active =  "Genes present in the CAZY database; glycosyltransferase; “glycoside hydrolase; xylan; monooxygenase; rhamnos*; cellulose; sialidase; *ose; acetylglucosaminidase; cellobiose; galact*; fructose; aldose; starch; mannose; mannan*; glucan; lyase; glycosyltransferase; glycosidase; pectin; SusD; SusC; fructokinase; galacto*; arabino*"
        antibiotic_resistance =  "Genes present in the ARDB; multidrug; “azole resistance”; antibiotic resistance”; TetR; “tetracycline resistance”; VanZ; betalactam*; beta-lactam; antimicrob*; lantibio*"
        
        pattern_dict = {}
        gene_classification = {}
        pattern_dict["Transposon"] = get(Transposon)
        pattern_dict["Plasmid"] = get(Plasmid)
        pattern_dict["Phage"] = get(Phage)
        pattern_dict["Other"] = get(Other_HGT_machinery)
        pattern_dict["CAZYmes"] = get(Carbohydrate_active)
        pattern_dict["Antibiotic resistance"] = get(antibiotic_resistance)

        # print (pat)
        f = open(self.gff)
        for line in f:
            if re.search("ncRNA", line):
                continue
            if re.search("product=(.*?);", line):
                product = re.search("product=(.*?);", line).group(1)
            elif re.search("product=(.*?)$", line):
                product = re.search("product=(.*?)$", line).group(1)
            else:
                continue
            gene_ID = self.get_name(line)
            for cla in pattern_dict:
                # print (cla, pattern_dict[cla])
                pattern = re.compile(pattern_dict[cla]) #, re.I
                if pattern.search(product):
                    # print (gene_ID, cla)
                    gene_classification[gene_ID] = cla
            # break
        f.close()
        with open("gene_classification", "wb") as fp:
            pickle.dump(gene_classification, fp)
        print ("gene classsification is done.")

    def read_classification(self):
        with open("gene_classification", "rb") as fp:
            self.gene_classification = pickle.load(fp)            

def get(x):
    x = x.split(";")
    # print (x)
    for i in range(len(x)):
        x[i] = x[i].strip()
        if x[i][0] == "“":
            x[i] = x[i][1:]
        if x[i][0] == "*":
            x[i] = x[i][1:]
        if x[i][-1] == "”":
            x[i] = x[i][:-1]
        if x[i][-1] == "*":
            x[i] = x[i][:-1]
    pat = ''
    for i in range(len(x)):
        pat += x[i] + "|"
    if pat[-1] == "|":
        pat = pat[:-1]
    return pat

def get_ref_len():
    ref_len_dict = {}
    fai = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.fai"
    for line in open(fai):
        array = line.split()
        ref_len_dict[array[0]] = int(array[1])
    return ref_len_dict

def examine_gene(ann):
    # convern_fea = ["Name", "product", "Pfam", "eggNOG", "KEGG"]
    convern_fea = ["Name", "COG"]
    with open("diff_genes", "rb") as fp:
        data = pickle.load(fp)  
    for ID in data:
        if ID == "NA":
            continue
        info = ann.gene_function[ID].strip()
        array = info.split(";")
        dict = {}
        for ele in array:
            arr = ele.split("=")
            dict[arr[0]] = arr[1]
        if "Name" not in dict:
            continue
        if dict["product"] != "hypothetical protein":
            # print (dict.values())
            for feature in convern_fea:
                if feature in dict:
                    print (dict[feature], end = ";\t")
                else:
                    print ("NA", end = ";\t")
            print ('')
 
def count_figures(ana, ann):
    # ana.HGT_count_with_samples()
    # ana.HGT_individual_with_samples()
    # ana.taxonomy_circos()
    # ana.taxonomy_circos_CRC()
    # ana.taxonomy_circos_cohort()
    # ana.bkp_num_histogram_cohort()
    ana.output_samples()
    # ana.taxonomy_count()
    # ana.taxonomy_count_cohort()
    # ana.taxonomy_count_group()
    # ana.prepare_circos(ann)

def main(): 
    gff = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.gff"
    ann = Annotation(gff)
    # ann.classify()
    ann.read_gff()
    # # ann.read_dict()
    ann.read_classification()
    # print ("annotation readed.")
    
    # breakpoint_genes = "sample_genes"
    # trans_genes = "transferred_genes"
    # print ("............")
    # draw_HGT_genome_pair(ann)
    # draw_HGT(ann)
    all_data, sample_num = get_data()
    print ("data is loaded.")
    net = Network(all_data)
    # net.compare_network()   
    # net.compare_network_test()
    net.infer_sale()

    # ana = Analyze(3, all_data)
    # ana.cal_average_bkp_num()
    # ana.inter_taxa_circos()
    # ana.differential_breakpoint_Lachnospiraceae(ann)
    # count_figures(ana, ann)  # For figure 1 in manuscript
    # ana.HGT_individual_with_samples()
    # ana.compare_depth()
    # ana.HGT_individual_with_samples()
    # ana.prepare_circos(ann)
    # ana.count_genome_density()
    # # # ana.taxonomy_barplot()
    # # # ana.taxonomy_count()
    # ana.taxonomy_circos()
    # ana.draw_bkp_plots(ann)
    # ana.differential_breakpoint(ann)
    # ana.compare_Lachnospiraceae_occurence()
    # ana.most_HGT_bkp(ann)
    # examine_gene(ann)

    # fun = Function()
    # print ("****")
    # fun.screen_Product_Lachnospiraceae(ann)
    # fun.screen_COG(ann)
    # fun.find_genes(all_data, ann)
    # fun.screen_class(ann)
    # fun.diff_gene_pathway(ann)

def draw_genome(focus_genome, focus_s, focus_e, feats):
    # out = open("/mnt/c/Users/swang66/Documents/genome_plot.csv", 'w')
    # print ("molecule, gene, start, end, strand, direction", file = out)
    
    f = open(gff)
    for line in f:
        array = line.split("\t")
        genome = array[0]
        if genome != focus_genome:
            continue
        info = array[8]
        match = re.search("Name=(.*?);", info)
        if not match:
            # continue
            name = array[2]
            
        else:
            name = match.group(1)
        if name == "CDS":
            continue
        g_type = array[2]
        
        start = int(array[3])
        end = int(array[4])
        if array[6] == "+":
            strand = "forward"
            direction = 1
        elif array[6] == "-":
            strand = "reverse"
            direction = -1
        if start >= focus_s and end <= focus_e:
            # print (genome, name, start, end, strand, direction, sep=",", file = out)
            feats.append([genome, start - focus_s, end-focus_s, array[6], name, g_type])          
    f.close()
    return feats

def draw_HGT_genome_pair(ann):
    with open("diff_bkps_genomes", "rb") as fp:
        genome_pairs = pickle.load(fp)  

    seqs, feats, links = [], [], []
    length = 100000
    line_len = 1
    region_dict = {"GUT_GENOME143131_1":13, "GUT_GENOME143505_1":14, "GUT_GENOME147678_1":13,\
     "GUT_GENOME096063_2":13, "GUT_GENOME144544_1":1, "GUT_GENOME147598_1":17}

    index = 0
    genome_dict = {}
    print (genome_pairs.keys())
    for edge in genome_pairs:
        print (edge)
        array = edge.split("&")
        genome_1 = array[0]
        genome_2 = array[1]

        one = region_dict[genome_1]
        two = region_dict[genome_2]
        
        plot_start_1 = length * one
        plot_end_1 = length * (one + 1)
        plot_start_2 = length * two
        plot_end_2 = length * (two + 1)


        if genome_1 not in genome_dict:
            seqs.append([genome_1, plot_end_1 - plot_start_1])
        if genome_2 not in genome_dict:
            seqs.append([genome_2, plot_end_2 - plot_start_2])
        genome_dict[genome_1] = 1
        genome_dict[genome_2] = 1
        feats = draw_genome(genome_1, plot_start_1, plot_end_1, feats)
        feats = draw_genome(genome_2, plot_start_2, plot_end_2, feats)
        # print (len(sort_select_edges[edge]))
        g1 = {}
        for bk_pair in genome_pairs[edge]:
            array = bk_pair.split("&")

            bk1 = array[0]
            genome_1 = bk1.split("|")[0]
            genome_1_pos = int(bk1.split("|")[1]) * window

            bk2 = array[1]
            genome_2 = bk2.split("|")[0]
            genome_2_pos = int(bk2.split("|")[1]) * window

            
            if genome_1_pos >= plot_start_1 and genome_1_pos <= plot_end_1 and \
                genome_2_pos >= plot_start_2 and genome_2_pos <= plot_end_2:
                links.append([genome_1, genome_1_pos-plot_start_1, genome_1_pos+line_len-plot_start_1,\
                 genome_2, genome_2_pos-plot_start_2, genome_2_pos+line_len-plot_start_2 ])  
                # print (genome_1_pos, genome_2_pos)
            c1 = int(genome_1_pos/100000)
            c2 = int(genome_2_pos/100000)
            k = str(c1) + "|" + str(c2) 
            # if k == '13|278':
            #     print (genome_1, genome_1_pos, genome_1_pos+5, genome_2, genome_2_pos, genome_2_pos+5)
            if k not in g1:
                g1[k] = 0
            g1[k] += 1
        sort_select_edges = sorted(g1.items(), key=lambda item: item[1], reverse = True)
        for i in range(50):
            print (edge, sort_select_edges[i])
        
        # index += 1
        # if index > 1:
        #     break  

        # break
    seqs = pd.DataFrame(seqs, columns = ["seq_id", "length"])
    feats = pd.DataFrame(feats, columns = ["seq_id", "start", "end", "strand", "name", "type"])
    links = pd.DataFrame(links, columns = ["seq_id", "start", "end", "seq_id2", "start2", "end2"])
    seqs.to_csv('/mnt/c/Users/swang66/Documents/seqs.csv', sep='\t')
    feats.to_csv('/mnt/c/Users/swang66/Documents/feats.csv', sep='\t')
    links.to_csv('/mnt/c/Users/swang66/Documents/links.csv', sep='\t')

def get_split_reads_cutoff(g): # g is the number of reads in the specific sample
    # prior_a = 25 # prior probability
    prior_b = 50000000#63333330  # prior probability
    given_n = 42648185 # mean number of reads among samples
    # # n = 182534663 # max number of reads 
    # given_r = 2 # the cutoff with n reads  4
    prior_a = 20 # prior probability
    given_r = 2 # the cutoff with n reads  4

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


class Micro_homo():

    def __init__(self):
        self.all_data, self.sample_num = None, None
        
        self.ref = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"
        self.workdir = "/mnt/d/breakpoints/HGT/micro_homo/"
        self.ref_fasta = Fasta(self.ref)        
        self.cutoff = 6
        self.min_score = 8
        self.all_bkp = {}
        self.all_bkp_sample_num = {}
        self.tole_diff = 10
        self.shortest_len = 4
        self.sam_number = 50000

        self.random_homo_seq_count = {}
        self.hgt_homo_seq_count = {}

    def get_reverse_complement_seq(self, sequence):
        sequence = sequence[::-1]
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
        string = sequence.translate(trantab)
        return string

    def extract_ref_seq(self, scaffold_name, start, end):
        if start < 1:
            start = 1
        # print (scaffold_name, start, end)
        return self.ref_fasta[scaffold_name][start:end].seq

    def prof_sort_bkp(self):
        self.all_data, self.sample_num = get_data()
        print ("data is loaded")
        sample_index = 0
        for sample in self.all_data:
            for bkp in sample.bkps:
                my_key = "|".join([bkp.from_ref, str(bkp.from_bkp), bkp.to_ref, str(bkp.to_bkp)])
                if my_key not in self.all_bkp:
                    self.all_bkp[my_key] = bkp
                    self.all_bkp_sample_num[my_key] = 0
                self.all_bkp_sample_num[my_key] += 1
            sample_index += 1
            # print (sample_index)
        print (len(self.all_bkp))
        # with open("all_of_bkp", "wb") as fp:
        #     pickle.dump(self.all_bkp, fp)
        # sorted_count = sorted(self.all_bkp_sample_num.items(), key=lambda item: item[1], reverse = True)
        # for i in range(50):
        #     my_key = sorted_count[i][0]
        #     sample_num = sorted_count[i][1]
        #     bkp = self.all_bkp[my_key]
        #     score = self.for_each_bkp(bkp, '')
        #     print (my_key, sample_num, score, bkp.from_strand, bkp.to_strand)

    def load_all_bkp(self):
        with open("all_of_bkp", "rb") as fp:
            self.all_bkp = pickle.load(fp)
        index = 0
        hit_num = 0
        bkp_key_list = list(self.all_bkp.keys())
        random.shuffle(bkp_key_list)
        for bkp_key in bkp_key_list[:self.sam_number]:
            bkp = self.all_bkp[bkp_key]
            # if bkp.cross_split_reads < 5:
            #     continue
            score = self.for_each_bkp(bkp, index)
            if score == -1:
                continue
            if score >= self.min_score:
                hit_num += 1
            index += 1
            if index % 1000 == 0:
                print (index, hit_num, hit_num/index)            
        print ("hit rate:", index, hit_num, hit_num/index)
        return hit_num/index

    def main(self):
        self.all_data, self.sample_num = get_data()
        print ("data is loaded")
        index = 0
        hit_num = 0
        sample_num = 0
        for sample in self.all_data:
            for bkp in sample.bkps:
                score = self.for_each_bkp(bkp, index)
                if score >= self.min_score:
                    hit_num += 1
                index += 1
                if index % 1000 == 0:
                    print (sample_num, index, hit_num, hit_num/index)
            sample_num += 1
            #     if index > 100:
            #         break
            # break
        print ("hit rate:", hit_num/index)

    def check_overlap(self, start_end_positions):
        flag = False
        # if start_end_positions[0][0] >= start_end_positions[1][0] and start_end_positions[0][0] <= start_end_positions[1][1]:
        #     flag = True
        # elif start_end_positions[0][1] >= start_end_positions[1][0] and start_end_positions[0][1] <= start_end_positions[1][1]:
        #     flag = True

        pos_diff = abs(start_end_positions[0][0] - start_end_positions[1][0])
        if pos_diff <= self.tole_diff: 
            flag = True
        return flag

    def get_mic_homo(self, seq1, seq2):
        seq1 = DNA(seq1)
        seq2 = DNA(seq2)
        # if seq1.frequencies("N")["N"] > 0 or seq2.frequencies("N")["N"] > 0:
        if seq1.frequencies("N")["N"] == len(seq1) or seq2.frequencies("N")["N"] == len(seq2):
            return -1
        elif len(seq1) != len(seq2):
            return -1
        alignment, score, start_end_positions = local_pairwise_align_ssw(seq1, seq2, gap_open_penalty=10, mismatch_score=-10)
        array = str(alignment).split("\n")
        match_seq1 = array[-2]
        match_seq2 = array[-1]
        
        score = len(match_seq1) # use match length as score
        vis_seq_1 = self.visulaize(start_end_positions[0], match_seq1)
        vis_seq_2 = self.visulaize(start_end_positions[1], match_seq2)
        # print (seq1)
        # print (seq2)
        # print (vis_seq_1)
        # print (vis_seq_2)

        if score >= self.min_score and self.check_overlap(start_end_positions):


            return score
        else:
            return 0

    def for_each_bkp(self, bkp, index):
        bkp.from_bkp -= 1
        from_seq = self.extract_ref_seq(bkp.from_ref, bkp.from_bkp-self.cutoff, bkp.from_bkp+self.cutoff)
        to_seq = self.extract_ref_seq(bkp.to_ref, bkp.to_bkp-self.cutoff, bkp.to_bkp+self.cutoff)
        if bkp.from_strand == "-":
            from_seq = self.get_reverse_complement_seq(from_seq)      
        if bkp.to_strand == "-":
            to_seq = self.get_reverse_complement_seq(to_seq)
        if re.search(">", from_seq) or re.search(">", to_seq):
            return -1
        # test = ["GUT_GENOME018982_31", "GUT_GENOME239728_21"]
        # if bkp.from_ref in test and bkp.to_ref in test:
        #     print (bkp.from_ref, bkp.to_ref, from_seq, to_seq, self.find_mh(from_seq, to_seq))
        # score = self.get_mic_homo(from_seq, to_seq)
        # return score
        if self.find_mh(from_seq, to_seq, "hgt"):
            return 100
        else:
            return 0

    def random_seq(self):
        chroms = list(self.ref_fasta.keys())
        chroms_num = len(chroms)
        hit_num = 0
        # for i in range(1000):
        index = 0
        for j in range(self.sam_number):
            seq_list = []
            for j in range(2):
                chrom = list(self.ref_fasta.keys())[np.random.randint(chroms_num)]
                chrom_len = len(self.ref_fasta[chrom])
                locus = np.random.randint(self.cutoff, chrom_len-self.cutoff)
                seq = self.extract_ref_seq(chrom, locus-self.cutoff, locus+self.cutoff)
                seq_list.append(seq)
                # print (chrom, locus)
            # score = self.get_mic_homo(seq_list[0], seq_list[1])
            if self.find_mh(seq_list[0], seq_list[1], "random"):
                score =  100
            else:
                # print (seq_list[0], seq_list[1])
                score =  0
            if score == -1:
                # print (score)
                continue
            if score >= self.min_score:
                hit_num += 1
            index += 1
            if index % 1000 == 0:
                print (index, hit_num, hit_num/(index))
        print ("random hit rate", hit_num, hit_num/(index))
        return hit_num/(index)

    def visulaize(self, interval, match_seq):
        vis_seq = ''
        for i in range(0, interval[0]):
            vis_seq += "-"
        vis_seq += match_seq
        # for i in range(interval[0], interval[1]+1):
        #     vis_seq += "*"
        for i in range(interval[1]+1, self.cutoff*2):
            vis_seq += "-"
        return vis_seq

    def output_fasta(self, from_seq, to_seq, index):
        seq1 = SeqRecord(Seq(from_seq),
                        id="from_seq")
        seq2 = SeqRecord(Seq(to_seq),
                        id="to_seq")
        SeqIO.write(seq1, self.workdir+"from_seq_%s.fasta"%(index), "fasta")
        SeqIO.write(seq2, self.workdir + "to_seq_%s.fasta"%(index), "fasta")
        command = "blastn -evalue 0.05 -word_size 4 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -dust no\
            -query %s -subject %s -out %s"%(self.workdir+"from_seq_%s.fasta"%(index),\
             self.workdir + "to_seq_%s.fasta"%(index), self.workdir + "blast_%s.out"%(index))
        print (command)

    def find_mh(self, seq1, seq2, method):
        flag = False
        seq1 = seq1.lower()
        seq2 = seq2.lower()
        # print (self.shortest_len)
        for i in range(len(seq1)-self.shortest_len+1):
            search_seq = seq1[i:i+self.shortest_len]
            mat = re.search(search_seq, seq2)
            # print (search_seq, seq2)
            if mat:
                if method == "random":
                    if search_seq not in self.random_homo_seq_count:
                        self.random_homo_seq_count[search_seq] = 0
                    self.random_homo_seq_count[search_seq] += 1

                if method == "hgt":
                    if search_seq not in self.hgt_homo_seq_count:
                        self.hgt_homo_seq_count[search_seq] = 0
                    self.hgt_homo_seq_count[search_seq] += 1

                mat_locus = mat.span()[0]
                # print (abs(mat_locus - i), self.tole_diff, "true")
                if abs(mat_locus - i) <= self.tole_diff :
                    # print (abs(mat_locus - i), "true")
                    flag =  True
        return flag

def microhomology_freq_compare():
    data = []
    for diff in range(1):
        for score in range(3, 11):
            mic = Micro_homo()
            mic.tole_diff = diff
            mic.shortest_len = score
            print ("score is", score, "diff is", diff)
            hgt_freq = mic.load_all_bkp()
            random_freq = mic.random_seq()
            data.append([diff, score, hgt_freq, "HGT"])
            data.append([diff, score, random_freq, "Random"])
            print ("-----------------------")
    df = pd.DataFrame(data, columns = ["diff", "length", "frequency", "group"])
    df.to_csv('/mnt/d/breakpoints/script/analysis/microhomo_freq.csv', sep='\t')

def count_seq_freq(mic, hgt_homo_seq_count):
    new_dict = {}
    print (len(hgt_homo_seq_count))
    for key in sorted(hgt_homo_seq_count):
        # print (key)
        if key not in new_dict:
            if mic.get_reverse_complement_seq(key) in new_dict:
                new_dict[mic.get_reverse_complement_seq(key)] += hgt_homo_seq_count[key]
            else:
                new_dict[key] = hgt_homo_seq_count[key]
        else:
            new_dict[key] += hgt_homo_seq_count[key]
            # print ("impossible")
    print (len(new_dict))
    return new_dict
        


def microhomology_seq_freq():
    # find homology pattern
    mic = Micro_homo()
    mic.tole_diff = 0
    mic.shortest_len = 6
    mic.sam_number = 50000
    data = []
    
    
    hgt_freq = mic.load_all_bkp()
    mic.hgt_homo_seq_count = count_seq_freq(mic, mic.hgt_homo_seq_count)
    sorted_count = sorted(mic.hgt_homo_seq_count.items(), key=lambda item: item[1], reverse = True)
    hgt_value_list = list(mic.hgt_homo_seq_count.values())
    sum_hgt_value = sum(hgt_value_list)
    # for key in sorted(mic.hgt_homo_seq_count):
    for item in sorted_count:
        key = item[0]
        random_num, hgt_num = 0, 0
        if key in mic.hgt_homo_seq_count:
            hgt_num = mic.hgt_homo_seq_count[key]
        if hgt_num > 5:
            print (key, round(hgt_num/sum_hgt_value, 6), hgt_num)
        data.append([key.upper(), round(hgt_num/sum_hgt_value, 6), hgt_num, "HGT"])

    # for key in sorted(mic.hgt_homo_seq_count):
    #     mic.random_homo_seq_count[key] = 0
    # random_freq = mic.random_seq()
    # mic.random_homo_seq_count = count_seq_freq(mic, mic.random_homo_seq_count)
    # random_value_list = list(mic.random_homo_seq_count.values())
    # sum_random_value = sum(random_value_list)
    # for key in sorted(mic.random_homo_seq_count):
    #     random_num, hgt_num = 0, 0
    #     if key in mic.random_homo_seq_count:
    #         random_num = mic.random_homo_seq_count[key]
    #     print (key, round(random_num/sum_random_value, 6), random_num)
    #     data.append([key.upper(), round(random_num/sum_random_value, 6), hgt_num, "random"])
    
    # df = pd.DataFrame(data, columns = ["seq", "frequency", "number", "group"])
    # df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files/micro_seq_pattern.csv', sep=',')
    

if __name__ == "__main__":
    print ("start...")
    gff = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.gff"

    plot_genomes = ["GUT_GENOME143131_1&GUT_GENOME143505_1", "GUT_GENOME143505_1&GUT_GENOME147678_1", \
    "GUT_GENOME096063_2&GUT_GENOME147678_1", 'GUT_GENOME096063_2&GUT_GENOME144544_1', 'GUT_GENOME144544_1&GUT_GENOME147598_1']


    # microhomology_freq_compare()
    # microhomology_seq_freq()
    
    main()
    # fa = Family()
    # fa.compare_cccur_in_each_cohort()
    # fa.get_marker_tree()
    # fa.plot_all_diff_hgts()
    # draw_genome("GUT_GENOME095938_1", 251800, 2189200)



    # GUT_GENOME143131_1&GUT_GENOME143505_1  22300
    # GUT_GENOME143505_1&GUT_GENOME147678_1  15952
    # GUT_GENOME096063_2&GUT_GENOME147678_1  11520
    # 'GUT_GENOME096063_2&GUT_GENOME144544_1', 7245
    # 'GUT_GENOME144544_1&GUT_GENOME147598_1', 4177
    print ("Finished")