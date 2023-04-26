import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import networkx as nx
from scipy.stats import mannwhitneyu
from scipy import stats
import scipy 
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
from random import shuffle
import sklearn
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler

from mechanism_taxonomy import Taxonomy

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

    def __init__(self, bkp_list, ID, pheno):
        self.bkps = bkp_list
        self.ID = ID  
        self.cohort = pheno[0]  
        self.disease = pheno[1]  
        self.full_disease = pheno[2].split(";")

    def get_HGT_matrix(self, level, edge_num):
        nodes_index = {}
        bkp_score = {}
        choose_edge = {}
        for bkp in self.bkps:
            edge = get_tag(bkp, level)
            # print (edge)
            support_ratio = bkp.cross_split_reads
            if edge not in bkp_score:
                bkp_score[edge] = support_ratio
            if support_ratio > bkp_score[edge]:
                bkp_score[edge] = support_ratio
        sort_bkp_score = sorted(bkp_score.items(), key=lambda item: item[1], reverse = True)
        choose_edge = {}
        total_edge_num = len(sort_bkp_score)
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
        if total_edge_num < edge_num:
            return [0, 0,0, 0,0,0], total_edge_num
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

class Data_load():

    def __init__(self):
        self.sample_obj_list = []
        
    def read_samples(self):

        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        # os.system(f"ls {tgs_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        # os.system(f"ls {wenkui_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            if len(my_bkps) > 0 and sra_id in phenotype_dict:
                sample = Sample(my_bkps, sra_id, phenotype_dict[sra_id])
                self.sample_obj_list.append(sample)
        print ("data is loaded.")
        
    def read_bkp(self, bkp_file):
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
        return my_bkps

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

class Marker():

    def __init__(self, group1, group2):
        self.group1 = group1
        self.group2 = group2
        self.all_HGTs = {}
        self.sample_count = {self.group1:0, self.group2:0}
        self.marker_num = marker_num
        self.markers = {}

    def extract_HGT(self, sample_obj_list):

        for sample in sample_obj_list:

            if len(sample.full_disease) != 1:
                continue
            if sample.disease == "control" and sample.full_disease[0] != "healthy":
                continue
            if sample.disease == '':
                continue
            if sample.disease ==  self.group1 or self.group1 in sample.full_disease:
                index = 0
            elif sample.disease ==self.group2 or self.group2 in sample.full_disease:
                index = 1
            else:
                continue
            if sample.disease not in self.sample_count:
                print (sample.ID, sample.disease, self.sample_count, sample.disease, sample.full_disease)
            self.sample_count[sample.disease] += 1

            sample_dict = {}
            for bkp in sample.bkps:
                if bkp.hgt_tag == "NA":
                    continue
                if bkp.hgt_tag in sample_dict:
                    continue
                if bkp.hgt_tag not in self.all_HGTs:
                    self.all_HGTs[bkp.hgt_tag] = {self.group1:0, self.group2:0}
                sample_dict[bkp.hgt_tag] = 1
                self.all_HGTs[bkp.hgt_tag][sample.disease] += 1
        filtered_HGT = {}
        print ("Bkp num in the two groups", len(self.all_HGTs))
        for hgt_tag in self.all_HGTs:
            if (self.all_HGTs[hgt_tag][self.group1] + self.all_HGTs[hgt_tag][self.group2])/(self.sample_count[self.group1]+self.sample_count[self.group2]) < cutoff:
                pass
            else:
                filtered_HGT[hgt_tag] = self.all_HGTs[hgt_tag]
        self.all_HGTs = filtered_HGT
        print ("Filtered bkp num in the two groups", len(self.all_HGTs))
        print ("%s num: %s, %s num %s."%(self.group1, self.sample_count[self.group1],self.group2, self.sample_count[self.group2]))

    def select_diff_HGT(self):
        hgt_p_value_dict = {}
        for hgt_tag in self.all_HGTs:
            g1_array = self.all_HGTs[hgt_tag][self.group1] * [1] + [0] * (self.sample_count[self.group1] - self.all_HGTs[hgt_tag][self.group1])
            g2_array = self.all_HGTs[hgt_tag][self.group2] * [1] + [0] * (self.sample_count[self.group2] - self.all_HGTs[hgt_tag][self.group2])
            U1, p = mannwhitneyu(g1_array, g2_array)
            hgt_p_value_dict[hgt_tag] = p
            # print (hgt_tag, p)
        sorted_hgt_p_value = sorted(hgt_p_value_dict.items(), key=lambda item: item[1], reverse = False)
        for i in range(self.marker_num):
            marker = sorted_hgt_p_value[i][0]
            p = sorted_hgt_p_value[i][1]
            self.markers[marker] = i
            # print ("marker", marker, p)

    def training(self, sample_obj_list):
        data, label = [], []
        for sample in sample_obj_list:
            if len(sample.full_disease) != 1:
                continue
            if sample.disease == '':
                continue
            if sample.disease ==  self.group1 or self.group1 in sample.full_disease:
                index = 0
            elif sample.disease ==self.group2 or self.group2 in sample.full_disease:
                index = 1
            else:
                continue
            
            marker_value = [0] * self.marker_num
            for bkp in sample.bkps:
                if bkp.hgt_tag in self.markers:
                    marker_value[self.markers[bkp.hgt_tag]] = 1
            data.append(marker_value)
            label.append(index)

        print ("before", len(label))
        #### Apply oversampling to balance the input data
        # oversample = SMOTE()
        # data, label = oversample.fit_resample(data, label)
        
        rus = RandomUnderSampler(random_state=42)
        data, label = rus.fit_resample(data, label)

        print ("after", len(label))


        rfc = RandomForestClassifier(class_weight='balanced', n_estimators=100)
        auc_scores = cross_val_score(rfc, data, label, cv=5, scoring='roc_auc')
        # Print the mean accuracy of the model across all folds
        print(self.group1, self.group2, auc_scores, "Mean AUC-ROC score:", auc_scores.mean())
        print ("<<<<<<<<<<<<<\n")
        return auc_scores.mean()




if __name__ == "__main__":

    abun_cutoff = 1e-7  #1e-7
    cutoff = 0.1
    level = 5
    marker_num = 20
    replication = 1

    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    phenotype_dict = read_phenotype()
    taxonomy = Taxonomy()
    dat = Data_load()
    dat.read_samples()

    group1 = "control"
    group2 = "IBD"

    data = []
    for marker_num in range(5, 100, 5):
        group_auc = []
        group_list = ["control", "CRC", "T2D",  "IBD"]
        # group_list = ["control", "CRC", "adenoma", "IGT", "T2D", "acute_diarrhoea",  "IBD"]
        for i in range(len(group_list)):
            for j in range(i+1, len(group_list)):
                replicate_result = []
                for z in range(replication):
                    group1 = group_list[i]
                    group2 = group_list[j]

                    mar = Marker(group1, group2)
                    shuffle(dat.sample_obj_list)
                    mar.extract_HGT(dat.sample_obj_list)
                    mar.select_diff_HGT()
                    mean_auc = mar.training(dat.sample_obj_list)
                    replicate_result.append(mean_auc)
                    data.append([marker_num, mean_auc])
                group_auc.append(np.mean(replicate_result))
        print ("#########", marker_num, np.mean(group_auc))
        data.append([marker_num, np.mean(group_auc)])
        print ("---------------\n")

    df = pd.DataFrame(data, columns = ["Feature_number", "AUC"])
    df.to_csv('/mnt/d/R_script_files/classifier_feature_num.csv', sep=',')

