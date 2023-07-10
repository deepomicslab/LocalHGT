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
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

from collections import defaultdict
from collections import Counter
import pickle
from Bio.KEGG import REST
from ete3 import Tree

from mechanism_taxonomy import Taxonomy
from analyze_transfer_gene import Annotation

level_list = ["phylum", "class", "order", "family", "genus", "species", "genome"]



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

        self.from_family = self.from_ref_lineage.split(";")[4]
        self.to_family = self.to_ref_lineage.split(";")[4]

        taxa1 = self.from_ref_lineage.split(";")[level]
        taxa2 = self.to_ref_lineage.split(";")[level]
        if taxa1[1:] == "__" or taxa2[1:] == "__":
            self.hgt_tag = "NA"
        else:
            self.hgt_tag = "&".join(sorted([taxa1, taxa2]))
        self.bk1_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) 
        self.bk2_tag = self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        
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
        self.cohort_set = set()
        
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
                self.cohort_set.add(phenotype_dict[sra_id][0])
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

    def __init__(self, group1, group2, sample_obj_list):
        self.group1 = group1
        self.group2 = group2

        self.near = 5000 #5000   ## up and down this distance from the breakpoint position

        self.all_HGTs = {}
        self.all_HGTs_family = defaultdict(list) # hgt_tag:[family name]
        self.all_breakpoints = {}
        self.diff_breakpoints = {}

        self.sample_count = None 
        self.markers = {}
        self.split_data_dict = {}
        self.sample_obj_list = sample_obj_list

        self.selected_samples = {}
        self.enough_sample_flag = self.select_sample()

        self.diff_bk = defaultdict(list)
        self.normal_bk = defaultdict(list)



    def extract_HGT(self):
        self.sample_count = {self.group1:0, self.group2:0}
        for sample in self.sample_obj_list:
            if sample.ID not in self.selected_samples:
                continue
            self.sample_count[sample.disease] += 1

            sample_dict = {}
            bkp_dict = {}
            for bkp in sample.bkps:
                if bkp.bk1_tag not in self.all_breakpoints:
                    self.all_breakpoints[bkp.bk1_tag] = {self.group1:0, self.group2:0}
                if bkp.bk2_tag not in self.all_breakpoints:
                    self.all_breakpoints[bkp.bk2_tag] = {self.group1:0, self.group2:0}

                if bkp.bk1_tag not in bkp_dict:
                    self.all_breakpoints[bkp.bk1_tag][sample.disease] += 1
                if bkp.bk2_tag not in bkp_dict:
                    self.all_breakpoints[bkp.bk2_tag][sample.disease] += 1

                bkp_dict[bkp.bk1_tag] = 1
                bkp_dict[bkp.bk2_tag] = 1

                if bkp.hgt_tag == "NA":
                    continue

                if bkp.hgt_tag not in sample_dict:
                    if bkp.hgt_tag not in self.all_HGTs:
                        self.all_HGTs[bkp.hgt_tag] = {self.group1:0, self.group2:0}               
                    self.all_HGTs[bkp.hgt_tag][sample.disease] += 1
                    self.all_HGTs_family[bkp.hgt_tag] = sorted([bkp.from_family, bkp.to_family])
                sample_dict[bkp.hgt_tag] = 1

        print ("pair num", len(self.all_HGTs), "breakpoint num", len(self.all_breakpoints))

    def select_diff_HGT(self, name): #genus pair
        
        data = []
        for hgt_tag in self.all_HGTs:
            # if self.all_HGTs_family[hgt_tag][0] !=  "f__Lachnospiraceae" and self.all_HGTs_family[hgt_tag][1] !=  "f__Lachnospiraceae":
            #     continue
            a = self.all_HGTs[hgt_tag][self.group1]
            b = self.sample_count[self.group1] - self.all_HGTs[hgt_tag][self.group1]
            c = self.all_HGTs[hgt_tag][self.group2]
            d = self.sample_count[self.group2] - self.all_HGTs[hgt_tag][self.group2]

            group1_freq = a/(a+b)
            group2_freq = c/(c+d)
            oddsratio, p_value = fisher_exact([[a, b], [c, d]])
            data.append([hgt_tag, p_value, oddsratio, a, group1_freq, group2_freq])

        df = pd.DataFrame(data, columns = ["genus_pair", "p_value", "oddsratio", "gp_num", self.group1, self.group2])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected

        filtered_df = df[df['p.adj'] < 0.05]

        self.print_diff_pair(filtered_df, name)

        group1_df = filtered_df[filtered_df[self.group1] > filtered_df[self.group2]]
        group2_df = filtered_df[filtered_df[self.group1] < filtered_df[self.group2]]
        print ("differential genus pairs enriched for each group", len(group1_df), len(group2_df))
        # print (filtered_df)
        print ("\n<<<<<<<<<<<<<<<%s enriched:"%(self.group1), len(group1_df))
        self.classify_breakpoints_pair(group1_df, f"/mnt/d/R_script_files/group1_{name}.csv")

        print ("\n<<<<<<<<<<<<<<<%s enriched:"%(self.group2), len(group2_df))
        self.classify_breakpoints_pair(group2_df, f"/mnt/d/R_script_files/group2_{name}.csv")

    def print_diff_pair_bk(self, diff_pairs):
        colors = ["#F7C530", "#95CC5E", "#D0DFE6","pink","#4169E1","#F0DB4F","#226F54","#FF7F50","#FFA500","#6EE2FF","lightgrey"]*4
        genus2family = get_taxa_dict()
        family_index = {}

        index_dict = {}
        for pair in diff_pairs:
            array = pair.split("&")
            for a in array:
                
                family = genus2family[a]
                index_dict[a] = family
                if family not in family_index:
                    family_index[family] = len(family_index)

        print ("family count", len(family_index))
        index_dict = dict(sorted(index_dict.items(), key=lambda x: x[1]))
        names = list(index_dict.keys())

        index_dict = {}
        my_color = []
        family_list = []
        for i in range(len(names)):
            index_dict[names[i]] = i
            family = genus2family[names[i]]
            color = colors[family_index[family]]
            family_list.append(family)
        my_color.append(color)
            
        print (index_dict)
        print (family_list)
        print (my_color)
        matrix = np.zeros((len(index_dict), len(index_dict)))
        for pair in diff_pairs:
            array = pair.split("&")
            matrix[index_dict[array[0]]][index_dict[array[1]]] = 1    
            matrix[index_dict[array[1]]][index_dict[array[0]]] = 1  
       
        df = pd.DataFrame(matrix, index=names, columns=names)  
        print (df)
        df.to_csv("/mnt/d/R_script_files/diff_pair_matrix.csv", sep=',')
        with open('/mnt/d/R_script_files/diff_pair_matrix_color.csv', 'w') as f:
            # iterate over the list and write each element to a new line in the file
            for item in my_color:
                f.write("%s," % item)

    def print_diff_pair(self, filtered_df, name):

        # colors = ["#F7C530", "#95CC5E", "#D0DFE6","pink","#4169E1","#F0DB4F","#226F54","#FF7F50","#FFA500","#6EE2FF","lightgrey"]*5
        # colors = ["red", "green", "blue", "purple", "yellow", "orange", "pink", "brown", "black", "gray"]*10
        colors = ["red", "green", "blue", "orange", "purple", "pink", "yellow", "brown", "gray", "black", 
                "white", "cyan", "magenta", "darkolivegreen", "navy", "teal", "maroon", "silver", "gold", "turquoise", 
                "indigo", "lavender", "coral", "violet", "beige", "khaki", "salmon", "tan", "fuchsia", 
                "lime", "aqua", "chartreuse", "crimson", "orchid", "plum"] * 4
        # print (colors)
        genus2family = get_taxa_dict()
        family_index = {}

        index_dict = {}
        link_times = defaultdict(int)
        for index, row in filtered_df.iterrows():
            # diff_pairs[row["genus_pair"]] = row
            pair = row["genus_pair"]
        # for pair in diff_pairs:
            array = pair.split("&")
            for a in array:
                family = genus2family[a]
                index_dict[a] = family
                link_times[a] += 1
                if family not in family_index:
                    family_index[family] = len(family_index)

        print ("involved genus count", len(index_dict))
        print ("involved family count", len(family_index))
        print (family_index)


        index_dict = dict(sorted(index_dict.items(), key=lambda x: x[1]))
        names = list(index_dict.keys())

        family_counts = defaultdict(int)
        index_dict = {}
        data = []
        hist_data = []
        for i in range(len(names)):
            index_dict[names[i]] = i
            family = genus2family[names[i]]
            family_counts[family] += 1
            color = colors[family_index[family]]
            print (names[i], 3, family, color)
            data.append([names[i], 3, family, color])
            hist_data.append([names[i], 0, 3, link_times[names[i]]])

        family_counts = dict(sorted(family_counts.items(), key=lambda x: x[1]))
        print ("family count", len(family_counts), family_counts)
        print ("genus link times",link_times )

        bed1, bed2 = [], []
        bed3, bed4 = [], []
        for index, row in filtered_df.iterrows():
            pair = row["genus_pair"]
            array = pair.split("&")

            if row[self.group1] > row[self.group2]:
                bed1.append([array[0], 1.5])
                bed2.append([array[1], 1.5])
            else:
                bed3.append([array[0], 1.5])
                bed4.append([array[1], 1.5])               
       
        df = pd.DataFrame(data, columns=["genus", "length", "family", "color"] )  
        df.to_csv("/mnt/d/R_script_files/diff_pair_genus_%s.csv"%(name), sep=',')

        df = pd.DataFrame(hist_data, columns=["genus", "start", "end", "value1"] )  
        df.to_csv("/mnt/d/R_script_files/diff_pair_genus_hist_%s.csv"%(name), sep=',', index=False)

        df = pd.DataFrame(bed1, columns=["chrom", "mid"] )  
        df.to_csv("/mnt/d/R_script_files/diff_pair_genus_bed1_%s.csv"%(name), sep=',', index=False)
    
        df = pd.DataFrame(bed2, columns=["chrom", "mid"] )  
        df.to_csv("/mnt/d/R_script_files/diff_pair_genus_bed2_%s.csv"%(name), sep=',', index=False)

        df = pd.DataFrame(bed3, columns=["chrom", "mid"] )  
        df.to_csv("/mnt/d/R_script_files/diff_pair_genus_bed3_%s.csv"%(name), sep=',', index=False)
    
        df = pd.DataFrame(bed4, columns=["chrom", "mid"] )  
        df.to_csv("/mnt/d/R_script_files/diff_pair_genus_bed4_%s.csv"%(name), sep=',', index=False)

    def get_tree(self):
        involved_genomes = {}
        for index, row in filtered_df.iterrows():
            diff_pairs[row["genus_pair"]] = row
            array = row["genus_pair"].split("&")
            involved_genomes.add(array[0])
            involved_genomes.add(array[1])
        # tree.prune(extracted_genome) # only keep the nodes involved

    def select_diff_breakpoint(self):
        data = []
        for bkp_tag in self.all_breakpoints:

            a = self.all_breakpoints[bkp_tag][self.group1]
            b = self.sample_count[self.group1] - self.all_breakpoints[bkp_tag][self.group1]
            c = self.all_breakpoints[bkp_tag][self.group2]
            d = self.sample_count[self.group2] - self.all_breakpoints[bkp_tag][self.group2]

            group1_freq = a/(a+b)
            group2_freq = c/(c+d)

            oddsratio, p_value = fisher_exact([[a, b], [c, d]])
            data.append([bkp_tag, p_value, oddsratio, a, group1_freq, group2_freq])

        df = pd.DataFrame(data, columns = ["breakpoint", "p_value", "oddsratio", "bk_num", self.group1, self.group2])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected

    
        filtered_df = df[df['p.adj'] < 0.05]
        # filtered_df = filtered_df[filtered_df[self.group1] > filtered_df[self.group2]]
        print ("differential breakpoints", len(filtered_df))
        diff_genomes = {}
        diff_breakpoints = {}

        for index, row in filtered_df.iterrows():

            pure_genome = "_".join(row["breakpoint"].split("_")[:-1]) 

            lineage = taxonomy.taxonomy_dict[pure_genome]
            print ( row["breakpoint"], lineage.split(";")[4], row["oddsratio"])
            family = lineage.split(";")[4]

            diff_genomes[pure_genome] = family
            diff_breakpoints[row["breakpoint"]] = row
        return diff_genomes, diff_breakpoints

    def classify_breakpoints(self):
        diff_genomes, diff_breakpoints = self.select_diff_breakpoint()

        diff_bkp_num, normal_bkp_num = 0, 0

        for sample in self.sample_obj_list:
            if sample.ID not in self.selected_samples:
                continue

            for bkp in sample.bkps:

                if bkp.bk1_tag in diff_breakpoints:
                    self.diff_bk[bkp.from_ref].append(bkp.from_bkp)
                    diff_bkp_num += 1
                else:
                    self.normal_bk[bkp.from_ref].append(bkp.from_bkp)
                    normal_bkp_num += 1

                if bkp.bk2_tag in diff_breakpoints:
                    self.diff_bk[bkp.to_ref].append(bkp.to_bkp)
                    diff_bkp_num += 1
                else:
                    self.normal_bk[bkp.to_ref].append(bkp.to_bkp)
                    normal_bkp_num += 1
        
        print ("genome num", len(self.diff_bk), len(self.normal_bk), "bkp num", diff_bkp_num, normal_bkp_num)

        diff_anno_dict = self.get_anno(self.diff_bk)
        normal_anno_dict = self.get_anno(self.normal_bk)

        # print (diff_anno_dict)
        # print (normal_anno_dict)

        enrichment_analysis_cog(diff_anno_dict["cog"], normal_anno_dict["cog"])

        print ("ko num", len(diff_anno_dict["ko"]), len(normal_anno_dict["ko"]))
        if len(diff_anno_dict["ko"]) > 0 and len(normal_anno_dict["ko"]) > 0:
            enrichment_analysis_kegg(diff_anno_dict["ko"], normal_anno_dict["ko"])

    def compare_all_breakpoints(self, all_kegg_diff): # directly compare all the breakpoints between groups 
        group1_bk = defaultdict(list)
        group2_bk = defaultdict(list)

        for sample in self.sample_obj_list:
            
            if sample.ID not in self.selected_samples:
                continue

            if sample.disease == self.group1:
                for bkp in sample.bkps:
                    group1_bk[bkp.from_ref].append(bkp.from_bkp)
                    group1_bk[bkp.to_ref].append(bkp.to_bkp)

            if sample.disease == self.group2:
                for bkp in sample.bkps:
                    group2_bk[bkp.from_ref].append(bkp.from_bkp)
                    group2_bk[bkp.to_ref].append(bkp.to_bkp)
        
        print ("genome num", len(group1_bk), len(group2_bk))

        diff_anno_dict = self.get_anno(group1_bk)
        normal_anno_dict = self.get_anno(group2_bk)

        # enrichment_analysis_cog(diff_anno_dict["cog"], normal_anno_dict["cog"])

        print ("ko num", len(diff_anno_dict["ko"]), len(normal_anno_dict["ko"]))
        if len(diff_anno_dict["ko"]) > 0 and len(normal_anno_dict["ko"]) > 0:
            filtered_df = enrichment_analysis_kegg(diff_anno_dict["ko"], normal_anno_dict["ko"])
            filtered_df.to_csv(all_kegg_diff, sep=',')

    def classify_breakpoints_pair(self, group1_df, kegg_output):

        diff_pair_bk = defaultdict(list)
        normal_pair_bk = defaultdict(list)

        diff_pairs = {}
        for index, row in group1_df.iterrows():
            diff_pairs[row["genus_pair"]] = row
            # print (row["genus_pair"], row[self.group1], row[self.group2])

        # diff_pairs = self.select_diff_HGT()
        diff_bkp_num, normal_bkp_num = 0, 0

        for sample in self.sample_obj_list:
            if sample.ID not in self.selected_samples:
                continue

            for bkp in sample.bkps:
                if bkp.hgt_tag in diff_pairs:
                    diff_pair_bk[bkp.from_ref].append(bkp.from_bkp)
                    diff_pair_bk[bkp.to_ref].append(bkp.to_bkp)
                    diff_bkp_num += 1
                else:
                    normal_pair_bk[bkp.from_ref].append(bkp.from_bkp)
                    normal_pair_bk[bkp.to_ref].append(bkp.to_bkp)
                    normal_bkp_num += 1 
        print ("pair num", len(diff_pairs), "bkp num", diff_bkp_num, normal_bkp_num)

        diff_anno_dict = self.get_anno(diff_pair_bk)
        normal_anno_dict = self.get_anno(normal_pair_bk)

        # if len(diff_anno_dict["cog"]) > 0 and len(normal_anno_dict["cog"]) > 0:
        #     enrichment_analysis_cog(diff_anno_dict["cog"], normal_anno_dict["cog"])
        # enrichment_analysis_scfa(diff_anno_dict["ko"], normal_anno_dict["ko"], self.group1, self.group2)

        # print ("ko num", len(diff_anno_dict["ko"]), len(normal_anno_dict["ko"]))
        if len(diff_anno_dict["ko"]) > 0 and len(normal_anno_dict["ko"]) > 0:
            filtered_df = enrichment_analysis_kegg(diff_anno_dict["ko"], normal_anno_dict["ko"])
            filtered_df.to_csv(kegg_output, sep=',')
            # self.examine_pathway(diff_anno_dict, diff_pair_bk)
            # self.examine_pathway(normal_anno_dict, normal_pair_bk)

    def examine_pathway(self, diff_anno_dict, bkp_dict):
        support_kos = get_support_kos(diff_anno_dict["ko"], "map00650")
        print ("%s kos support the pathway."%(len(support_kos)))
        support_genome = {}
        support_family = defaultdict(int)
        gene_dict = defaultdict(int)
        breakpoint_tag = {}

        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_site = site
                        locate_insert_flag = True
                        # break
                        # get all sites 
                        tag = genome + "&" + str(int(site/bin_size)) 
                        breakpoint_tag[tag] = 1

                if not locate_insert_flag:
                    continue

                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                flag = False
                if "KEGG" not in gene_anno_dict:
                    KEGG_list = []
                else:
                    KEGG_list = gene_anno_dict["KEGG"].split(",")
                for ko in KEGG_list:
                    if ko in support_kos:
                        flag = True

                if flag:  
                    support_genome[genome] = 1
                    pure_genome = "_".join(genome.split("_")[:-1])
                    lineage = taxonomy.taxonomy_dict[pure_genome]
                    family = lineage.split(";")[4]
                    support_family[family] += 1
                    if "Name" in gene_anno_dict:
                        gene_name = gene_anno_dict["Name"]
                    elif "gene" in gene_anno_dict:
                        gene_name = gene_anno_dict["gene"]
                    else:
                        gene_name = "NA"

                    gene_dict[gene_name] += 1
                    print (family, genome, gene_interval[0], gene_interval[1], locate_site, gene_name, sep = "\t")

        print ("genome", len(support_genome), len(support_family), len(gene_dict))
        print (support_family)
        print (gene_dict)

        # self.get_support_bkp(breakpoint_tag)

    def get_support_bkp(self, breakpoint_tag):

        for sample in self.sample_obj_list:
            if sample.ID not in self.selected_samples:
                continue

            for bkp in sample.bkps:
                if bkp.bk1_tag in breakpoint_tag or bkp.bk2_tag in breakpoint_tag:
                    print (sample.ID, sample.disease, bkp.from_ref, bkp.to_ref, bkp.from_family, bkp.to_family)

    def get_anno(self, bkp_dict):
        
        anno_dict = {"ko":[], "cog":[]}

        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]

                if "KEGG" not in gene_anno_dict:
                    KEGG_list = []
                else:
                    KEGG_list = gene_anno_dict["KEGG"].split(",")

                if "COG" not in gene_anno_dict:
                    cog = ""
                else:
                    cog = gene_anno_dict["COG"]

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    anno_dict["ko"] += KEGG_list
                    for i in range(len(cog)):
                        anno_dict["cog"] += [cog[i]]

        return anno_dict

    def select_sample(self):
        group1_num, group2_num = 0, 0
        for sample in self.sample_obj_list:

            if len(sample.full_disease) != 1:
                continue
            if sample.disease == "control" and sample.full_disease[0] != "healthy":
                continue
            if sample.disease == '':
                continue

            if focus_cohort != "all":
                if sample.cohort != focus_cohort:
                    continue

            if sample.disease ==  self.group1 or self.group1 in sample.full_disease:
                index = 0
                group1_num += 1
            elif sample.disease ==self.group2 or self.group2 in sample.full_disease:
                index = 1
                group2_num += 1
            else:
                continue 
            self.selected_samples[sample.ID] = index
        print ("sample counts: total num, group1 num, and group2 num", len(self.selected_samples), group1_num, group2_num)

        if group1_num == 0 or group2_num == 0:
            return False
        else:
            return True

def get_pathway_name_class(pathway_id):
    pathway_info = REST.kegg_get(f'path:{pathway_id}').read()
    pathway_name = pathway_info.split('\n')[1].split(';')[0].strip()
    pathway_name = pathway_name.replace("NAME", '')
    pathway_name = pathway_name.replace('"', '')
    pathway_name = pathway_name.strip()

    class_info = ["Other", "Other"]
    # print (pathway_info)
    for element in pathway_info.split('\n'):
        if re.search("CLASS", element):
            element = element.replace("CLASS", '').strip()
            class_info = element.split(";")
    # if class_info[0] == "Other":
    #     print (pathway_info)

    first_class = class_info[0].strip()
    second_class = class_info[1].strip()

    return pathway_name, first_class, second_class

def enrichment_analysis_cog(my_list, background_list):
    my_dict = Counter(my_list)
    background_dict = Counter(background_list)
    data = []
    # for category in my_dict:
    for category in my_dict:
        if category == "R" or category == "Y":
            continue
        if category in my_dict:
            a = my_dict[category]
        else:
            a = 0
        b = len(my_list) - a
        if category in background_dict:
            c = background_dict[category]
        else:
            c = 0
        d = len(background_list) - c

        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
        # if p_value >= 0.05:
        #     continue
        # print (category, p_value, oddsratio, a, b, c, d) 
        data.append([category, p_value, oddsratio, a])

    df = pd.DataFrame(data, columns = ["category", "p_value", "fold", "gene_num"])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected

    filtered_df = df[df['p.adj'] < 0.05]
    print (filtered_df)
    print ("differential COG", len(filtered_df))

def enrichment_analysis_kegg(input_ko_ids, background_ko_ids):
    
    input_counts = get_pathways(input_ko_ids)
    background_counts = get_pathways(background_ko_ids)
    print ("pathway count num", len(input_counts), len(background_counts))
    data = []
    # Perform a Fisher's exact test for each pathway
    for pathway_id in set(input_counts.keys()) | set(background_counts.keys()):
        if pathway_id[:2] == "ko":
            continue

        a = input_counts[pathway_id]
        b = len(input_ko_ids) - a
        c = background_counts[pathway_id]
        d = len(background_ko_ids) - c
        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
        
        # if p_value >= 0.05:
        #     continue

        pathway_name, first_class, second_class = get_pathway_name_class(pathway_id)

        # if pathway_id == "map00650":
        #     print (pathway_name, pathway_id, p_value, a, oddsratio)

        data.append([pathway_name, pathway_id, p_value, a, oddsratio, first_class, second_class])

    df = pd.DataFrame(data, columns = ["pathway_name", "pathway_id", "p_value", "gene_num", "fold", "first_class", "second_class"])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected

    for index, row in df.iterrows():
        # diff_pairs[row["genus_pair"]] = row
        if row["pathway_id"] == "map00650":
            print ("<<<<<<<<<<<<\n", row)

    filtered_df = df[df['p.adj'] < 0.05]
    print (filtered_df)
    print ("differential pathway", len(filtered_df))
    return filtered_df

def enrichment_analysis_scfa(input_ko_ids, background_ko_ids, group1, group2):
    input_anno_dict = get_scfa(input_ko_ids)
    background_anno_dict = get_scfa(background_ko_ids)
    # print (input_ko_ids[:5])
    print (input_anno_dict)
    data = []
    for scfa_type in input_anno_dict:
        if scfa_type not in background_anno_dict:
            continue
        a = input_anno_dict[scfa_type]
        b = len(input_ko_ids) - a
        c =  background_anno_dict[scfa_type]
        d = len(background_ko_ids) - c

        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
        print (scfa_type, oddsratio, p_value, a/(a+b), c/(c+d), a, c)
        data.append([scfa_type, oddsratio, p_value, a/(a+b), c/(c+d)])

    df = pd.DataFrame(data, columns = ["scfa", "fold", "p_value", "freq1", "freq2"])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected

    # new_data = []
    # for index, row in df.iterrows():
    #     # print (row)
    #     print (row["scfa"], row["fold"], row["p.adj"])
    #     # if row["scfa"] != "scfa":
    #     new_data.append([row["scfa"], row["freq1"], group1, row["p.adj"]])
    #     new_data.append([row["scfa"], row["freq2"], group2, row["p.adj"]])
    # df = pd.DataFrame(new_data, columns = ["scfa", "freq", "group", "p.adj"])

    df.to_csv("/mnt/d/R_script_files/scfa_freq.csv", sep=',', index=False)



def get_scfa(input_ko_ids):
    anno_dict = defaultdict(int)
    for ko in input_ko_ids:
        if ko[3:] in scfa_dict:
            anno_dict["SCFA"] +=1
            anno_dict[scfa_dict[ko[3:]]] += 1
    return anno_dict


def get_pathways(input_ko_ids):


    # Count the number of input and background KO IDs in each pathway
    input_counts = defaultdict(int)
    for ko_id in input_ko_ids:
        if ko_id[:3] == "ko:":
            ko_id = ko_id[3:]
        for pathway_id in ko_pathway_dict[ko_id]:
            input_counts[pathway_id] += 1
    return input_counts

def get_support_kos(input_ko_ids, pathway_id):
    support_ko = {}
    for ko_id in input_ko_ids:
        if ko_id[:3] == "ko:":
            new_ko_id = ko_id[3:]
        if pathway_id in ko_pathway_dict[new_ko_id]:
            support_ko[ko_id] = 1
    return support_ko

def get_taxa_dict():
    genus2family = {}
    for lineage in taxonomy.taxonomy_dict.values():
        genus = lineage.split(";")[5]
        family = lineage.split(";")[4]
        genus2family[genus] = family
    return genus2family

def get_SCFA_genes():
    scfa_file = "/mnt/d/HGT/seq_ana/scfa/12864_2021_7944_MOESM1_ESM.csv"
    scfa_dict = {}
    f = open(scfa_file)
    f.readline()
    f.readline()

    for line in f:
        array = line.strip().split(",")
        if array[0] == "scfa":
            array[0] = "SCFA"
        elif array[0] == "succinate":
            array[0] = "Succinate"

        scfa_dict[array[1]] = array[0]
    # print (scfa_dict)
    return scfa_dict

if __name__ == "__main__":

    abun_cutoff = 1e-7  #1e-7
    cutoff = 0.1
    level = 5

    bin_size = 1000

    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"
    gff = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.gff"
    phenotype_dict = read_phenotype()
    taxonomy = Taxonomy()
    dat = Data_load()
    dat.read_samples()
    print (dat.cohort_set)
    # tree = Tree("/mnt/d/HGT/time_lines/distribution/bac120_iqtree.nwk")

    f = open('/mnt/d/HGT/seq_ana/ko_pathway_dict.pickle', 'rb')
    ko_pathway_dict = pickle.load(f)
    print ("pathway dict loaded", len(ko_pathway_dict))

    annotation = Annotation(gff)
    annotation.read_gff()

    scfa_dict = get_SCFA_genes()

    shuffle(dat.sample_obj_list)

    focus_cohort = "all"
    # focus_cohort = "YachidaS_2019"

    combination_dict = {}
    data = []
    group_auc = []
    group_list = ["CRC", "control", "T2D",  "IBD"]
    # group_list = ["control", "CRC", "adenoma", "IGT", "T2D", "acute_diarrhoea",  "IBD"]

    # for i in range(len(group_list)):
    #     for j in range(i+1, len(group_list)):
    #         group1 = group_list[i]
    #         group2 = group_list[j]

    group1 = "acute_diarrhoea"
    group2 = "control"
    combination = group1 + " vs. " + group2
    name = group1 + "_" + group2
    # kegg_diff_pair = f"/mnt/d/R_script_files/{name}.csv"
    print ("\n\n<<<<<<<<<<<<<<<<<<<<<<<", combination)

    mar = Marker(group1, group2, dat.sample_obj_list)
    # if not mar.enough_sample_flag:
    #     print ("skip the comparison, because one group is empty.")
    #     continue
    mar.extract_HGT()
    mar.select_diff_HGT(name)
    # mar.classify_breakpoints()
    # all_kegg_diff = f"/mnt/d/R_script_files/{name}_all.csv"
    # mar.compare_all_breakpoints(all_kegg_diff)

        #     break
        # break