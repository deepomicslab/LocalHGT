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
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
import powerlaw

from mechanism_taxonomy import Taxonomy


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
        if self.disease == "acute_diarrhoea":
            self.disease = "diarrhoea"
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

    def judge_scale_free(self, level, edge_num):
        HGT_matrix, total_edge_num =  self.get_HGT_matrix(level, edge_num)
        HGT_matrix = nx.from_numpy_matrix(HGT_matrix)
        p1, p2, p3 = infer_scale_free(HGT_matrix)
        return p1, p2, p3, total_edge_num

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

    def get_cohort_num(self):
        cohort_dict = {}
        for sample in self.sample_obj_list:
            if sample.cohort not in cohort_dict:
                cohort_dict[sample.cohort] = 0
            cohort_dict[sample.cohort] += 1
        print (cohort_dict, "total num", sum(list(cohort_dict.values())))
        
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

class Network():

    def __init__(self, all_data):
        self.data = all_data

    def compare_network(self): 
        group1 = "CRC"
        group2 = "IBD"
        properties = ['density', 'transitivity', 'algebraic_connectivity', 'assortativity', 'Node', "Edge"]
        data = []
        num_count = [{}, {}]
        edge_num_list = [10, 12, 20, 30, 40, 50]
        for level in range(1, 7):
            edge_num = edge_num_list[level-1]
            cohort_sam_num = {}
            cohort_base = {}
            properties_dict = {}
            for pro in properties:
                properties_dict[pro] = [[],[]]
            for sample in self.data:

                if sample.disease ==  group1 or group1 in sample.full_disease:
                    index = 0
                    num_count[0][sample.ID] = 1
                elif sample.disease == group2 or group2 in sample.full_disease:
                    index = 1
                    num_count[1][sample.ID] = 1
                else:
                    continue

                pro_list, total_edge_num = sample.get_network_properties(level, edge_num)
                if total_edge_num < edge_num:
                    continue
                for i in range(len(properties)):
                    value = pro_list[i]#/sample.bases
                    origin = pro_list[i]
                    properties_dict[properties[i]][index].append(value)
                    data.append([properties[i], value, sample.disease, sample.cohort, level_list[level-1], origin])
            for i in range(len(properties)):
                # U1, p = mannwhitneyu(properties_dict[properties[i]][0], properties_dict[properties[i]][1])
                U1, p = scipy.stats.ranksums(properties_dict[properties[i]][0], properties_dict[properties[i]][1])
                print (len(properties_dict[properties[i]][0]), len(properties_dict[properties[i]][1]), level_list[level-1], properties[i], p)#, np.mean(properties_dict[properties[i]][0]), np.mean(properties_dict[properties[i]][1]), sep = "\t")
        df = pd.DataFrame(data, columns = ["Property", "Value", "Group", "Cohort", "Level", "Origin"])
        df.to_csv('/mnt/d/R_script_files/network_comparison_normalized.csv', sep=',')
        print ("%s num is %s, %s num is %s."%(group1, len(num_count[0]), group2, len(num_count[1])))

    def compare_network_mul_group(self):        
        # properties = ['density', 'transitivity', 'algebraic_connectivity', 'assortativity', 'Node', "Edge"]
        properties = ['assortativity', 'transitivity', 'algebraic_connectivity', 'density' ]
        # properties = ['transitivity']
        data = []
        
        num_count = [{}, {}, {}, {}]
        edge_num_list = [10, 12, 20, 30, 40, 50]
        for level in range(1, 2):
            group_dict = {}
            group_pro_dict = {}
            edge_num = edge_num_list[level-1]
            for sample in self.data:

                if len(sample.full_disease) != 1:
                    continue
                if sample.disease == '':
                    continue
                if sample.disease == "control" and sample.full_disease[0] != "healthy":
                    continue

                pro_list, total_edge_num = sample.get_network_properties(level, edge_num)
                if total_edge_num < edge_num:
                    continue
                if sample.disease not in group_dict:
                    group_dict[sample.disease] = 0
                    group_pro_dict[sample.disease] = {}
                group_dict[sample.disease] += 1
                for i in range(len(properties)):
                    value = pro_list[i]#/sample.bases
                    origin = pro_list[i]
                    data.append([properties[i], value, sample.disease, sample.cohort, level_list[level-1], origin])

                    if properties[i] not in group_pro_dict[sample.disease]:
                        group_pro_dict[sample.disease][properties[i]] = []
                    group_pro_dict[sample.disease][properties[i]].append(value)

            print (level, group_dict)
            group_list = list(group_pro_dict.keys())
            for proper in group_pro_dict["control"]:
                matrix_dat = []
                p_values = []
                for i in range(len(group_list)):
                    row = [group_list[i]]
                    for j in range(0, len(group_list)):
                        group1 = group_list[i]
                        group2 = group_list[j]
                        # U1, p = mannwhitneyu(group_pro_dict[group][proper], group_pro_dict["control"][proper])
                        U1, p = scipy.stats.ranksums(group_pro_dict[group1][proper], group_pro_dict[group2][proper])
                        # p = '{:.1e}'.format(p)
                        if i >= j:
                            row.append(p)
                            p_values.append(float(p))
                        else:
                            row.append("NA")
                        # if i > j:
                        #     matrix_dat.append([proper, group1, group2, p])
                        # else:
                        #     matrix_dat.append([proper, group1, group2, "NA"])
                        print (proper, group1, group2, p)
                    matrix_dat.append(row)

                # p_values is a list/array of unadjusted p-values
                print (p_values)
                rejected, adjusted_pvalues, _, _ = multipletests(p_values, method='fdr_bh')
                # rejected, adjusted_pvalues, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
                print (adjusted_pvalues)
                z = 0
                new_matrix_dat = []
                for row in matrix_dat:
                    for e in range(1, len(row)):
                        if row[e] != "NA":
                            new_p = adjusted_pvalues[z]
                            row[e] = '{:.1e}'.format(new_p)
                            z += 1
                    new_matrix_dat.append(row)
                matrix_dat = new_matrix_dat

                df = pd.DataFrame(matrix_dat, columns = ["group"] + group_list)
                df.to_csv('/mnt/d/R_script_files/network_comparison_matrix_%s.csv'%(proper), sep=',')
                        
                print ("<<<<<<<<<<<<<<<<")
        df = pd.DataFrame(data, columns = ["Property", "Value", "Group", "Cohort", "Level", "Origin"])
        df.to_csv('/mnt/d/R_script_files/network_comparison_normalized.csv', sep=',')

    def infer_sale(self):
        # detect scale free
        data = []
        new_data = []
        # edge_num_list = [10, 10, 30, 40, 100, 120]
        f = open("/mnt/d/R_script_files/scale_free_count.txt", 'w')
        edge_num_list = [10, 12, 20, 30, 40, 50]
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
            print (level, level_list[level-1], scale_free_num, network_num, scale_free_num/network_num)
            new_data.append([level, level_list[level-1], scale_free_num, network_num, scale_free_num/network_num])

        f.close()
        # df = pd.DataFrame(data, columns = ["ratio", "Comparison", "Group", "Cohort", "Level"])
        # df.to_csv('/mnt/c/Users/user/Desktop/HGT/HGT_R_plot_files//scale_free.csv', sep=',')
        df = pd.DataFrame(new_data, columns = ["level_index", "level", "scale_free_num", "network_num", "Frequency"])
        df.to_csv('/mnt/d/R_script_files//scale_free.csv', sep=',')

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

class FR_matrix(Data_load):

    def __init__(self, group):
        self.group = group
        self.all_bkps = {}
        self.level = 6
        self.node_index = {}
        self.read_samples()
        self.sample_num = len(self.all_bkps)
        self.matrix = None
        self.matrix_sample_dict ={} 
        self.nodes_list = []

    def get_node_index(self):
        for sra_id in self.all_bkps:
            for bkp in self.all_bkps[sra_id]:
                from_taxon = bkp.from_ref_lineage.split(";")[self.level]
                to_taxon = bkp.to_ref_lineage.split(";")[self.level]
                for taxa in [from_taxon, to_taxon]:
                    if taxa not in self.node_index:
                        self.node_index[taxa] = len(self.node_index)
        print ("node number", len(self.node_index))
        self.get_node_list()
    
    def get_matrix(self):
        #= np.zeros((len(self.node_index), len(self.node_index)))    
        for sra_id in self.all_bkps:
            for bkp in self.all_bkps[sra_id]:
                from_taxon = bkp.from_ref_lineage.split(";")[self.level]
                to_taxon = bkp.to_ref_lineage.split(";")[self.level]  

                # from_index =  self.node_index[from_taxon]
                # to_index =  self.node_index[to_taxon]
                if from_taxon not in self.matrix_sample_dict:
                    self.matrix_sample_dict[from_taxon] = defaultdict(set)
                if to_taxon not in self.matrix_sample_dict:
                    self.matrix_sample_dict[to_taxon] = defaultdict(set)

                self.matrix_sample_dict[from_taxon][to_taxon].add(sra_id)
                self.matrix_sample_dict[to_taxon][from_taxon].add(sra_id)
        
        self.matrix = np.zeros((len(self.node_index), len(self.node_index)))  
        for taxon_1 in self.matrix_sample_dict:
            for taxon_2 in self.matrix_sample_dict[taxon_1]:
                support_sample_num = len(self.matrix_sample_dict[taxon_1][taxon_2])
                support_freq = round(support_sample_num/self.sample_num, 3)

                taxon_1_index = self.node_index[taxon_1]
                taxon_2_index = self.node_index[taxon_2]

                self.matrix[taxon_1_index][taxon_2_index] = support_freq
                self.matrix[taxon_2_index][taxon_1_index] = support_freq

        df = pd.DataFrame(self.matrix, index=self.nodes_list, columns=self.nodes_list)
        df.to_csv("/mnt/d/HGT/seq_ana/species_matrix_%s.csv"%(self.group))
        
    def get_node_list(self):
        self.nodes_list = [""] * len(self.node_index)
        for node in self.node_index:
            self.nodes_list[self.node_index[node]] = node

    def read_samples(self):

        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        # os.system(f"ls {tgs_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        # os.system(f"ls {wenkui_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            sample = Sample(my_bkps, sra_id, phenotype_dict[sra_id])
            if self.group != "all":
                if sample.disease == "control" and sample.full_disease[0] != "healthy":
                    continue
                if sample.disease != self.group:
                    continue
            if len(my_bkps) > 0 and sra_id in phenotype_dict:
                self.all_bkps[sra_id] = my_bkps

        print ("data is loaded, number of sample is %s in %s"%(len(self.all_bkps), self.group))

    def get_genus_phylum(self):
        genus_phylum_pair = {}
        taxa_count = defaultdict(int)
        for sra_id in self.all_bkps:
            sample_dict = defaultdict(int)
            for bkp in self.all_bkps[sra_id]:

                from_taxon = bkp.from_ref_lineage.split(";")[self.level]
                if len(from_taxon) == 3:
                    from_phylum = "p__"
                else:
                    from_phylum = bkp.from_ref_lineage.split(";")[1]

                to_taxon = bkp.to_ref_lineage.split(";")[self.level] 
                if len(to_taxon) == 3:
                    to_phylum = "p__"
                else:
                    to_phylum = bkp.to_ref_lineage.split(";")[1]

                genus_phylum_pair[from_taxon] = from_phylum
                genus_phylum_pair[to_taxon] = to_phylum

                sample_dict[from_taxon] += 1
                sample_dict[to_taxon] += 1
            for taxon in sample_dict:
                taxa_count[taxon] += 1
        
        sorted_taxa_count = sorted(taxa_count.items(), key=lambda item: item[1], reverse = True)
        top_taxa = dict(sorted_taxa_count[:600])

        dark2=['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D']
        set1=['#E41A1C', '#377EB8','#FFFF33', '#984EA3', '#FF7F00',  '#A65628', '#F781BF', '#999999', '#4DAF4A']
        color_palette = dark2 + set1
        index_dict = {'Firmicutes_A': 0, 'Bacteroidota': 1, 'Firmicutes': 2, 'Proteobacteria': 3, 'Actinobacteriota': 4, 'Firmicutes_C': 5, \
        'Verrucomicrobiota': 6, 'Firmicutes_B': 7, "Cyanobacteria":8, "Fusobacteriota":9}      
        data = []
        phylum_set = set()
        for taxon in genus_phylum_pair:
            if taxon not in top_taxa or len(taxon) == 3:
                continue
            phylum = genus_phylum_pair[taxon]
            if phylum[3:] not in index_dict:
                phylum = "other"
            phylum_set.add(phylum)
            data.append([taxon, phylum])

        df = pd.DataFrame(data, columns=[level_list[self.level-1], "phylum"])
        df.to_csv("/mnt/d/HGT/seq_ana/%s_phylum_pair.csv"%(level_list[self.level-1]))
        print ("phylum count", len(phylum_set))
        print (phylum_set)


if __name__ == "__main__":

    abun_cutoff = 1e-7  #1e-7
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"
    phenotype_dict = read_phenotype()
    taxonomy = Taxonomy()

    dat = Data_load()
    dat.read_samples()
    dat.get_cohort_num()
    net = Network(dat.sample_obj_list)
    net.compare_network_mul_group()
    # net.infer_sale()

    # group_list = ["control", "CRC", "adenoma", "IGT", "T2D", "acute_diarrhoea",  "IBD"]
    # for group in group_list:
    #     fr = FR_matrix(group)
    #     fr.get_node_index()
    #     fr.get_matrix()

    # fr = FR_matrix("all")
    # fr.get_genus_phylum()