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
        # group1 = "control"
        # group2 = "IBD"
        # group3 = "CRC"
        # group4 = "T2D"
        
        properties = ['density', 'transitivity', 'algebraic_connectivity', 'assortativity', 'Node', "Edge"]
        data = []
        num_count = [{}, {}, {}, {}]
        edge_num_list = [10, 12, 20, 30, 40, 50]
        for level in range(1, 4):
            group_dict = {}
            edge_num = edge_num_list[level-1]
            for sample in self.data:

                if len(sample.full_disease) != 1:
                    continue
                if sample.disease == '':
                    continue

                pro_list, total_edge_num = sample.get_network_properties(level, edge_num)
                if total_edge_num < edge_num:
                    continue
                if sample.disease not in group_dict:
                    group_dict[sample.disease] = 0
                group_dict[sample.disease] += 1
                for i in range(len(properties)):
                    value = pro_list[i]#/sample.bases
                    origin = pro_list[i]
                    data.append([properties[i], value, sample.disease, sample.cohort, level_list[level-1], origin])
            print (level, group_dict)
        df = pd.DataFrame(data, columns = ["Property", "Value", "Group", "Cohort", "Level", "Origin"])
        df.to_csv('/mnt/d/R_script_files/network_comparison_normalized.csv', sep=',')
        


    def infer_sale(self):
        # detect scale free
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

if __name__ == "__main__":
    abun_cutoff = 1e-7  #1e-7
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    phenotype_dict = read_phenotype()
    taxonomy = Taxonomy()
    dat = Data_load()
    dat.read_samples()
    net = Network(dat.sample_obj_list)
    net.compare_network_mul_group()