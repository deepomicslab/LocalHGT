import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import pickle
 
def read_meta():
    
    sra_sample_dict = {}

    for line in open(meta_data):
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

def read_design():
    
    sample_individual_dict = {}
    sample_time_point = {}

    for line in open(design_file):
        if line.strip() == '':
            continue
        array = line.strip().split()
        if array[0] != 'Sample':
            sample_id = array[0]
            individual = array[3]
            sample_individual_dict[sample_id] = individual
            sample_time_point[sample_id] = int(array[4])
    return sample_individual_dict, sample_time_point

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

        self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        self.abundance = None
        self.split_abundance = None

        self.read = None

# class HGT_event(object):
#     self __init__(self, list):

class Match():

    def __init__(self):
        self.cohort_data = {}
        self.all_hgt = {}
        self.hgt_count = {} # record the hgt exists in how many samples
        self.sample_array_dict = {}
        self.hgt_array_dict = {}
        self.HGT_list = []
        self.paired_bkp = {}
        self.can_match_bkp = {}
        self.can_match_bkp_info = {}
        self.edge_weight_sample = {} # edge weight by the coexist in single sample
        self.correlation_matrix = {} # edge weigth by correlation in all samples
        self.matched_bkp_pairs = set()
        
    def read_samples(self):

        all_acc_file = result_dir + "/acc.list"
        os.system(f"ls {result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            
            td_id = sra_sample_dict[sra_id]
            if td_id[:2] != "CD":
                continue
            # print (sample_individual_dict[td_id])
            # if int(sample_individual_dict[td_id]) not in [1, 2, 3, 4,5, 6, 7, 8, 9, 10]:
            #     continue
            
            my_bkps = self.read_bkp(acc_file)
            self.cohort_data[sra_id] = my_bkps
        self.filter_hgt()
        
    def read_bkp(self, bkp_file):
        my_bkps = []
        f = open(bkp_file)
        all_rows = csv.reader(f)
        sample_dict = {}
        total_HGT_split_num = 0
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
                    break
                eb = Acc_Bkp(row)
                eb.abundance = eb.cross_split_reads/reads_num
                if eb.from_ref_genome == eb.to_ref_genome:
                    continue
                if eb.cross_split_reads < split_cutoff:
                    continue
                if eb.abundance < abun_cutoff:
                    continue
                total_HGT_split_num += eb.cross_split_reads

                # focus = ['GUT_GENOME000147_82', 'GUT_GENOME096083_15']
                # if not (eb.from_ref in focus and eb.to_ref in focus):
                #     continue

                # if eb.hgt_tag not in self.all_hgt:
                sample_dict[eb.hgt_tag] = 1
                my_bkps.append(eb)
        for eb in my_bkps:
            eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
        f.close()
        for hgt_tag in sample_dict:
            if hgt_tag not in self.hgt_count:
                self.hgt_count[hgt_tag] = 0
            self.hgt_count[hgt_tag] += 1
        return my_bkps

    def filter_hgt(self):
        for hgt_tag in self.hgt_count:
            if self.hgt_count[hgt_tag] >= sample_cutoff:
                hgt_index = len(self.all_hgt)
                self.all_hgt[hgt_tag] = hgt_index
        print ("All HGT", len(self.all_hgt))

    def get_correlation(self):
        for hgt_tag in self.all_hgt:
            self.hgt_array_dict[hgt_tag] = [0] * len(self.cohort_data)
        
        sample_index_dict = {}
        index = 0
        for sra_id in self.cohort_data:
            sample_index_dict[sra_id] = index
            index += 1

        for sra_id in self.cohort_data:
            index = sample_index_dict[sra_id]

            for hgt in self.cohort_data[sra_id]:
                if hgt.hgt_tag not in self.all_hgt:
                    continue
                self.hgt_array_dict[hgt.hgt_tag][index] = 1
        
        for hgt_tag_1 in self.hgt_array_dict:

            for hgt_tag_2 in self.hgt_array_dict:
                if hgt_tag_1 == hgt_tag_2:
                    continue
                if hgt_tag_1 not in self.can_match_bkp or hgt_tag_2 not in self.can_match_bkp:
                    continue
                res = stats.spearmanr(self.hgt_array_dict[hgt_tag_1], self.hgt_array_dict[hgt_tag_2])
                # if res.correlation > 0:
                #     correlation_dict[hgt_tag_2] = res.correlation
                bkp_pair = self.get_bkp_pair_name(hgt_tag_1, hgt_tag_2)
                self.correlation_matrix[bkp_pair] = res.correlation

    def get_precise_HGT(self):

        data = []
        for sra_id in self.cohort_data:
            sample_dict = {}
            for hgt in self.cohort_data[sra_id]:
                if hgt.hgt_tag in self.paired_bkp:
                    # print (sra_id, hgt.from_ref, hgt.from_bkp, hgt.to_ref, hgt.to_bkp)
                    sample_dict[hgt.hgt_tag] = hgt
            # for hgt_tag in self.paired_bkp:
            #     hgt_tag_2 = self.paired_bkp[hgt_tag]
            for bkp_pair in self.matched_bkp_pairs:
                hgt_tag = bkp_pair[0]
                hgt_tag_2 = bkp_pair[1]
                if hgt_tag in sample_dict and hgt_tag_2 in sample_dict:
                    # print (sample_dict[hgt_tag].from_ref, sample_dict[hgt_tag].from_bkp, sample_dict[hgt_tag].to_ref, sample_dict[hgt_tag].to_bkp, \
                    # sample_dict[hgt_tag_2].from_ref, sample_dict[hgt_tag_2].from_bkp, sample_dict[hgt_tag_2].to_ref, sample_dict[hgt_tag_2].to_bkp)

                    array_1 = [sample_dict[hgt_tag].from_ref, sample_dict[hgt_tag].from_bkp, sample_dict[hgt_tag].to_ref, sample_dict[hgt_tag].to_bkp]
                    array_2 = [sample_dict[hgt_tag_2].from_ref, sample_dict[hgt_tag_2].from_bkp, sample_dict[hgt_tag_2].to_ref, sample_dict[hgt_tag_2].to_bkp]
                    
                    if array_1[0] != array_2[0]:
                        array_2 = [array_2[2], array_2[3], array_2[0], array_2[1]]

                    if abs(array_1[1] - array_2[1]) < abs(array_1[3] - array_2[3]):
                        insertion = [array_1[0], array_1[1]]
                        deletion = [array_1[2], array_1[3], array_2[3]]
                    else:
                        insertion = [array_1[2], array_1[3]]
                        deletion = [array_1[0], array_1[1], array_2[1]]
                    if deletion[2] < deletion[1]:
                        a = deletion[1]
                        deletion[1] = deletion[2]
                        deletion[2] = a
                    # print (sra_id, insertion, deletion)
                    data.append([sra_id] + insertion + deletion)
        df = pd.DataFrame(data, columns = ["sample", "receptor", "insert_locus", "donor", "delete_start", "delete_end"])
        df.to_csv(identified_hgt, sep=',')

    def matching(self):
        # Create a graph with nodes on the left and right
        G = nx.Graph()
        bkp_list = list(self.can_match_bkp)
        possible_match_num = 0
        print ("bkp num", len(bkp_list))
        for i in range(len(bkp_list)):
            for j in range(i+1, len(bkp_list)):
        # for bkp1 in bkp_list:    
        #     for bkp2 in bkp_list:   
        #         if bkp1 == bkp2:
        #             continue
                bkp1, bkp2 = bkp_list[i], bkp_list[j]
                if self.check_point_share(bkp1, bkp2):
                    possible_match_num += 1
                    bkp_pair = self.get_bkp_pair_name(bkp1, bkp2)
                    edge_weight_by_sample = 0
                    edge_weight_by_correlation = 0
                    if bkp_pair in self.edge_weight_sample and self.edge_weight_sample[bkp_pair] > 0:
                        edge_weight_by_sample = 1
                    if self.correlation_matrix[bkp_pair] > 0.7:
                        edge_weight_by_correlation = self.correlation_matrix[bkp_pair]
                    gene_length = self.cal_gene_length(bkp1, bkp2) # /bin
                    edge_weight = (edge_weight_by_correlation + edge_weight_by_sample)/gene_length/bin_size
                    # edge_weight = (edge_weight_by_sample)/gene_length/bin_size
                    # print (bkp_pair.split("&"), edge_weight)
                    if edge_weight > 0:
                        G.add_edge(bkp1, bkp2, weight = edge_weight)
                # else:
                #     G.add_edge(bkp1, bkp2, weight=0)

        G.add_nodes_from(bkp_list)
        for cc in nx.connected_components(G):
            G_lcc = G.subgraph(cc)
            # print("Nodes in largest connected component:", len(G_lcc.nodes))
            matching = nx.algorithms.matching.max_weight_matching(G_lcc, weight='weight')
            # print("matching")
            self.matched_bkp_pairs = self.matched_bkp_pairs | matching

        print("N.O. of matched bkp pairs:", len(self.matched_bkp_pairs), "Possible match:", possible_match_num/2)
        for bkp_pair in self.matched_bkp_pairs:
            self.paired_bkp[bkp_pair[0]] = bkp_pair[1]
            self.paired_bkp[bkp_pair[1]] = bkp_pair[0]

    def check_point_share(self, hgt_tag_1, hgt_tag_2):
        array_1 = hgt_tag_1.split("&")
        array_2 = hgt_tag_2.split("&")

        if not (array_1[0] == array_2[0] and array_1[2] == array_2[2]) and not (array_1[0] == array_2[2] and array_1[2] == array_2[0]):
            return False
        if array_1[0] != array_2[0]:
            array_2 = [array_2[2], array_2[3], array_2[0], array_2[1]]
        if (array_1[1] == array_2[1] and abs(int(array_1[3]) - int(array_2[3])) > 3) or (array_1[3] == array_2[3] and abs(int(array_1[1]) - int(array_2[1])) > 3):
            return True
        return False

    def cal_gene_length(self, hgt_tag_1, hgt_tag_2):
        array_1 = hgt_tag_1.split("&")
        array_2 = hgt_tag_2.split("&")

        if array_1[0] != array_2[0]:
            array_2 = [array_2[2], array_2[3], array_2[0], array_2[1]]
        if array_1[1] == array_2[1]:
            return abs(int(array_1[3]) - int(array_2[3]))
        if array_1[3] == array_2[3] :
            return abs(int(array_1[1]) - int(array_2[1]))

    def check_match_in_sample(self):
        for sra_id in self.cohort_data:
            ref_pair_dict = {}
            genome_pair_dict = {}
            ref_genome_relation = {}
            for hgt in self.cohort_data[sra_id]:
                if hgt.hgt_tag not in self.all_hgt:
                    continue
                ref_pair = "&".join(sorted([hgt.from_ref, hgt.to_ref]))
                genome_pair = "&".join(sorted([hgt.from_ref_genome, hgt.to_ref_genome]))
                ref_genome_relation[ref_pair] = genome_pair

                if ref_pair not in ref_pair_dict:
                    ref_pair_dict[ref_pair] = []
                ref_pair_dict[ref_pair].append(hgt.hgt_tag)

                if genome_pair not in genome_pair_dict:
                    genome_pair_dict[genome_pair] = []
                genome_pair_dict[genome_pair].append(hgt.hgt_tag)

            for ref_pair in ref_pair_dict:
                # if len(genome_pair_dict[ref_genome_relation[ref_pair]]) == 2 and len(ref_pair_dict[ref_pair]) == 2 and self.check_point_share(ref_pair_dict[ref_pair][0], ref_pair_dict[ref_pair][1]):
                if len(ref_pair_dict[ref_pair]) == 2 and self.check_point_share(ref_pair_dict[ref_pair][0], ref_pair_dict[ref_pair][1]):
                   # print (sra_id, ref_pair, ref_pair_dict[ref_pair])
                    bkp_pair = self.get_bkp_pair_name(ref_pair_dict[ref_pair][0], ref_pair_dict[ref_pair][1])
                    if bkp_pair not in self.edge_weight_sample:
                        self.edge_weight_sample[bkp_pair] = 0
                    self.edge_weight_sample[bkp_pair] += 1
    
    def get_bkp_pair_name(self, bkp1, bkp2):
        bkp_pair = "%".join(sorted([bkp1, bkp2]))
        return bkp_pair

    def classify_bkp(self):
        bkp_list = list(self.all_hgt.keys())
        # for bkp1 in bkp_list:    
        #     for bkp2 in bkp_list:   
        #         if bkp1 == bkp2:
        #             continue
        for i in range(len(bkp_list)):
            for j in range(i+1, len(bkp_list)):
                bkp1, bkp2 = bkp_list[i], bkp_list[j]

                if self.check_point_share(bkp1, bkp2):
                    if bkp1 not in self.can_match_bkp:
                        self.can_match_bkp[bkp1] = 0
                        # self.can_match_bkp_info[bkp1] = []
                    self.can_match_bkp[bkp1] += 1
                    # self.can_match_bkp_info[bkp1] += [bkp2]

                    if bkp2 not in self.can_match_bkp:
                        self.can_match_bkp[bkp2] = 0
                        # self.can_match_bkp_info[bkp2] = []
                    self.can_match_bkp[bkp2] += 1
                    # self.can_match_bkp_info[bkp2] += [bkp1]
        print ("%s bkp after filtering, and %s bkp have possible matched object."%(len(bkp_list), len(self.can_match_bkp)))
        single_match_num = 0
        alt_match_num = 0
        for bkp in self.can_match_bkp:
            if self.can_match_bkp[bkp] > 1:
                alt_match_num += 1
                # print ("multiple pairs", bkp, self.can_match_bkp_info[bkp])
            else:
                single_match_num += 1
        with open(saved_can_match_bkp, 'wb') as f:
            pickle.dump(self.can_match_bkp, f)
        print ("%s has one matched, %s has more than one."%(single_match_num, alt_match_num))
    
    def load_can_match(self):
        # open the file in read binary mode
        with open(saved_can_match_bkp, 'rb') as f:
            # load the dictionary object from the file using pickle
            self.can_match_bkp = pickle.load(f)       

    def draw(self, G):
        # # Compute the bipartite layout
        pos = nx.bipartite_layout(G, linked_bkp)

        # Draw the bipartite graph
        nx.draw(G, pos=pos, with_labels=True, node_color='Blue', edge_color='gray')
        plt.savefig('/mnt/d/HGT/time_lines/bipartite_graph.pdf')

if __name__ == "__main__":

    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    design_file = "/mnt/d/HGT/time_lines/sample_design.tsv"
    result_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    # hgt_table = "/mnt/d/HGT/time_lines/SRP366030.HGT.table.csv"
    identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    saved_can_match_bkp = "/mnt/d/HGT/time_lines/SRP366030.can_match.pickle"

    bin_size = 100
    split_cutoff = 0  #10
    sample_cutoff = 0  # 8
    abun_cutoff = 1e-7  #1e-7


    # for bin_size in [100, 200, 500, 1000]:
    #     for split_cutoff in range(6, 20, 2):
    #         for sample_cutoff in range(2, 11, 2):
    sra_sample_dict = read_meta()
    sample_individual_dict, sample_time_point = read_design()
    # print (sample_individual_dict)
    tim = Match()
    tim.read_samples()
    print ("load is done.")
    tim.classify_bkp()
    # tim.load_can_match()
    tim.check_match_in_sample()
    tim.get_correlation()
    tim.matching()
    tim.get_precise_HGT()


    # print (bin_size, split_cutoff, sample_cutoff)
