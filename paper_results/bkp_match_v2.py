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
from pyfaidx import Fasta
from sklearn.cluster import DBSCAN
import random
 
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
        self.bk_1_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size))
        self.bk_2_tag = self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        self.abundance = None
        self.split_abundance = None

        self.read = None

# class HGT_event(object):
#     self __init__(self, list):

class Match():

    def __init__(self):
        self.cohort_data = {}
        self.sample_array_dict = {}
        self.hgt_array_dict = {}
        self.HGT_list = []
        self.paired_bkp = {}
        self.can_match_bkp = {}
        self.can_match_bkp_info = {}
        self.edge_weight_sample = {} # edge weight by the coexist in single sample
        self.correlation_matrix = {} # edge weigth by correlation in all samples
        self.matched_bkp_pairs = set()
        self.ref_fasta = Fasta(database)
        
    def read_samples(self):

        all_acc_file = result_dir + "/acc.list"
        os.system(f"ls {result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]   

            # if sra_id != "SRR18490938":
            #     continue

            # td_id = sra_sample_dict[sra_id]
            # if td_id[:2] != "CD":
            #     continue
            # print (sample_individual_dict[td_id])
            
            my_bkps = self.read_bkp(acc_file)
            self.cohort_data[sra_id] = my_bkps
        
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
                # if eb.cross_split_reads < split_cutoff:
                #     continue
                # if eb.abundance < abun_cutoff:
                #     continue
                total_HGT_split_num += eb.cross_split_reads

                # focus = ['GUT_GENOME000147_82', 'GUT_GENOME096083_15']
                # if not (eb.from_ref in focus and eb.to_ref in focus):
                #     continue
                sample_dict[eb.hgt_tag] = 1
                my_bkps.append(eb)
        for eb in my_bkps:
            eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
        f.close()
        return my_bkps

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

    def delete_direction(self, pos_list):
        ## pos_list: [pos1, direction1, strand1, pos2, dir2, strand2]
        if pos_list[0] > pos_list[3]:
            pos_list = pos_list[3:] + pos_list[:3]
        delete_start =  pos_list[0]
        delete_end = pos_list[3]
        dir_flag = True
        if pos_list[1] != "tail" or pos_list[4] != "head":
            dir_flag = False
        return delete_start, delete_end, dir_flag

    def check_if_match(self, bkp_obj_1, bkp_obj_2, ID):
        flag = False
        if not ((bkp_obj_1.from_ref == bkp_obj_2.from_ref and bkp_obj_1.to_ref == bkp_obj_2.to_ref) or (bkp_obj_1.to_ref == bkp_obj_2.from_ref and bkp_obj_1.from_ref == bkp_obj_2.to_ref)):
            return False
        if bkp_obj_1.from_ref == bkp_obj_2.from_ref and abs(bkp_obj_1.from_bkp - bkp_obj_2.from_bkp) < max_diff:
            if bkp_obj_1.to_ref == bkp_obj_2.to_ref and abs(bkp_obj_1.to_bkp - bkp_obj_2.to_bkp) > max_diff:
                receptor = bkp_obj_1.from_ref
                insert_pos = bkp_obj_1.from_bkp
                donor = bkp_obj_1.to_ref
                delete_start, delete_end, dir_flag = self.delete_direction([bkp_obj_1.to_bkp, bkp_obj_1.to_side, bkp_obj_1.to_strand, bkp_obj_2.to_bkp, bkp_obj_2.to_side, bkp_obj_2.to_strand])
                flag = True and dir_flag

        elif bkp_obj_1.to_ref == bkp_obj_2.from_ref and abs(bkp_obj_1.to_bkp - bkp_obj_2.from_bkp) < max_diff:
            if bkp_obj_1.from_ref == bkp_obj_2.to_ref and abs(bkp_obj_1.from_bkp - bkp_obj_2.to_bkp) > max_diff:
                receptor = bkp_obj_1.to_ref
                insert_pos = bkp_obj_1.to_bkp
                donor = bkp_obj_1.from_ref
                delete_start, delete_end, dir_flag = self.delete_direction([bkp_obj_1.from_bkp, bkp_obj_1.from_side, bkp_obj_1.from_strand, bkp_obj_2.to_bkp, bkp_obj_2.to_side, bkp_obj_2.to_strand])
                flag = True and dir_flag

        elif bkp_obj_1.from_ref == bkp_obj_2.to_ref and abs(bkp_obj_1.from_bkp - bkp_obj_2.to_bkp) < max_diff:
            if bkp_obj_1.to_ref == bkp_obj_2.from_ref and abs(bkp_obj_1.to_bkp - bkp_obj_2.from_bkp) > max_diff:
                receptor = bkp_obj_1.from_ref
                insert_pos = bkp_obj_1.from_bkp
                donor = bkp_obj_1.to_ref
                delete_start, delete_end, dir_flag = self.delete_direction([bkp_obj_1.to_bkp, bkp_obj_1.to_side, bkp_obj_1.to_strand, bkp_obj_2.from_bkp, bkp_obj_2.from_side, bkp_obj_2.from_strand])
                flag = True and dir_flag

        elif bkp_obj_1.to_ref == bkp_obj_2.to_ref and abs(bkp_obj_1.to_bkp - bkp_obj_2.to_bkp) < max_diff:
            if bkp_obj_1.from_ref == bkp_obj_2.from_ref and abs(bkp_obj_1.from_bkp - bkp_obj_2.from_bkp) > max_diff:
                receptor = bkp_obj_1.to_ref
                insert_pos = bkp_obj_1.to_bkp
                donor = bkp_obj_1.from_ref
                delete_start, delete_end, dir_flag = self.delete_direction([bkp_obj_1.from_bkp, bkp_obj_1.from_side, bkp_obj_1.from_strand, bkp_obj_2.from_bkp, bkp_obj_2.from_side, bkp_obj_2.from_strand])
                flag = True and dir_flag

        flag = flag and bkp_obj_1.if_reverse == bkp_obj_2.if_reverse
        # if flag:
        #     if len(self.pair_dict[donor+ "&" + str(int(delete_start/bin_size))]) >1 and len(self.pair_dict[donor+ "&" + str(int(delete_end/bin_size))]) > 1:
        #         flag = False
        if flag: 
            if delete_end - delete_start < min_hgt_len:
                flag = False
        if flag:
            event_data = [ID, receptor,insert_pos,donor,delete_start,delete_end, bkp_obj_1.if_reverse]

            # match_num = self.remove_ambiguity(event_data)
            match_num = self.remove_ambiguity_pop(event_data)
            # match_num = 2
            print (ID, receptor,insert_pos,donor,delete_start,delete_end, match_num)
            if match_num <= 2:
                result_data.append(event_data)
            else:
                flag = False
        return flag

    def remove_ambiguity(self, event_data):
        ## ID, receptor,insert_pos,donor,delete_start,delete_end
        ID = event_data[0]
        pos = []
        for bkp in self.cohort_data[ID]:
            if bkp.from_ref == event_data[1] and abs(bkp.from_bkp - event_data[2]) < max_diff: # ins 
                if bkp.to_ref == event_data[3]:
                # if bkp.to_ref.split("_")[1] == event_data[3].split("_")[1]:
                    pos.append(bkp.to_bkp)
            if bkp.to_ref == event_data[1] and abs(bkp.to_bkp - event_data[2]) < max_diff: # ins
                if bkp.from_ref == event_data[3]:
                # if bkp.from_ref.split("_")[1] == event_data[3].split("_")[1]:
                    pos.append(bkp.from_bkp)
        # if len(pos) > 2:
        # print (pos)
        return len(pos)

    def remove_ambiguity_pop(self, event_data):
        ## ID, receptor,insert_pos,donor,delete_start,delete_end
        ID = event_data[0]
        pos = []
        sample_list = list(self.cohort_data.keys()) 
        
        random.shuffle(sample_list)
        select_samples = sample_list[:200] + [ID]
        select_samples = list(set(select_samples)) # get unique samples

        for ID in select_samples:
            for bkp in self.cohort_data[ID]:
                if bkp.from_ref == event_data[1] and abs(bkp.from_bkp - event_data[2]) < max_diff: # ins 
                    if bkp.to_ref == event_data[3]:
                    # if bkp.to_ref.split("_")[1] == event_data[3].split("_")[1]:
                        pos.append(bkp.to_bkp)
                elif bkp.to_ref == event_data[1] and abs(bkp.to_bkp - event_data[2]) < max_diff: # ins
                    if bkp.from_ref == event_data[3]:
                    # if bkp.from_ref.split("_")[1] == event_data[3].split("_")[1]:
                        pos.append(bkp.from_bkp)
                else:
                    continue
        # if len(pos) > 2:
        # print (pos)
        if len(pos) == 0:
            return 0
        dbscan = DBSCAN(eps=bin_size, min_samples=1)
        dbscan.fit(np.array(pos).reshape(-1, 1))
        cluster_num = max(dbscan.labels_) + 1
        # print the cluster labels
        # print(dbscan.labels_)
        return cluster_num

    def cal_gene_length(self, hgt_tag_1, hgt_tag_2):
        array_1 = hgt_tag_1.split("&")
        array_2 = hgt_tag_2.split("&")

        if array_1[0] != array_2[0]:
            array_2 = [array_2[2], array_2[3], array_2[0], array_2[1]]
        if array_1[1] == array_2[1]:
            return abs(int(array_1[3]) - int(array_2[3]))
        if array_1[3] == array_2[3] :
            return abs(int(array_1[1]) - int(array_2[1]))
    
    def draw(self, G):
        # # Compute the bipartite layout
        pos = nx.bipartite_layout(G, linked_bkp)

        # Draw the bipartite graph
        nx.draw(G, pos=pos, with_labels=True, node_color='Blue', edge_color='gray')
        plt.savefig('/mnt/d/HGT/time_lines/bipartite_graph.pdf')

    def check_if_bkp_at_ends(self, bkp):
        # check if  the breakpoints locate at genome ends.
        valid_flag = True #

        from_ref_len = len(self.ref_fasta[bkp.from_ref])
        if bkp.from_bkp < window/2 or from_ref_len - bkp.from_bkp <window/2:
            valid_flag = False

        to_ref_len = len(self.ref_fasta[bkp.to_ref])
        if bkp.to_bkp < window/2 or to_ref_len - bkp.to_bkp <window/2:
            valid_flag = False

        return valid_flag

    def match_each_sample(self, ID):
        possible_match_num = 0
        sample_bkps = self.cohort_data[ID]
        
        valid_num = 0
        for bkp in sample_bkps:
            # if bkp.abundance >= abun_cutoff:
            #     valid_num += 1
            if bkp.cross_split_reads > 1:
                valid_num += 1
        print ("%s: raw bkp num is %s, filtered bkp num is %s."%(ID, len(sample_bkps), valid_num))

        for i in range(len(sample_bkps)):
            if not self.check_if_bkp_at_ends(sample_bkps[i]): # the bkp should not in reference ends
                continue
            # if sample_bkps[i].abundance < abun_cutoff:
            #     continue
            if sample_bkps[i].cross_split_reads < 2:
                continue
            for j in range(i+1, len(sample_bkps)):
                if not self.check_if_bkp_at_ends(sample_bkps[j]): # the bkp should not in reference ends
                    continue
                # if sample_bkps[j].abundance < abun_cutoff:
                #     continue
                if sample_bkps[j].cross_split_reads < 2:
                    continue
                bkp1 = sample_bkps[i].hgt_tag
                bkp2 = sample_bkps[j].hgt_tag
                flag = self.check_if_match(sample_bkps[i], sample_bkps[j], ID)
                if flag:
                    # print (bkp1, bkp2)
                    possible_match_num += 1
        print ("%s has %s bkps and its match num is %s"%(ID, len(sample_bkps), possible_match_num))



if __name__ == "__main__":

    bin_size = 100
    window = 200
    split_cutoff = 0  #10
    sample_cutoff = 0  # 8
    abun_cutoff = 5e-8   #1e-7 #1e-7  5e-8
    result_data = []
    max_diff = 20
    min_hgt_len = 500

    result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"
    # result_dir = "/mnt/d/HGT/seq_ana/homo_filter/"
    identified_hgt = "/mnt/d/HGT/seq_ana/identified_event.csv"
    database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"

    # meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    # design_file = "/mnt/d/HGT/time_lines/sample_design.tsv"
    # result_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    # identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    # database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna"
    # sra_sample_dict = read_meta()
    # sample_individual_dict, sample_time_point = read_design()

    # meta_data = "//mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/SRP366030.csv.txt"
    # design_file =  "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid//sample_design.tsv"
    # result_dir = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/hgt/result/"
    # identified_hgt = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/Hybrid/match/SRP366030.identified_event.csv"
    # database = "/mnt/delta_WS_1/wangshuai/02.HGT/detection/reference/UHGG_reference.formate.fna"


    tim = Match()
    tim.read_samples()
    print ("load is done.")
    sample_list = list(tim.cohort_data.keys())
    # sample_list = ["SRR18491248", "SRR18490939", "SRR18491317", "SRR18491328", "SRR18491254"]
    # sample_list = ["SRR18491277"] # , "SRR18490984"
    # print (sample_list)
    sample_index = 0
    for sample in sample_list:
        # if sample != "CCIS15704761ST-4-0":
        #     continue
        print ("processed %s sample."%(sample_index))
        tim.match_each_sample(sample)
        sample_index += 1
    # tim.match_each_sample('SRR18491328')

    df = pd.DataFrame(result_data, columns = ["sample", "receptor", "insert_locus", "donor", "delete_start", "delete_end", "reverse_flag"])
    df.to_csv(identified_hgt, sep=',')

