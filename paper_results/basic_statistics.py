
import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from ete3 import Tree

from mechanism_taxonomy import Taxonomy
from HGT_network import read_phenotype

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

        self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        self.pair_tag = "&".join(sorted([self.from_ref_genome, self.to_ref_genome]))
        self.abundance = None
        self.split_abundance = None

        self.read = None

class Basic_count():

    def __init__(self):
        self.cohort_data = {}
        
    def read_samples(self):
        cohort_count = {}
        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")

        # os.system(f"ls {tgs_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        # os.system(f"ls {wenkui_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            if len(my_bkps) > 0:
                self.cohort_data[sra_id] = my_bkps
                if sra_id not in phenotype_dict:
                    os.system("rm %s"%(acc_file))
                    continue
                cohort = phenotype_dict[sra_id][0]
                if cohort not in cohort_count:
                    cohort_count[cohort] = 0
                cohort_count[cohort] += 1
        print ("data loaded", cohort_count, sum(list(cohort_count.values())))
    
    def read_bkp(self, bkp_file):
        my_bkps = []
        f = open(bkp_file)
        all_rows = csv.reader(f)
        sample_dict = {}
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
                sample_dict[eb.hgt_tag] = 1
                my_bkps.append(eb)
        for eb in my_bkps:
            eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
        f.close()
        return my_bkps

    def count(self):
        # bkp number distribution across all the samples
        # count mean and median value of a sample
        data  = []
        bkp_count_list = []
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            bkp_count_list.append(len(sample_bkp_list))
            data.append([sample, len(sample_bkp_list)])
        sorted_bkp_count_list = sorted(bkp_count_list)
        print ("sample num is %s, mean bkp count is %s, median bkp count is %s, minimum bkp count is %s,\
         max bkp count is %s."%(len(self.cohort_data), np.mean(bkp_count_list), np.median(bkp_count_list),\
          sorted_bkp_count_list[0], sorted_bkp_count_list[-1]))

        df = pd.DataFrame(data, columns = ["Sample", "Bkp_count"])
        df.to_csv('/mnt/d/R_script_files/basic_stasitcis_count.csv', sep=',')

    def count_genome_pair(self):
        # genome_pair_dict = {}
        # for sample in self.cohort_data:
        #     sample_bkp_list = self.cohort_data[sample]
        #     sample_genome_pair = {}
        #     for bkp in sample_bkp_list:
        #         if bkp.pair_tag not in sample_genome_pair:
        #             sample_genome_pair[bkp.pair_tag] = 0
        #         sample_genome_pair[bkp.pair_tag] += 1
        #     for pair_tag in sample_genome_pair:
        #         if pair_tag not in genome_pair_dict:
        #             genome_pair_dict[pair_tag] = []
        #         genome_pair_dict[pair_tag].append(sample_genome_pair[pair_tag]/len(sample_bkp_list))

        genome_pair_dict = {}
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            sample_genome_pair = {}
            for bkp in sample_bkp_list:
                if bkp.pair_tag not in genome_pair_dict:
                    genome_pair_dict[bkp.pair_tag] = 0
                genome_pair_dict[bkp.pair_tag] += 1
        return genome_pair_dict
      
    def compare_intra_inter(self, level):
        # count the intra-taxa HGT frequency in each sample
        inter_freq = []
        intra_freq = []
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            inter_num, intra_num = 0, 0
            for bkp in sample_bkp_list:
                s1 = get_genome_taxa(bkp.from_ref_genome, level)
                s2 = get_genome_taxa(bkp.to_ref_genome, level)
                if s1[1:] == '__' or  s2[1:] == '__':
                    continue
                if s1 == s2 :
                    intra_num += 1
                else:
                    inter_num += 1
            if intra_num + inter_num == 0:
                print (sample.ID)
            inter_freq.append(inter_num/(intra_num + inter_num))
            intra_freq.append(intra_num/(intra_num + inter_num))
            intra_freq_data.append([intra_num/(intra_num + inter_num), level_list[level-1]])
        U1, p = mannwhitneyu(inter_freq, intra_freq)
        print (level, p, np.mean(intra_freq), np.median(intra_freq))

    def sort_taxa_by_freq(self, level):
        # sort taxa by its mean frequency in each sample
        total_freq_dict = {}
        for sample in self.cohort_data:
            sample_bkp_list = self.cohort_data[sample]
            sample_count = {}
            for bkp in sample_bkp_list:
                s1 = get_genome_taxa(bkp.from_ref_genome, level)
                s2 = get_genome_taxa(bkp.to_ref_genome, level)
                if s1[1:] == '__' or  s2[1:] == '__':
                    continue
                if s1 not in sample_count:
                    sample_count[s1] = 0
                if s2 not in sample_count:
                    sample_count[s2] = 0
                sample_count[s1] += 1
                sample_count[s2] += 1

            sample_freq = get_freq(sample_count)
            for taxa in sample_freq:
                if taxa not in total_freq_dict:
                    total_freq_dict[taxa] = []
                total_freq_dict[taxa].append(sample_freq[taxa])
        sample_num = len(self.cohort_data)
        sorted_dict = count_mean_freq(total_freq_dict, sample_num)
        return sorted_dict

    def count_inter_taxa_HGT(self):
        select_num = 10

        for level in range(1, 7):
            # if level >= 4:
            #     select_num = 20
            sorted_count = self.sort_taxa_by_freq(level)
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
                for sample in self.cohort_data:
                    for bkp in self.cohort_data[sample]:
                        from_tax = get_genome_taxa(bkp.from_ref_genome, level)
                        to_tax = get_genome_taxa(bkp.to_ref_genome, level)

                        if from_tax[1:] == '__' or  to_tax[1:] == '__':
                            continue

                        if from_tax not in index_dict:
                            index_dict[from_tax] = len(class_list)-1
                        if to_tax not in index_dict:
                            index_dict[to_tax] = len(class_list)-1

                        if index_dict[from_tax] == i:
                            j = index_dict[to_tax]
                            if j < i:
                                one_class[j] += 1
                        if index_dict[to_tax] == i:
                            j = index_dict[from_tax]
                            if j < i:
                                one_class[j] += 1
                data.append(one_class)
            for j in range(len(head_name)):
                head_name[j] = "_".join(head_name[j].split())
            df = pd.DataFrame(data, columns = head_name)
            df.index = head_name
            df.to_csv('/mnt/d/R_script_files/inter_taxa_files/inter_taxa_%s.csv'%(level_list[level-1]), sep='\t')
            print (level_list[level-1], "done")       
            
def get_genome_taxa(genome, level):
    if level < 7:
        # g1 = get_pure_genome(genome)
        taxa = taxonomy.taxonomy_dict[genome].split(";")[level]
        return taxa
    else:
        return genome

def get_freq(count_dict):
    total_count = sum(list(count_dict.values()))
    freq_dict = {}
    for taxa in count_dict:
        freq_dict[taxa] = count_dict[taxa]/total_count
    return freq_dict

def count_mean_freq(raw_dict, sample_num):
    mean_freq_dict = {}
    for taxa in raw_dict:
        if len(raw_dict[taxa])/sample_num < 0.1: # the taxa must exists in more than 10% of samples
            continue
        mean_freq_dict[taxa] = np.mean(raw_dict[taxa])
    sorted_dict = sorted(mean_freq_dict.items(), key=lambda item: item[1], reverse = True)
    print ("<<<<<<<<<<<<<<<<", len(sorted_dict))

    return sorted_dict

def prepare_tree_bk(): # rename genome to species name
    sorted_dict = ba.sort_taxa_by_freq(8) # 8 means genome level
    print ("total genome number is", len(sorted_dict))
    extracted_genome = []
    freq_dict = {}
    for i in range(300):
        extracted_genome.append(sorted_dict[i][0])
        freq_dict[sorted_dict[i][0]] = float(sorted_dict[i][1])
    # Load a tree structure from a newick file.
    f = open("/mnt/d/HGT/time_lines/distribution/hgt_tree.nwk", 'w')
    a = open("/mnt/d/HGT/time_lines/distribution/tree_annotation.txt", 'w')
    b = open("/mnt/d/HGT/time_lines/distribution/bar_annotation.txt", 'w')
    c = open("/mnt/d/HGT/time_lines/distribution/connection.txt", 'w')
    t = Tree("/mnt/d/HGT/time_lines/distribution/bac120_iqtree.nwk")
    t.prune(extracted_genome) # only keep the nodes with top HGT freq
    dark2=['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D']
    set1=['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']
    color_palette = dark2 + set1
    phylum_dict = {}
    index_dict = {'Firmicutes_A': 0, 'Bacteroidota': 1, 'Firmicutes': 2, 'Proteobacteria': 3, 'Actinobacteriota': 4, 'Firmicutes_C': 5, 'Verrucomicrobiota': 6}
    print ("TREE_COLORS\nSEPARATOR SPACE\nDATA\n", file = a)
    print ("DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,my_data\nCOLOR,#45818e\nDATA\n", file = b)
    node_name_dict = {}

    for node in t.traverse("postorder"):
        print (node.name)
        if node.name == '':
            continue
        genome_name = node.name
        species = get_genome_taxa(node.name, 6)
        phylum = get_genome_taxa(node.name, 1)[3:]
        if phylum not in phylum_dict:
            phylum_dict[phylum] = len(phylum_dict)
        if species[1:] == '__':
            species = node.name
        else:
            species = species[3:]
            species = "_".join(species.split())
        node.name = species
        node_name_dict[genome_name] = node.name
        print (f"{species} range {color_palette[index_dict[phylum]]} {phylum}", file = a)
        print (f"{node.name},{freq_dict[genome_name]},label1", file = b)
    # print (t)
    print (t.write(), file = f)
    f.close()
    a.close()
    b.close()
    c.close()
    print (phylum_dict)

def read_file_to_string(file_path):
    with open(file_path, 'r') as f:
        file_contents = f.read()
    return file_contents

def prepare_tree(): # just use genome name
    gradient_template = "/mnt/d/HGT/time_lines/distribution/dataset_gradient_template.txt"
    # Load a tree structure from a newick file.
    f = open("/mnt/d/HGT/time_lines/distribution/hgt_tree.nwk", 'w')
    a = open("/mnt/d/HGT/time_lines/distribution/tree_annotation.txt", 'w')
    b = open("/mnt/d/HGT/time_lines/distribution/bar_annotation.txt", 'w')
    c = open("/mnt/d/HGT/time_lines/distribution/connection.txt", 'w')
    g = open("/mnt/d/HGT/time_lines/distribution/gradient_annotation.txt", 'w')
    t = Tree("/mnt/d/HGT/time_lines/distribution/bac120_iqtree.nwk")


    sorted_dict = ba.sort_taxa_by_freq(8) # 8 means genome level
    print ("total genome number is", len(sorted_dict))
    extracted_genome = []
    extracted_genome_dict = {}
    genome_num = len(sorted_dict)
    freq_dict = {}
    for i in range(genome_num): # choose node number
        if len(t.search_nodes(name=sorted_dict[i][0])) == 0: # genome not in a tree
            continue
        extracted_genome.append(sorted_dict[i][0])
        extracted_genome_dict[sorted_dict[i][0]] = 1
        freq_dict[sorted_dict[i][0]] = float(sorted_dict[i][1])

    print ("total genome number in the tree is", len(extracted_genome))
    t.prune(extracted_genome) # only keep the nodes with top HGT freq
    dark2=['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D']
    set1=['#E41A1C', '#377EB8','#FFFF33', '#984EA3', '#FF7F00',  '#A65628', '#F781BF', '#999999', '#4DAF4A']
    color_palette = dark2 + set1
    index_dict = {'Firmicutes_A': 0, 'Bacteroidota': 1, 'Firmicutes': 2, 'Proteobacteria': 3, 'Actinobacteriota': 4, 'Firmicutes_C': 5, 'Verrucomicrobiota': 6}
    print ("TREE_COLORS\nSEPARATOR SPACE\nDATA\n", file = a)
    print ("DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,my_data\nCOLOR,#45818e\nDATA\n", file = b)
    print (read_file_to_string(gradient_template), file = g)
    node_name_dict = {}

    print (extracted_genome[:10])
    for node in t.traverse("postorder"):
        # print (node.name)
        if node.name == '':
            continue
        genome_name = node.name
        phylum = get_genome_taxa(node.name, 1)[3:]
        if phylum not in index_dict:
            index_dict[phylum] = len(index_dict)
        print (f"{genome_name} range {color_palette[index_dict[phylum]]} {phylum}", file = a)
        print (f"{node.name},{freq_dict[genome_name]},label1", file = b)
        print (f"{node.name} {np.log(freq_dict[genome_name])}", file = g)
    # print (t)
    print (t.write(), file = f)
    count_pair(ba.cohort_data, extracted_genome_dict, c)

    f.close()
    a.close()
    b.close()
    c.close()
    g.close()
    print (index_dict, len(index_dict))
    
def count_pair(cohort_data, extracted_genome_dict, connection_flag):
    print ("DATASET_CONNECTION\nSEPARATOR COMMA\nDATASET_LABEL,example connections\nCOLOR,#ff0ff0\nDRAW_ARROWS,0\nARROW_SIZE,20\nLOOP_SIZE,100\nMAXIMUM_LINE_WIDTH,10\nCURVE_ANGLE,0\nCENTER_CURVES,1\nALIGN_TO_LABELS,1\nDATA\n", file =connection_flag )
    # set1=['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']
    set1=["orange", "red", "blue", "green", "purple", "pink", "yellow"]
    pair_dict = {}
    print ("sample num:", len(cohort_data))
    for sample in cohort_data:
        for bkp in cohort_data[sample]:
            from_tax = get_genome_taxa(bkp.from_ref_genome, 8)
            to_tax = get_genome_taxa(bkp.to_ref_genome, 8)
            if from_tax not in extracted_genome_dict or to_tax not in extracted_genome_dict:
                continue
            pair_name = "&".join(sorted([from_tax, to_tax]))
            if pair_name not in pair_dict:
                pair_dict[pair_name] = 0
            pair_dict[pair_name] += 1
    for pair_name in pair_dict:
        array = pair_name.split("&")
        g1 = array[0]
        g2 = array[1]
        phylum1 = get_genome_taxa(g1, 1)[3:]
        phylum2 = get_genome_taxa(g2, 1)[3:]
        class1 = get_genome_taxa(g1, 2)[3:]
        class2 = get_genome_taxa(g2, 2)[3:]
        order1 = get_genome_taxa(g1, 3)[3:]
        order2 = get_genome_taxa(g2, 3)[3:]
        family1 = get_genome_taxa(g1, 4)[3:]
        family2 = get_genome_taxa(g2, 4)[3:]
        genus1 = get_genome_taxa(g1, 5)[3:]
        genus2 = get_genome_taxa(g2, 5)[3:]
        if genus1 != '' and genus1 == genus2:
            color = set1[0]
        elif family1 != '' and family1 == family2:
            color = set1[1]
        elif order1 != '' and order1 == order2:
            color = set1[2]
        elif class1 != '' and class1 == class2:
            color = set1[3]
        elif phylum1 != '' and phylum1 == phylum2:
            color = set1[4]
        else:
            color = set1[5]
        
        print (f"{g1},{g2},{pair_dict[pair_name]},{color},normal,", file = connection_flag)

def cal_corr(genome_pair_dict):
    t = Tree("/mnt/d/HGT/time_lines/distribution/bac120_iqtree.nwk")
    data = []
    dist_list, freq_list = [], []
    dist_bin_dict = {}

    total_num = sum(list(genome_pair_dict.values()))
    genome_pair_freq = {}
    for genome_pair in genome_pair_dict:
        freq = genome_pair_dict[genome_pair]/total_num  # the frequency of HGT between genome 1 and genome 2

        genome_pair_freq[genome_pair] = freq
        genome_arr = genome_pair.split("&")
        genome1 = genome_arr[0]
        genome2 = genome_arr[1]
        ### cal the phylogenetic distance between genome 1 and genome 2

        if len(t.search_nodes(name=genome1)) == 0 or len(t.search_nodes(name=genome2)) == 0: # genome not in a tree
            continue

        # Get the distance between two nodes
        node1 = t.get_leaves_by_name(genome1)[0]
        node2 = t.get_leaves_by_name(genome2)[0]
        distance = node1.get_distance(node2)

        # print (genome1, genome2, distance, freq)
        dist_list.append(distance)
        freq_list.append(freq)
        data.append([genome1, genome2, distance, freq])

        dist_bin = int(distance/0.5) 
        if dist_bin not in dist_bin_dict:
            dist_bin_dict[dist_bin] = 0
        dist_bin_dict[dist_bin] += genome_pair_dict[genome_pair]
        # break

    res = stats.spearmanr(dist_list, freq_list)
    correlation = res.correlation
    print (res, "spearmanr correlation", correlation)
    df = pd.DataFrame(data, columns = ["genome1", "genome2", "distance", "frequency"])
    df.to_csv('/mnt/d/R_script_files/distance_frequency.csv', sep=',')

    res = stats.spearmanr(list(dist_bin_dict.keys()), list(dist_bin_dict.values()))
    correlation = res.correlation
    print (res, "spearmanr correlation", correlation)

    res = stats.spearmanr((np.array(list(dist_bin_dict.keys())))*0.5, list(dist_bin_dict.values()))
    correlation = res.correlation
    print (res, "spearmanr correlation", correlation)

    data = []
    for dist_bin in  dist_bin_dict:
        data.append([dist_bin, 0.25+(dist_bin)*0.5, dist_bin_dict[dist_bin], dist_bin_dict[dist_bin]/total_num])
    df = pd.DataFrame(data, columns = ["index", "distance", "count", "frequency"])
    df.to_csv('/mnt/d/R_script_files/distance_frequency_bin.csv', sep=',')


if __name__ == "__main__":
    bin_size = 100
    abun_cutoff = 1e-7  #1e-7

    taxonomy = Taxonomy()
    phenotype_dict = read_phenotype()


    # hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    # tgs_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    # wenkui_dir = "/mnt/d/breakpoints/HGT/CRC/wenkui/"


    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"

    ba = Basic_count()
    ba.read_samples()
    # ba.count()  
    # ba.count_inter_taxa_HGT()

    # ################## cal the correlation between HGT frequency and phylogenetic distance
    genome_pair_dict = ba.count_genome_pair()
    cal_corr(genome_pair_dict)





    # ######## just sort taxa by HGT freq
    # taxa_sort_data = []
    # for level in range(1, 7):
    #     sorted_dict = ba.sort_taxa_by_freq(level)
    #     top_sum = 0
    #     for i in range(5):
    #         print (i, sorted_dict[i][0], sorted_dict[i][1])
    #         taxa_sort_data.append([sorted_dict[i][0], sorted_dict[i][1], level_list[level-1]])
    #         top_sum += sorted_dict[i][1]
    #     taxa_sort_data.append([level_list[level-1][0]+"__other", 1-top_sum, level_list[level-1]])

    # df = pd.DataFrame(taxa_sort_data, columns = ["Taxa", "Frequency", "Level"])
    # df.to_csv('/mnt/d/R_script_files/taxa_sort.csv', sep=',')

    #### prepare count tree
    # prepare_tree()

    # ######## get intra-taxa HGT freq
    # intra_freq_data = []
    # for level in range(1, 6):
    #     ba.compare_intra_inter(level)
    # #     ba.sort_inter_taxa_by_freq(level)

    # df = pd.DataFrame(intra_freq_data, columns = ["Frequency", "Level"])
    # df.to_csv('/mnt/d/R_script_files/intra_freq.csv', sep=',')


    ######## get inter-taxa HGT count
    # ba.count_inter_taxa_HGT()

    