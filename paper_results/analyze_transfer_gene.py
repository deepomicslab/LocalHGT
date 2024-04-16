#!/usr/bin/env python3

"""
Analyze the function of HGT-related genes and analyze the HGT transfer patterns
"""

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
from collections import Counter, defaultdict
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from Bio import SeqIO
from collections import Counter

COG_annotation = ["A: RNA processing and modification", "B: Chromatin structure and dynamics", "C: Energy production and conversion", "D: Cell cycle control, cell division, chromosome partitioning", "E: Amino acid transport and metabolism", "F: Nucleotide transport and metabolism", "G: Carbohydrate transport and metabolism", "H: Coenzyme transport and metabolism", "I: Lipid transport and metabolism", "J: Translation, ribosomal structure and biogenesis", "K: Transcription", "L: Replication, recombination and repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: Posttranslational modification, protein turnover, chaperones", "P: Inorganic ion transport and metabolism", "Q: Secondary metabolites biosynthesis, transport and catabolism", "R: General function prediction only", "S: Function unknown", "T: Signal transduction mechanisms", "U: Intracellular trafficking, secretion, and vesicular transport", "V: Defense mechanisms", "W: Extracellular structures", "Y: Nuclear structure", "Z: Cytoskeleton"]

def get_COG_dict():
    COG_dict = {}
    for anno in COG_annotation:
        arr = anno.split(":")
        name = arr[0]
        COG_dict[name] = anno  #arr[1].strip()

    COG_profile_dict = {}
    raw_dict = {"Metabolism":"QPIHFEGC", "Cellular processes and signaling":"XOUZNMTVDWY", "Information storage and Processing":"ABLKJ", "Function unknown":"RS"}
    for key in raw_dict:
        for i in range(len(raw_dict[key])):
            COG_profile_dict[raw_dict[key][i]] = key

    return COG_dict, COG_profile_dict

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
        self.abundance = None
        self.split_abundance = None

        self.read = None
        self.score = float(list[11])

        # self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        # self.pair_tag = "&".join(sorted([self.from_ref_genome, self.to_ref_genome]))

class Event(object):

    def __init__(self, array):
        self.sample = array[1]
        self.ins_genome = array[2]
        self.ins_genome_pure = "_".join(self.ins_genome.split("_")[:-1])
        self.ins_pos = int(array[3])
        self.del_genome = array[4]
        self.del_start = int(array[5])
        self.del_end = int(array[6])
        self.reverse_flag = array[7]
        self.IS_flag = None
        self.Transposon_flag = None

        bin_size = 100
        self.tag = self.ins_genome + "&" + str(round(self.ins_pos/bin_size)) + "&" + self.del_genome + "&" + \
            str(round(self.del_start/bin_size)) + "&" + str(round(self.ins_pos/bin_size))

    def check_IS(self, min_gene_frac, annotation):
        transfer_interval = [self.del_start, self.del_end]
        if self.del_genome in annotation.gene_annotation:
            self.IS_flag = True
            self.Transposon_flag = True
            # print ("Event", self.del_genome, transfer_interval)
            gene_intervals = annotation.gene_annotation[self.del_genome]["intervals"]
            for gene_interval in gene_intervals:
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], [transfer_interval], min_gene_frac)
                if not locate_transfer_flag:
                    continue
                gene_anno_dict = annotation.gene_annotation[self.del_genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "product" not in gene_anno_dict:
                    gene_anno_dict["product"] = "NA"
                if not re.search("IS[0-9]", gene_anno_dict["product"]):
                    self.IS_flag = False
                product_classification = annotation.classify_product(gene_anno_dict["product"])
                if product_classification != "Transposon":
                    self.Transposon_flag = False
                # print (gene_interval, gene_anno_dict["product"], self.IS_flag, self.Transposon_flag, sep = "\t")
        else:
            self.IS_flag = False
        # print ("Event", self.del_genome, transfer_interval, self.IS_flag)
        # print ("<<<<<<<<<<<<<<<<<<<<\n")

def get_gene_lengths(identified_hgt):
    HGT_event_list = []
    contig_lengths = []
    for line in open(identified_hgt):
        array = line.strip().split(",")
        if array[1] == "sample":
            continue
        gene_len = int(array[6]) - int(array[5])
        event = Event(array)
        HGT_event_list.append(event)
        # if gene_len == 1543502:
        #     print (array)
        contig_lengths.append(gene_len)

    print (np.mean(contig_lengths), np.median(contig_lengths), len(contig_lengths), min(contig_lengths), max(contig_lengths), sum(contig_lengths))
    print (sorted(contig_lengths)[-10:])

    data = []
    for a in contig_lengths:
        data.append([a, "mm"])
    df = pd.DataFrame(data, columns = ["Length", "mm"])
    df.to_csv("/mnt/d/R_script_files//gene_len.csv", sep=',')   

    print ("Finish reading events, num is", len(HGT_event_list))
    return HGT_event_list

class Annotation():

    def __init__(self, gff):
        self.gff = gff
        self.near = 100
        self.min_gene_frac = 0.5
        self.gene_annotation = {}
        self.gene_classification = {}
        self.pattern_dict = {}
        self.init_pattern_dict()

    def read_gff(self):

        f = open(self.gff)
        for line in f:
            array = line.split("\t")
            genome = array[0]
            g_type = array[2]
            detail = array[8].strip()
            start = int(array[3])
            end = int(array[4])
            if genome not in self.gene_annotation:
                self.gene_annotation[genome] = {}
                self.gene_annotation[genome]["intervals"] = []
            self.gene_annotation[genome]["intervals"].append([start, end])
            self.gene_annotation[genome][str(start)+ "_"+str(end)] = self.understand_gene(detail)
        f.close() 

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
            gene_ID, product = '', ''
            gene_anno_dict = self.gene_annotation[genome][str(inter[0]) + "_" + str(inter[1])]
            if "ID" in gene_anno_dict:
                gene_ID = gene_anno_dict["ID"]
            if "product" in gene_anno_dict:
                product = gene_anno_dict["product"]

            CDS_length = int(inter[1]) - int(inter[0])
            if inter[0] >= gene_interval[0] and inter[0] <= gene_interval[1] and (gene_interval[1] - inter[0])/CDS_length > self.min_gene_frac :
                genes_around.append(product)
            elif inter[0] <= gene_interval[0] and inter[1] >= gene_interval[0] and (inter[1] - gene_interval[0])/CDS_length > self.min_gene_frac:
                genes_around.append(product)

        if len(genes_around) > 0:
            # print (genome, gene_interval, gene_interval[1] - gene_interval[0], CDS_length, genes_around)
            return genes_around
        else:
            return []

    def understand_gene(self, detail):
        array = detail.split(";")
        anno_dict = {}
        for arr in array:
            category = arr.split("=")
            anno_dict[category[0]] = category[1]
        return anno_dict

    def init_pattern_dict(self):
        Transposon = "transpos*; insertion; resolv*; Tra[A-Z]; Tra[0-9]; IS[0-9]; conjugate transposon"
        Plasmid = "resolv*; relax*; conjug*; trb; mob*; plasmid; type IV; toxin; chromosome partitioning; chromosome segregation"
        Phage ="capsid; phage; tail; head; tape measure; antitermination; antiterminatio"
        Other_HGT_machinery=" integrase; excision*; exonuclease; recomb; toxin; CRISPR; restrict*; resolv*; topoisomerase; reverse transcrip"
        
        Carbohydrate_active =  "Genes present in the CAZY database; glycosyltransferase; glycoside hydrolase; xylan; monooxygenase; rhamnos*; cellulose; sialidase; *ose; acetylglucosaminidase; cellobiose; galact*; fructose; aldose; starch; mannose; mannan*; glucan; lyase; glycosyltransferase; glycosidase; pectin; SusD; SusC; fructokinase; galacto*; arabino*"
        antibiotic_resistance =  "Genes present in the ARDB; multidrug; azole resistance; antibiotic resistance; TetR; tetracycline resistance; VanZ; betalactam*; beta-lactam; antimicrob*; lantibio*"
        # hypothetical_protein = "hypothetical protein"
        
        gene_classification = {}
        self.pattern_dict["Transposon"] = get(Transposon)
        self.pattern_dict["Plasmid"] = get(Plasmid)
        self.pattern_dict["Phage"] = get(Phage)
        self.pattern_dict["Other_HGT_mechanisms"] = get(Other_HGT_machinery)
        self.pattern_dict["CAZYmes"] = get(Carbohydrate_active)
        self.pattern_dict["Antibiotic resistance"] = get(antibiotic_resistance)
        # self.pattern_dict["hypothetical protein"] = get(hypothetical_protein)
    
    def classify_product_bk(self, product):
        product_classification = 'unclassified'
        for cla in self.pattern_dict:
            # print (cla, pattern_dict[cla])
            pattern = re.compile(self.pattern_dict[cla]) #, re.I
            if pattern.search(product):
                product_classification = cla
        return product_classification

    def classify_product(self, product):
        # Transposon pattern
        transposon_pattern = re.compile(r"transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon")

        # Plasmid pattern
        plasmid_pattern = re.compile(r"relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation")

        # Phage pattern
        phage_pattern = re.compile(r"capsid|phage|tail|head|tape measure|antiterminatio")

        # Other HGT mechanisms pattern
        hgt_pattern = re.compile(r"integrase|excision\S*|exonuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip")

        # Carbohydrate active pattern
        carbohydrate_pattern = re.compile(r"glycosyltransferase|glycoside hydrolase|xylan|monooxygenase|rhamnos\S*|cellulose|sialidase|\S*ose($|\s|\-)|acetylglucosaminidase|cellobiose|galact\S*|fructose|aldose|starch|mannose|mannan\S*|glucan|lyase|glycosyltransferase|glycosidase|pectin|SusD|SusC|fructokinase|galacto\S*|arabino\S*")

        # Antibiotic resistance pattern
        antibiotic_pattern = re.compile(r"azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*")

        product_classification = 'unclassified'

        # Search for matches
        if plasmid_pattern.search(product):
            # print("Found a plasmid!")
            product_classification = "plasmid"
        if phage_pattern.search(product):
            # print("Found a phage!")
            product_classification = "phage"
        if transposon_pattern.search(product):
            # print("Found a transposon!")
            product_classification = "transposon"
        if hgt_pattern.search(product):
            # print("Found an HGT mechanism!")
            product_classification = "Other_HGT_mechanisms"
        if carbohydrate_pattern.search(product):
            # print("Found a carbohydrate active enzyme!")
            product_classification = "CAZYmes"
        if antibiotic_pattern.search(product):
            # print("Found an antibiotic resistance gene!")
            product_classification = "ARG"

        return product_classification
    
    def if_IS(self, product):
        pattern = re.compile("IS[0-9]")
        if pattern.search(product):
            return True
        else:
            return False

class Transfer_times():

    def __init__(self):
        self.HGT_event_dict = {}
        self.max_diff = 50

    def read_events(self, identified_hgt):
        i = 0
        for line in open(identified_hgt):
            array = line.strip().split(",")
            if array[1] == "sample":
                continue
            event = Event(array)
            sample = array[1]
            # if phenotype_dict[sample][0] == "Time-series":
            #     continue
            # print (phenotype_dict[sample][0])
            i += 1

            if sample not in self.HGT_event_dict:
                self.HGT_event_dict[sample] = []
            self.HGT_event_dict[sample].append(event) 
        print ("sample num ", len(self.HGT_event_dict), "kept HGTs num is", i)

    def main_count_times(self):
        final_data = []
        final_data = self.count_times(final_data, "different")
        final_data = self.count_times(final_data, "same")
        df = pd.DataFrame(final_data, columns = ["times", "count", "frequency", "group"])
        df.to_csv("/mnt/d/R_script_files/transfer_times.csv", sep=',')  

    def all_count_times(self): 
        ## count how many times the sequence is transferred in all samples
        
        trans_times_dict_sample_count = defaultdict(set)
        record = defaultdict(int)
        trans_times_dict = defaultdict(set)
        multiple_sites_sequence = {} # record the transferred sequences that have multiple sites

        data = []
        window = 100
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:
                segment_tag = "&".join([event.del_genome, str(round(event.del_start/window)), str(round(event.del_end/window))])  #, event.ins_genome_pure
                insertion_site = event.ins_genome + "&" + str(round(event.ins_pos/window))
                trans_times_dict[segment_tag].add(insertion_site)

        site_number = []
        for segment_tag in trans_times_dict:
            # print (len(trans_times_dict[segment_tag]))
            site_number.append(len(trans_times_dict[segment_tag]))
            if len(trans_times_dict[segment_tag]) > 1:
                multiple_sites_sequence[segment_tag] = 1
        site_number_count = Counter(site_number)
        print (len(trans_times_dict), sum(site_number), site_number_count)

        data = []
        over_ten = 0
        for count in range(1, max(list(site_number_count.keys()))+1):
            if count < 11:
                data.append([count, site_number_count[count]])
            else:
                over_ten += site_number_count[count]
        print ("multiple rario", (len(trans_times_dict)-site_number_count[1])/len(trans_times_dict), (len(trans_times_dict)-site_number_count[1]))
        data.append([">10", over_ten])
        df = pd.DataFrame(data, columns = ["count", "number"])
        df.to_csv("/mnt/d/R_script_files/transfer_times_seq.csv", sep=',')  


        individual_multiple_trans = {}
        individual_multiple_trans_recipient = {}

        individual_multiple_segment_dict = {} # record the transferred sequence that have multiple target sites in an individual
        for sample in self.HGT_event_dict:
            sample_set = defaultdict(set)
            for event in self.HGT_event_dict[sample]:
                segment_tag = "&".join([event.del_genome, str(round(event.del_start/window)), str(round(event.del_end/window))])  #, event.ins_genome_pure
                insertion_site = event.ins_genome + "&" + str(round(event.ins_pos/window))
                sample_set[segment_tag].add(insertion_site)

            for segment in sample_set:
                if len(sample_set[segment]) > 1:
                    individual_multiple_segment_dict[segment] = 1
        print (len(individual_multiple_segment_dict), len(individual_multiple_segment_dict)/len(trans_times_dict))

        individual_multiple_segment_dict = {} # record the transferred sequence that have multiple target sites of a same recipient in an individual
        for sample in self.HGT_event_dict:
            sample_set = {}
            for event in self.HGT_event_dict[sample]:
                segment_tag = "&".join([event.del_genome, str(round(event.del_start/window)), str(round(event.del_end/window))])  #, event.ins_genome_pure
                insertion_site = event.ins_genome + "&" + str(round(event.ins_pos/window))
                if segment_tag not in sample_set:
                    sample_set[segment_tag] = {}
                if event.ins_genome_pure not in sample_set[segment_tag]:
                    sample_set[segment_tag][event.ins_genome_pure] = set()
                sample_set[segment_tag][event.ins_genome_pure].add(insertion_site)

            for segment in sample_set:
                for recipient_genome in sample_set[segment]:
                    if len(sample_set[segment][recipient_genome]) > 1:
                        individual_multiple_segment_dict[segment] = 1
        print (len(individual_multiple_segment_dict), len(individual_multiple_segment_dict)/len(trans_times_dict))


        HGTs_multiple_samples = defaultdict(set) # record how many individuals contain the HGTs related to the multiple-sites transferred sequences
        print ("No. of the transferred sequences that have multiple sites", len(multiple_sites_sequence))
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:
                segment_tag = "&".join([event.del_genome, str(round(event.del_start/window)), str(round(event.del_end/window))])  #, event.ins_genome_pure
                if segment_tag not in multiple_sites_sequence:
                    continue
                insertion_site = event.ins_genome + "&" + str(round(event.ins_pos/window))
                # trans_times_dict[segment_tag].add(insertion_site)
                HGT_event_tag = segment_tag + "&" + insertion_site
                if sample in sra2individual_ID:  ## avoid the impact of time-series data
                    print (sample)
                    sample = sra2individual_ID[sample]
                    print (sample)
                HGTs_multiple_samples[HGT_event_tag].add(sample)
        single = 0
        for HGT_event_tag in HGTs_multiple_samples:
            if len(HGTs_multiple_samples[HGT_event_tag]) == 1:
                single += 1
        print (len(HGTs_multiple_samples), single, single/len(HGTs_multiple_samples))
        
    def count_times(self, final_data, trans_type):
        ## count how many times the sequence is transferred in each sample
        
        trans_times_dict_sample_count = defaultdict(set)
        record = defaultdict(int)

        data = []
        for sample in self.HGT_event_dict:
            trans_times_dict = {}
            for event in self.HGT_event_dict[sample]:
                recorded_already = False
                for segment_tag in trans_times_dict:
                    array = segment_tag.split("&")
                    if trans_type == "different":
                        if array[0] == event.del_genome and abs(int(event.del_start) - int(array[1])) < self.max_diff and \
                            abs(int(event.del_end) - int(array[2])) < self.max_diff :#and array[3] == event.ins_genome_pure:
                            recorded_already = True
                            trans_times_dict[segment_tag] += 1
                            break

                    elif trans_type == "same":
                        if array[0] == event.del_genome and abs(int(event.del_start) - int(array[1])) < self.max_diff and \
                            abs(int(event.del_end) - int(array[2])) < self.max_diff and array[3] == event.ins_genome_pure:
                            recorded_already = True
                            trans_times_dict[segment_tag] += 1
                            break
                    else:
                        print ("wrong trans_type", trans_type)

                if not recorded_already:
                    segment_tag = "&".join([event.del_genome, str(event.del_start), str(event.del_end), event.ins_genome_pure])
                    trans_times_dict[segment_tag] = 1
            #print (sample, list(trans_times_dict.values()))
            for segment_tag in trans_times_dict:
                if trans_times_dict[segment_tag] > 1:
                    trans_times_dict_sample_count[segment_tag].add(sample)
                    record[segment_tag] += trans_times_dict[segment_tag]
            data += list(trans_times_dict.values())
        num = len(data)
        multiple_trans_num = 0
        counts = {}
        for item in data:
            if item in counts:
                counts[item] += 1
            else:
                counts[item] = 1
            if item > 1 :
                multiple_trans_num += 1
        
        for item in counts:
            final_data.append([item, counts[item], counts[item]/num, trans_type])
        print ("%s of segments are transferred multiple times."%(multiple_trans_num/num), multiple_trans_num, num-multiple_trans_num, num, sum(data)) 
        single_trans_num = num - multiple_trans_num
        print (sum(data)-single_trans_num, single_trans_num, (sum(data)-single_trans_num)/sum(data) )

        multiple_persons_HGT_num = 0
        all_HGT_num = 0
        for segment_tag in trans_times_dict_sample_count:
            HGT_num = record[segment_tag] 
            if len(trans_times_dict_sample_count[segment_tag]) > 1:
                multiple_persons_HGT_num += HGT_num
            all_HGT_num += HGT_num
        print ("uniq samples", multiple_persons_HGT_num, all_HGT_num, all_HGT_num - multiple_persons_HGT_num, 1-multiple_persons_HGT_num/all_HGT_num)
        return final_data

    def check_gene_type(self, multi_segment_tag):
        annotation = Annotation(gff)
        annotation.read_gff()

        contig_lengths = []
        IS_num = 0
        for segment_tag in multi_segment_tag:
            array = segment_tag.split("&")
            print (array)
            IS_flag = False
            product_list = annotation.given_seg(array[0], [int(array[1]), int(array[2])])
            contig_lengths.append(int(array[2])-int(array[1]))
            for product in product_list:
                if annotation.if_IS(product):
                    IS_flag = True
            if IS_flag:
                IS_num += 1
        print (IS_num, len(multi_segment_tag), IS_num/len(multi_segment_tag))
        print (np.mean(contig_lengths), np.median(contig_lengths), len(contig_lengths), min(contig_lengths), max(contig_lengths), sum(contig_lengths))
        print (sorted(contig_lengths)[-10:])

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

def count_product():
    HGT_event_list = get_gene_lengths(identified_hgt) ## get the gene length distribution of transferred seqs
    annotation = Annotation(gff)
    annotation.read_gff()

    all_products = []
    for event in HGT_event_list:
        all_products += annotation.given_seg(event.del_genome, [event.del_start, event.del_end])
        # break
    
    prod_type_dict = {}
    for product in all_products:
        prod_type = annotation.classify_product(product)
        # print (product, prod_type)
        if prod_type not in prod_type_dict:
            prod_type_dict[prod_type] = 0
        prod_type_dict[prod_type] += 1
    print (prod_type_dict)

def count_transfer_times():
    # get_gene_lengths(identified_hgt)  #  get transfer seq length distribution
    trans = Transfer_times()  
    trans.read_events(identified_hgt)
    # trans.main_count_times()
    trans.all_count_times()

def merge_intervals(intervals):
    # Sort the intervals by their start time
    intervals.sort(key=lambda x: x[0])
    
    # Initialize an empty list to store the merged intervals
    merged = []
    
    # Iterate through all the intervals
    for interval in intervals:
        # If the merged list is empty or the current interval doesn't overlap with the last merged interval
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            # If the current interval overlaps with the last merged interval, merge them
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    
    # Return the merged intervals
    return merged

class Extract_KO():

    def __init__(self, HGT_event_dict):
        self.HGT_event_dict = HGT_event_dict
        self.min_gene_frac = min_gene_frac
        self.transfer_regions = {}
        self.insert_sites = {}
        # self.no_transfer_regions = {}
        # self.no_ins_regions= {}
        self.transfer_kos = []
        self.no_transfer_kos = []
        self.insert_kos = []
        self.no_insert_kos = []
        self.near = 0
    
    def classify_regions(self):
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:

                event.check_IS(self.min_gene_frac, annotation)
                if remove_transposon_flag:
                    if event.Transposon_flag:
                        continue

                if event.del_genome not in self.transfer_regions:
                    self.transfer_regions[event.del_genome] = []
                if event.ins_genome not in self.insert_sites:
                    self.insert_sites[event.ins_genome] = []
                self.transfer_regions[event.del_genome].append([event.del_start, event.del_end])
                self.insert_sites[event.ins_genome].append(event.ins_pos)
        self.merge_transfer_region()
    
    def merge_transfer_region(self):
        for genome in self.transfer_regions:
            self.transfer_regions[genome] = merge_intervals(self.transfer_regions[genome])

    def classify_kos(self):  # gene in transfer region, or not
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")
                #### check if the gene locates in transfer region
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                if locate_transfer_flag:
                    self.transfer_kos += KEGG_list
                else:
                    self.no_transfer_kos += KEGG_list
                # print (KEGG_list, locate_transfer_flag)
            # break
        print (len(self.transfer_kos), len(self.no_transfer_kos))
        print_data(self.transfer_kos, transfer_ko_file)
        print_data(self.no_transfer_kos, no_transfer_ko_file)

    def collect_all_ko(self):
        all_ko = []
        for genome in annotation.gene_annotation:
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")
                all_ko += KEGG_list
        print_data(all_ko, all_kos_file)

    def classify_kos_insert(self):  # gene in insert site, or not
        for genome in self.insert_sites:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(self.insert_sites[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0]- self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    self.insert_kos += KEGG_list
                else:
                    self.no_insert_kos += KEGG_list
                # print (KEGG_list, locate_insert_flag)
            # break
        print (len(self.insert_kos), len(self.no_insert_kos))
        print_data(self.insert_kos, insert_ko_file)
        print_data(self.no_insert_kos, no_insert_ko_file)

    def classify_bkp_kos(self):
        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        bkp_dict = {}

        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]

            if only_healthy: # only check healthy persons
                pheno = phenotype_dict[sra_id]
                disease = pheno[1]  
                full_disease = pheno[2].split(";")

                if disease != "control" or full_disease[0] != "healthy":
                    continue

            my_bkps = self.read_bkp(acc_file)
            for bkp in my_bkps:
                if bkp.from_ref not in bkp_dict:
                    bkp_dict[bkp.from_ref] = []
                bkp_dict[bkp.from_ref].append(bkp.from_bkp)
                if bkp.to_ref not in bkp_dict:
                    bkp_dict[bkp.to_ref] = []
                bkp_dict[bkp.to_ref].append(bkp.to_bkp)
        
        bkp_ko = []
        no_bkp_ko = []
        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                if re.search("IS[0-9]", gene_anno_dict["product"]):  # IS element
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    bkp_ko += KEGG_list
                else:
                    no_bkp_ko += KEGG_list
        print_data(bkp_ko, bkp_ko_file)
        print_data(no_bkp_ko, no_bkp_ko_file)
    
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

class Extract_COG(Extract_KO):

    def __init__(self, HGT_event_dict):
        Extract_KO.__init__(self, HGT_event_dict)
        self.transfer_cog = []
        self.no_transfer_cog = []
        self.insert_cog = []
        self.no_insert_cog = []
        self.bkp_cog = []
        self.no_bkp_cog = []
        self.data = []

    def classify_cog(self):  # gene in transfer region, or not
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                # print (gene_anno_dict)
                if "COG" not in gene_anno_dict:
                    continue
                cog = gene_anno_dict["COG"]
                #### check if the gene locates in transfer region
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                if locate_transfer_flag:
                    for i in range(len(cog)):
                        self.transfer_cog += [cog[i]]
                else:
                    for i in range(len(cog)):
                        self.no_transfer_cog += [cog[i]]
                # print (KEGG_list, locate_transfer_flag)
            # break
        print (len(self.transfer_cog), len(self.no_transfer_cog))
        self.data = enrichment_analysis(self.transfer_cog, self.no_transfer_cog, "transfer", self.data)

    def classify_cog_insert(self):  # gene in insert site, or not
        for genome in self.insert_sites:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(self.insert_sites[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "COG" not in gene_anno_dict:
                    continue
                cog = gene_anno_dict["COG"]

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    for i in range(len(cog)):
                        self.insert_cog += [cog[i]]
                else:
                    for i in range(len(cog)):
                        self.no_insert_cog += [cog[i]]

        print (len(self.insert_cog), len(self.no_insert_cog))
        self.data = enrichment_analysis(self.insert_cog, self.no_insert_cog, "insert", self.data)

    def classify_bkp_cog(self):
        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        bkp_dict = {}

        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            for bkp in my_bkps:
                if bkp.from_ref not in bkp_dict:
                    bkp_dict[bkp.from_ref] = []
                bkp_dict[bkp.from_ref].append(bkp.from_bkp)
                if bkp.to_ref not in bkp_dict:
                    bkp_dict[bkp.to_ref] = []
                bkp_dict[bkp.to_ref].append(bkp.to_bkp)
        
        bkp_ko = []
        no_bkp_ko = []
        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]

                # if re.search("IS[0-9]", gene_anno_dict["product"]):  # IS element
                #     continue
                if remove_transposon_flag:
                    if annotation.classify_product(gene_anno_dict["product"]) == "Transposon":
                        continue

                if "COG" not in gene_anno_dict:
                    continue
                cog = gene_anno_dict["COG"]

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    for i in range(len(cog)):
                        self.bkp_cog += [cog[i]]
                else:
                    for i in range(len(cog)):
                        self.no_bkp_cog += [cog[i]]
        print (len(self.bkp_cog), len(self.no_bkp_cog))
        self.data = enrichment_analysis(self.bkp_cog, self.no_bkp_cog, "BKP", self.data)

    def main(self):
        self.classify_regions()
        self.classify_bkp_cog()
        self.classify_cog()
        self.classify_cog_insert()

        df = pd.DataFrame(self.data, columns = ["category", "category_detail", "p_value", "fold", "gene_num", "profile", "locus_type"])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected
        df.to_csv(cog_enrich, sep=',')
        print ("enriched COG num", len(self.data))

class Extract_product(Extract_KO):

    def __init__(self, HGT_event_dict):
        Extract_KO.__init__(self, HGT_event_dict)
        self.transfer_product = []
        self.no_transfer_product = []
        self.insert_product = []
        self.no_insert_product = []
        self.bkp_product = []
        self.no_bkp_product = []
        self.data = []

    def classify_product(self):  # gene in transfer region, or not
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                # print (gene_anno_dict)
                prod_category = annotation.classify_product(gene_anno_dict["product"])
                if prod_category == "unclassified":
                    continue

                #### check if the gene locates in transfer region
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                if locate_transfer_flag:
                    self.transfer_product += [prod_category]
                else:
                    self.no_transfer_product += [prod_category]
                # print (KEGG_list, locate_transfer_flag)
            # break
        print (len(self.transfer_product), len(self.no_transfer_product))
        self.data = enrichment_analysis_product(self.transfer_product, self.no_transfer_product, "transfer", self.data)

    def classify_product_insert(self):  # gene in insert site, or not
        for genome in self.insert_sites:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(self.insert_sites[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                prod_category = annotation.classify_product(gene_anno_dict["product"])

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    self.insert_product += [prod_category]
                else:
                    self.no_insert_product += [prod_category]

        print (len(self.insert_product), len(self.no_insert_product))
        self.data = enrichment_analysis_product(self.insert_product, self.no_insert_product, "insert", self.data)

    def classify_bkp_product(self):
        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        bkp_dict = {}

        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            for bkp in my_bkps:
                if bkp.from_ref not in bkp_dict:
                    bkp_dict[bkp.from_ref] = []
                bkp_dict[bkp.from_ref].append(bkp.from_bkp)
                if bkp.to_ref not in bkp_dict:
                    bkp_dict[bkp.to_ref] = []
                bkp_dict[bkp.to_ref].append(bkp.to_bkp)
        
        bkp_ko = []
        no_bkp_ko = []
        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                prod_category = annotation.classify_product(gene_anno_dict["product"])
                # if re.search("IS[0-9]", gene_anno_dict["product"]):  # IS element
                #     continue
                if remove_transposon_flag:
                    if annotation.classify_product(gene_anno_dict["product"]) == "Transposon":
                        continue

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    self.bkp_product += [prod_category]
                else:
                    self.no_bkp_product += [prod_category]
        print (len(self.bkp_product), len(self.no_bkp_product))
        self.data = enrichment_analysis_product(self.bkp_product, self.no_bkp_product, "bkp", self.data)

    def main(self):
        self.classify_regions()
        self.classify_bkp_product()
        self.classify_product()
        self.classify_product_insert()

        df = pd.DataFrame(self.data, columns = ["category", "p_value", "fold", "gene_num", "locus_type"])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected
        df.to_csv(product_enrich, sep=',')
        # print ("enriched product num", len(self.data))

class Extract_seq(Extract_KO):

    def __init__(self, HGT_event_dict):
        Extract_KO.__init__(self, HGT_event_dict)
        self.ref_fasta = Fasta(database)
        # self.records_dict = SeqIO.index(database, "fasta")

    def transfer_cds(self):  # gene in transfer region, or not
        self.classify_regions()
        f = open(transfer_cds_fasta, 'w')
        h = open(no_transfer_cds_fasta, "w")

        cds_index = 0
        conserve_cds_index = 0
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                transfer_seq = self.ref_fasta[genome][gene_interval[0]:gene_interval[1]].seq
                # transfer_seq = self.records_dict[genome].seq[gene_interval[0]-1:gene_interval[1]]
                if genome + "_" + str(gene_interval[0]) + "_" + str(gene_interval[1]) == "GUT_GENOME000269_2_220262_220672":
                    print ("Bug here", transfer_seq)
                if locate_transfer_flag:
                    cds_index += 1
                    cds_ID = genome + "_" + str(gene_interval[0]) + "_" + str(gene_interval[1]) + "_" + str(cds_index) 
                    print (">%s"%(cds_ID), file = f)
                    print (transfer_seq, file = f)
                else:
                    conserve_cds_index += 1
                    cds_ID = genome + "_" + str(gene_interval[0]) + "_" + str(gene_interval[1]) + "_" + str(conserve_cds_index) 
                    print (">%s"%(cds_ID), file = h)
                    print (transfer_seq, file = h)                    
        f.close()
        h.close()

def get_contig_lengths(fasta_file):   
    # Initialize an empty dictionary to store the contig lengths
    contig_lengths = {}
    
    # Parse the fasta file and loop over the records
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Store the contig ID and length in the dictionary
        contig_lengths[record.id] = len(record.seq)
    
    return contig_lengths

def read_blastn(blastn_file, fasta_file):
    contig_lengths = get_contig_lengths(fasta_file)
    # print (contig_lengths)
    map_cds = set()
    # Open the blastn results file
    with open(blastn_file, "r") as f:
        # Loop over each line in the file
        for line in f:
            # Split the line into fields
            fields = line.strip().split("\t")
            # Calculate the alignment length and proportion
            align_len = int(fields[3])
            query_name = fields[0]
            query_len = contig_lengths[query_name]
            align_prop = align_len / query_len
            # Check if the alignment proportion is larger than 0.5
            if align_prop > 0.5:
                # Extract the query name
                query_name = fields[0]
                # Print the query name
                # print(query_name)
                map_cds.add(query_name)
    print (len(map_cds), len(map_cds)/len(contig_lengths))
    return len(map_cds), len(contig_lengths)

def blast_main(fasta_file, db):
    command = f"blastn -db {db} -query {fasta_file} -outfmt 6 -out {fasta_file}.out -num_threads 8"
    # print (command)
    os.system(command)
    mapped_num, all_num = read_blastn(f"{fasta_file}.out", fasta_file)
    return mapped_num, all_num

def enrichment_analysis(my_list, background_list, locus_type, data):
    my_dict = Counter(my_list)
    background_dict = Counter(background_list)
    # data = []
    # for category in my_dict:
    for category in COG_dict:
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
        print (category, p_value, oddsratio, a, b, c, d) 
        data.append([category, COG_dict[category], p_value, oddsratio, a, COG_profile_dict[category], locus_type])
    return data

def enrichment_analysis_product(my_list, background_list, locus_type, data):
    my_dict = Counter(my_list)
    background_dict = Counter(background_list)
    # data = []
    for prod_category in background_dict:
        if prod_category in my_dict:
            a = my_dict[prod_category]
        else:
            a = 0
        b = len(my_list) - a
        if prod_category in background_dict:
            c = background_dict[prod_category]
        else:
            c = 0
        d = len(background_list) - c

        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
        print (locus_type, prod_category, p_value, oddsratio, a, b, c, d) 
        data.append([prod_category, p_value, oddsratio, a, locus_type])
    return data

def check_overlap(A, intervals, min_gene_frac):
    # Calculate the length of interval A
    A_length = A[1] - A[0]
    
    # Iterate through all the intervals in the list
    for interval in intervals:
        # Calculate the overlap between interval A and the current interval
        overlap = min(A[1], interval[1]) - max(A[0], interval[0])
        
        # If there is an overlap and the length of the overlap is longer than half of A, return True
        if overlap > 0 and overlap > A_length * min_gene_frac:
            return True
    
    # If there is no overlap that meets the criteria, return False
    return False

def print_data(my_set, file):
    f = open(file, 'w')
    for element in my_set:
        if element[:3] == "ko:":
            element = element[3:]
        print(element, end = "\n", file = f)
    f.close()   

class Classify(): # check the composition of transferred sequences

    def __init__(self, HGT_event_dict):
        self.HGT_event_dict = HGT_event_dict
        self.min_gene_frac = min_gene_frac

    def main(self):
        window = 100
        record_segment_tag = {}
        total_num, IS_num, trans_num = 0, 0, 0
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:
                segment_tag = "&".join([event.del_genome, str(round(event.del_start/window)), str(round(event.del_end/window))])  #, event.ins_genome_pure
                if segment_tag in record_segment_tag:
                    continue
                event.check_IS(self.min_gene_frac, annotation)
                if event.IS_flag:
                    IS_num += 1
                if event.Transposon_flag:
                    trans_num += 1
                total_num += 1
                record_segment_tag[segment_tag] = 1

            # break
        print (IS_num, total_num, IS_num/total_num, trans_num, trans_num/total_num, IS_num/trans_num)

def read_meta(meta_data):
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

def get_individual_ID(): ## for the hybrid dataset
    sra2individual_ID = {}
    sra_sample_dict = read_meta(meta_data)
    sample_individual_dict, sample_time_point = read_design()
    for sra in sra_sample_dict:
        sample_id = sra_sample_dict[sra]
        individual_id = sample_individual_dict[sample_id]
        sra2individual_ID[sra] = individual_id
    return sra2individual_ID

def count_uniq_event_ratio(HGT_event_dict): # cal the ratio of uniq HGT events which occur only in one sample
    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    sra_sample_dict = read_meta(meta_data)

    uniq_event_dict = defaultdict(set)
    for sample in HGT_event_dict:
        for event in HGT_event_dict[sample]:
            if sample in sra_sample_dict:
                sample = sra_sample_dict[sample]
            uniq_event_dict[event.tag].add(sample)
    all_event_num = len(uniq_event_dict)
    uniq_event_num = 0
    for event_tag in uniq_event_dict:
        if len(uniq_event_dict[event_tag]) == 1:
            uniq_event_num += 1
    print (uniq_event_num, all_event_num, uniq_event_num/all_event_num)

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
    abun_cutoff = 1e-7 
    min_gene_frac = 0.5

    # identified_hgt = "/mnt/d/HGT/seq_ana/bk/identified_event.csv"
    # identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    identified_hgt = "/mnt/d/HGT/seq_ana/identified_event.csv" # HGT event results detected by LocalHGT, can be downloaded from https://doi.org/10.5281/zenodo.10906354
    gff = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.gff" # gene annotation for UHGG, can be get from  https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/
    database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna" # UHGG v1 reference
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/" # HGT breakpoint results detected by LocalHGT, can be downloaded from https://doi.org/10.5281/zenodo.10906354
    all_kos_file = "/mnt/d/HGT/seq_ana/all_kos.txt" # to save all KOs
    transfer_cds_fasta = "/mnt/d/HGT/seq_ana/transfer_cds.fasta" # to save the CDS in transferred sequence
    no_transfer_cds_fasta = "/mnt/d/HGT/seq_ana/no_transfer_cds.fasta" # to save the CDS in non-transferred sequence

    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"  # meta data of the time-lines cohort
    design_file = "/mnt/d/HGT/time_lines/sample_design.tsv" # record the time point and corresponding individual of the time-lines cohort
    sra2individual_ID = get_individual_ID()


    transfer_ko_file = "/mnt/d/HGT/seq_ana/transfer_ko.txt" # to save the gene KO in transferred sequence
    no_transfer_ko_file = "/mnt/d/HGT/seq_ana/no_transfer_ko.txt" # to save the gene KO in non-transferred sequence

    insert_ko_file = "/mnt/d/HGT/seq_ana/insert_ko.txt" # to save the gene KO in insertion breakpoints of HGT events
    no_insert_ko_file = "/mnt/d/HGT/seq_ana/no_insert_ko.txt"

    bkp_ko_file = "/mnt/d/HGT/seq_ana/bkp_ko.txt" # to save the gene KO in all HGT breakpoints
    no_bkp_ko_file = "/mnt/d/HGT/seq_ana/no_bkp_ko.txt"

    cog_enrich = "/mnt/d/R_script_files/cog_enrich.csv"  # to save the enriched COG categories
    product_enrich = "/mnt/d/R_script_files/product_enrich.csv" # to save the enriched gene product categories
    # cog_insert = "/mnt/d/R_script_files/cog_insert.csv"

    COG_dict, COG_profile_dict = get_COG_dict()
    phenotype_dict = read_phenotype()

    remove_transposon_flag = False
    only_healthy = False  # only check healthy persons

    # annotation = Annotation(gff) # load gene annotation
    # annotation.read_gff()

    # count_transfer_times() # analyze gene transfer pattern
    get_gene_lengths(identified_hgt)  # get transferred sequence length distribution 

    # trans = Transfer_times()
    # trans.read_events(identified_hgt)

    # count_uniq_event_ratio(trans.HGT_event_dict)

    # extract_seq = Extract_seq(trans.HGT_event_dict)
    # extract_seq.transfer_cds()

    # extract = Extract_KO(trans.HGT_event_dict)
    # extract.classify_bkp_kos()
    # # extract.collect_all_ko()
    # extract.classify_regions()
    # extract.classify_kos()
    # extract.classify_kos_insert()

    # data = []
    # cog = Extract_COG(trans.HGT_event_dict)
    # cog.main()

    # data = []
    # prod = Extract_product(trans.HGT_event_dict)
    # prod.main()

    # classify = Classify(trans.HGT_event_dict) # count insertion sequence element and transposon frequency in the transferred sequences
    # classify.main()


    ###############to analyze if phage and plasmid is enriched in transferred sequences, not used anymore
    # phage_db = "/mnt/d/HGT/seq_ana/BlastDB/allprophage_DB"
    # plasmid_db = "/mnt/d/HGT/seq_ana/database/plsdb.fna"  
    # phage_trans, all_trans = blast_main(transfer_cds_fasta, plasmid_db)
    # phage_no, all_no = blast_main(no_transfer_cds_fasta, plasmid_db) 
    # a = phage_trans
    # b = all_trans - a
    # c = phage_no
    # d = all_no -c 
    # oddsratio, p_value = fisher_exact([[a, b], [c, d]])
    # print (oddsratio, p_value)







        
