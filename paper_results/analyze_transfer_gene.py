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

    # data = []
    # for a in contig_lengths:
    #     data.append([a, "mm"])
    # df = pd.DataFrame(data, columns = ["Length", "mm"])
    # df.to_csv("/mnt/d/R_script_files//gene_len.csv", sep=',')   

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
        Plasmid = "resolv*; relax*; conjug*; trb; mob*; plasmid; “type IV”; toxin; “chromosome partitioning”; “chromosome segregation”"
        Phage ="capsid; phage; tail; head; “tape measure”; antitermination"
        Other_HGT_machinery=" integrase; excision*; exonuclease; recomb; toxin; CRISPR; restrict*; resolv*; topoisomerase; reverse transcrip"
        Carbohydrate_active =  "Genes present in the CAZY database; glycosyltransferase; “glycoside hydrolase; xylan; monooxygenase; rhamnos*; cellulose; sialidase; *ose; acetylglucosaminidase; cellobiose; galact*; fructose; aldose; starch; mannose; mannan*; glucan; lyase; glycosyltransferase; glycosidase; pectin; SusD; SusC; fructokinase; galacto*; arabino*"
        antibiotic_resistance =  "Genes present in the ARDB; multidrug; “azole resistance”; antibiotic resistance”; TetR; “tetracycline resistance”; VanZ; betalactam*; beta-lactam; antimicrob*; lantibio*"
        hypothetical_protein = "hypothetical protein"
        
        gene_classification = {}
        self.pattern_dict["Transposon"] = get(Transposon)
        self.pattern_dict["Plasmid"] = get(Plasmid)
        self.pattern_dict["Phage"] = get(Phage)
        self.pattern_dict["Other"] = get(Other_HGT_machinery)
        self.pattern_dict["CAZYmes"] = get(Carbohydrate_active)
        self.pattern_dict["Antibiotic resistance"] = get(antibiotic_resistance)
        self.pattern_dict["hypothetical protein"] = get(hypothetical_protein)
    

    def classify_product(self, product):
        product_classification = '-'
        for cla in self.pattern_dict:
            # print (cla, pattern_dict[cla])
            pattern = re.compile(self.pattern_dict[cla]) #, re.I
            if pattern.search(product):
                product_classification = cla
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
        for line in open(identified_hgt):
            array = line.strip().split(",")
            if array[1] == "sample":
                continue
            event = Event(array)

            sample = array[1]
            if sample not in self.HGT_event_dict:
                self.HGT_event_dict[sample] = []
            self.HGT_event_dict[sample].append(event) 
    
    def count_times(self):
        ## count how many times the sequence is transferred in each sample
        final_data = []

        data = []
        for sample in self.HGT_event_dict:
            trans_times_dict = {}
            for event in self.HGT_event_dict[sample]:
                recorded_already = False
                for segment_tag in trans_times_dict:
                    array = segment_tag.split("&")
                    if array[0] == event.del_genome and abs(int(event.del_start) - int(array[1])) < self.max_diff and \
                        abs(int(event.del_end) - int(array[2])) < self.max_diff :#and array[3] == event.ins_genome_pure:
                        recorded_already = True
                        trans_times_dict[segment_tag] += 1
                        break
                if not recorded_already:
                    segment_tag = "&".join([event.del_genome, str(event.del_start), str(event.del_end), event.ins_genome_pure])
                    trans_times_dict[segment_tag] = 1
            #print (sample, list(trans_times_dict.values()))
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
            final_data.append([item, counts[item], counts[item]/num, "different"])
        print ("%s of segments are transferred multiple times."%(multiple_trans_num/num), multiple_trans_num, num)  


        data = []
        multi_segment_tag = []
        for sample in self.HGT_event_dict:
            trans_times_dict = {}
            for event in self.HGT_event_dict[sample]:
                recorded_already = False
                for segment_tag in trans_times_dict:
                    array = segment_tag.split("&")
                    if array[0] == event.del_genome and abs(int(event.del_start) - int(array[1])) < self.max_diff and \
                        abs(int(event.del_end) - int(array[2])) < self.max_diff and array[3] == event.ins_genome_pure:
                        recorded_already = True
                        trans_times_dict[segment_tag] += 1
                        break
                if not recorded_already:
                    segment_tag = "&".join([event.del_genome, str(event.del_start), str(event.del_end), event.ins_genome_pure])
                    trans_times_dict[segment_tag] = 1
            #print (sample, list(trans_times_dict.values()))
            for segment_tag in trans_times_dict:
                if trans_times_dict[segment_tag] > 1:  ## the segment exists at least that times in one sample
                    multi_segment_tag.append(segment_tag)
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
            final_data.append([item, counts[item], counts[item]/num, "same"])
        print ("%s of segments are transferred multiple times."%(multiple_trans_num/num), multiple_trans_num, num)   

        df = pd.DataFrame(final_data, columns = ["times", "count", "frequency", "group"])
        df.to_csv("/mnt/d/R_script_files/transfer_times.csv", sep=',')  
        
        self.check_gene_type(multi_segment_tag)

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
    trans = Transfer_times()
    trans.read_events(identified_hgt)
    trans.count_times()

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
        self.min_gene_frac = 0.5
        self.transfer_regions = {}
        self.insert_sites = {}
        # self.no_transfer_regions = {}
        # self.no_ins_regions= {}
        self.transfer_kos = []
        self.no_transfer_kos = []
        self.insert_kos = []
        self.no_insert_kos = []
    
    def classify_regions(self):
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:
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

    def classify_kos(self):
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


if __name__ == "__main__":

    # identified_hgt = "/mnt/d/HGT/seq_ana/identified_event.csv"
    identified_hgt = "/mnt/d/HGT/time_lines/SRP366030.identified_event.csv"
    gff = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.gff"

    transfer_ko_file = "/mnt/d/HGT/seq_ana/transfer_ko.txt"
    no_transfer_ko_file = "/mnt/d/HGT/seq_ana/no_transfer_ko.txt"

    annotation = Annotation(gff)
    annotation.read_gff()

    trans = Transfer_times()
    trans.read_events(identified_hgt)

    extract = Extract_KO(trans.HGT_event_dict)
    extract.classify_regions()
    extract.classify_kos()







        
