"""
collect all the ko ids from the UHGG database
obtain the corresponding pathways of each ko
save the (ko:pathway_list) relationship in the dict

input : gff file
output: dict (ko:pathway_list)

"""


from Bio.KEGG import REST
from collections import defaultdict
from scipy.stats import fisher_exact
import pandas as pd
import pickle

def extract_ko_gff(gff_file):
    ko_list = set()
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):  # Skip comment lines
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                attributes = fields[8].split(';')
                for attr in attributes:
                    if attr.startswith('ID='):
                        gene_id = attr[3:]
                    elif attr.startswith('KEGG='):
                        if attr == "KEGG=-":
                            continue
                        ko_hub = attr[5:]
                        new_arr = ko_hub.split(",")
                        for ko_id in new_arr:
                            ko_list.add(ko_id[3:])
                        # print(gene_id, new_arr)
    return list(ko_list)

def get_pathways(input_ko_ids):
    # Retrieve the KEGG pathways for the input and background KO IDs
    input_pathways = defaultdict(set)
    i = 0
    for ko_id in input_ko_ids:
        if i % 500 == 0:
            print (i)
        pathway_file = REST.kegg_link('pathway', ko_id).read()
        for line in pathway_file.rstrip().split('\n'):
            if line.strip() == "":
                continue
            fields = line.split('\t')
            #print (fields)
            pathway_id = fields[1].split(':')[-1]
            input_pathways[ko_id].add(pathway_id)
        i += 1
    return input_pathways


gff = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.gff"
UHGG_uniq_ko_list = extract_ko_gff(gff)
print("UHGG KO num:", len(UHGG_uniq_ko_list))

input_pathways = get_pathways(UHGG_uniq_ko_list)

# Open a file for writing in binary mode
with open('/mnt/d/HGT/seq_ana/ko_pathway_dict.pickle', 'wb') as f:
    # Use the pickle module to serialize the dictionary and write it to the file
    pickle.dump(input_pathways, f)
