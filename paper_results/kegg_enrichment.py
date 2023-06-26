from Bio.KEGG import REST
from collections import defaultdict
from scipy.stats import fisher_exact
import pandas as pd
import pickle
import re
from statsmodels.stats.multitest import multipletests

# Define the input and background KO IDs

def read_list(file_path):
    with open(file_path, "r") as f:
        return [line.strip() for line in f.readlines()]

def get_pathways(input_ko_ids):
    # Count the number of input and background KO IDs in each pathway
    input_counts = defaultdict(int)
    for ko_id in input_ko_ids:
        for pathway_id in ko_pathway_dict[ko_id]:
            input_counts[pathway_id] += 1
    return input_counts

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

def enrichment_analysis(input_counts, background_counts):
    data = []
    # Perform a Fisher's exact test for each pathway
    for pathway_id in set(input_counts.keys()) | set(background_counts.keys()):
        a = input_counts[pathway_id]
        b = len(input_ko_ids) - a
        c = background_counts[pathway_id]
        d = len(background_ko_ids) - c
        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
        
        if p_value >= 0.01:
            continue

        pathway_name, first_class, second_class = get_pathway_name_class(pathway_id)

        if oddsratio > 1:
            print(f'{pathway_name} ({first_class}): enriched, p={p_value:.3f}, input={a}/{len(input_ko_ids)}, background={c}/{len(background_ko_ids)}, fold_enrichment={oddsratio:.3f}')
        else:
            print(f'{pathway_name} ({first_class}): depleted, p={p_value:.3f}, input={a}/{len(input_ko_ids)}, background={c}/{len(background_ko_ids)}, fold_enrichment={oddsratio:.3f}')

        data.append([pathway_name, pathway_id, p_value, a, oddsratio, first_class, second_class])

    df = pd.DataFrame(data, columns = ["pathway_name", "pathway_id", "p_value", "gene_num", "fold", "first_class", "second_class"])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected
    df.to_csv(output, sep=',')
    print ("enriched pathway num", len(data))


if __name__ == "__main__":

    with open('/mnt/d/HGT/seq_ana/ko_pathway_dict.pickle', 'rb') as f:
        ko_pathway_dict = pickle.load(f)


    input_ko_ids = read_list("/mnt/d/HGT/seq_ana/transfer_ko.txt")
    background_ko_ids = read_list("/mnt/d/HGT/seq_ana/no_transfer_ko.txt")
    output = "/mnt/d/R_script_files/transfer_enriched.csv"

    # input_ko_ids = read_list("/mnt/d/HGT/seq_ana/insert_ko.txt")
    # background_ko_ids = read_list("/mnt/d/HGT/seq_ana/no_insert_ko.txt")
    # output = "/mnt/d/R_script_files/insert_enriched.csv"

    # input_ko_ids = read_list("/mnt/d/HGT/seq_ana/bkp_ko.txt")
    # background_ko_ids = read_list("/mnt/d/HGT/seq_ana/no_bkp_ko.txt")
    # output = "/mnt/d/R_script_files/bkp_enriched.csv"

    # background_ko_ids = read_list("/mnt/d/HGT/seq_ana/all_kos.txt")

    # print (input_ko_ids[:5])
    input_counts = get_pathways(input_ko_ids)
    background_counts = get_pathways(background_ko_ids)
    enrichment_analysis(input_counts, background_counts)



