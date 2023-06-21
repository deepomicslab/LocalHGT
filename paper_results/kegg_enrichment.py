from Bio.KEGG import REST
from collections import defaultdict
from scipy.stats import fisher_exact
import pandas as pd

# Define the input and background KO IDs

def read_list(file_path):
    with open(file_path, "r") as f:
        return [line.strip() for line in f.readlines()]

def get_pathways(input_ko_ids):
    # Retrieve the KEGG pathways for the input and background KO IDs
    input_pathways = defaultdict(set)
    for ko_id in input_ko_ids:
        pathway_file = REST.kegg_link('pathway', ko_id).read()
        for line in pathway_file.rstrip().split('\n'):
            if line.strip() == "":
                continue
            fields = line.split('\t')
            #print (fields)
            pathway_id = fields[1].split(':')[-1]
            input_pathways[ko_id].add(pathway_id)

    # Count the number of input and background KO IDs in each pathway
    input_counts = defaultdict(int)
    for ko_id in input_ko_ids:
        for pathway_id in input_pathways[ko_id]:
            input_counts[pathway_id] += 1

    return input_pathways, input_counts

input_ko_ids = read_list("/mnt/d/HGT/seq_ana/transfer_ko.txt")
background_ko_ids = read_list("/mnt/d/HGT/seq_ana/no_transfer_ko.txt")
input_pathways, input_counts = get_pathways(input_ko_ids)
background_pathways, background_counts = get_pathways(background_ko_ids)

print (input_ko_ids[:5])

data = []
# Perform a Fisher's exact test for each pathway
for pathway_id in set(input_counts.keys()) | set(background_counts.keys()):
    a = input_counts[pathway_id]
    b = len(input_ko_ids) - a
    c = background_counts[pathway_id]
    d = len(background_ko_ids) - c
    p_value = fisher_exact([[a, b], [c, d]])[1]

    pathway_name = REST.kegg_get(f'path:{pathway_id}').read().split('\n')[1].split(';')[0].strip()
    pathway_name = pathway_name.replace("NAME", '')
    pathway_name = pathway_name.replace('"', '')
    pathway_name = pathway_name.strip()

    input_ratio = float(a)/len(input_ko_ids)
    background_ratio = float(c)/len(background_ko_ids)

    if p_value < 0.05 and input_ratio > background_ratio:
        print(f'{pathway_name} ({pathway_id}): enriched, p={p_value:.3f}, input={a}/{len(input_ko_ids)}, background={c}/{len(background_ko_ids)}')

        data.append([pathway_name, pathway_id, p_value, a])

df = pd.DataFrame(data, columns = ["pathway_name", "pathway_id", "p_value", "KO_num"])
df.to_csv("/mnt/d/HGT/seq_ana/transfer_enriched.csv", sep=',')
