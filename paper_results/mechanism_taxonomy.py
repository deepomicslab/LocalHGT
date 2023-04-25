import os, pickle


class Taxonomy():
    def __init__(self):
        self.taxonomy_dict = {}
        self.get()
        
    def read_UHGG(self):
        df = pd.read_table(UHGG_meta) 
        for index, row in df.iterrows():
            # self.read_UHGG[row["Genome"]] = row["Species_rep"]
            genome = row["Genome"]
            lineage = row["Lineage"]
            self.taxonomy_dict[genome] = lineage
        # for line in open(UHGG_meta):
        #     array = line.split()

    def get(self):
        save_file = "/mnt/d/breakpoints/script/analysis/taxonomy_dict.pkl"
        if not os.path.isfile(save_file):
            self.read_UHGG()
            with open(save_file, 'wb') as f:
                pickle.dump(self.taxonomy_dict, f)
        else:
            with open(save_file, 'rb') as f:
                self.taxonomy_dict = pickle.load(f)


def calculate_frequency(my_list):
    # create an empty dictionary to store the frequency counts
    freq_counts = {}

    # iterate over each element in the list
    for element in my_list:
        # if the element is already in the dictionary, increment its count
        if element in freq_counts:
            freq_counts[element] += 1
        # otherwise, add it to the dictionary with a count of 1
        else:
            freq_counts[element] = 1
    for element in freq_counts:
        freq_counts[element] = round(freq_counts[element]/len(my_list), 2)
    # return the frequency counts
    return freq_counts

def read_mechanism_result():
    record_dict = {}
    for line in open(mechanism_result):
        array = line.split()
        if array[0] != "event":
            continue
        del_genome = array[2]
        pure_del_genome = "_".join(del_genome.split("_")[:-1])
        mechanism = array[6] 
        lineage = taxonomy.taxonomy_dict[pure_del_genome]

        taxa = lineage.split(";")[2]
        num = mecha2num[mechanism]
        # print (lineage, mechanism)

        if taxa not in record_dict:
            record_dict[taxa] = []
        record_dict[taxa].append(mechanism)
    # print (record_dict)
    for taxa in record_dict:
        freq_counts = calculate_frequency(record_dict[taxa])
        print (taxa, freq_counts)


def read_mechanism_each_sample():
    record_dict = {}
    a, b, c = [], [], [] #"NHEJ":0, "alt-EJ":0, "TEI":0
    for line in open(mechanism_result):
        array = line.strip().split()
        if array[0] != "sample_freq":
            continue
        a.append(float(array[2]))
        b.append(float(array[3]))
        c.append(float(array[4]))
    print (sorted(a)[0], sorted(a)[-1])
    print (sorted(b)[0], sorted(b)[-1])
    print (sorted(c)[0], sorted(c)[-1])
        




if __name__ == "__main__":

    mecha2num = {"NHEJ":0, "alt-EJ":1, "TEI":2}
    UHGG_meta = "/mnt/d/breakpoints/HGT/UHGG/genomes-all_metadata.tsv"
    mechanism_result = "/mnt/d/HGT/time_lines/mechanism_result.txt"

    # taxonomy = Taxonomy()
    # read_mechanism_result()
    read_mechanism_each_sample()