import csv

class Acc_Bkp(object):
    def __init__(self, list):
        self.from_ref = list[0]
        self.to_ref = list[2]
        self.from_bkp = int(list[1])
        self.to_bkp = int(list[3])
        self.if_reverse = list[6]
        self.from_side = list[4]
        self.to_side = list[5]
        self.from_ref_genome = "_".join(self.from_ref.split("_")[:-1])
        # self.from_ref_lineage = taxonomy.taxonomy_dict[self.from_ref_genome]
        self.to_ref_genome = "_".join(self.to_ref.split("_")[:-1])
        # self.to_ref_lineage = taxonomy.taxonomy_dict[self.to_ref_genome]
        self.score = float(list[9])

class Find_trans_gene():

    def __init__(self, bkp_file):
        self.bkp_file = bkp_file
        self.bkps = []
        self.read_bkp()
         
    def read_bkp(self):
        f = open(self.bkp_file)
        all_rows = csv.reader(f)
        for row in all_rows:
            if row[0] == "from_ref":
                continue
            eb = Acc_Bkp(row)
            if eb.from_ref_genome != eb.to_ref_genome:
                self.bkps.append(eb)
        f.close()

    def test(self):
        genome_pair = {}
        for bkp in self.bkps:
            refs = [bkp.from_ref, bkp.to_ref]
            sorted_refs = sorted(refs)
            pair_name = "&".join(sorted_refs)
            if sorted_refs[0] == refs[0]:
                position = [bkp.from_bkp, bkp.to_bkp]
            else:
                position = [bkp.to_bkp, bkp.from_bkp]
            if pair_name not in genome_pair:
                genome_pair[pair_name] = []
            genome_pair[pair_name].append(position)
        for pair_name in genome_pair:
            if len(genome_pair[pair_name]) == 2:
                print (pair_name, genome_pair[pair_name])


bkp_file = "/mnt/d/breakpoints/script/analysis/new_result/ERR2726420.acc.csv"
fts = Find_trans_gene(bkp_file)
fts.test()

