import csv
import pysam
import random
from pyfaidx import Fasta
import skbio
from skbio import DNA, TabularMSA
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import StripedSmithWaterman
import numpy as np
import argparse
import sys
import re
import pandas as pd
import os
import pickle
from scipy.stats import mannwhitneyu
from scipy import stats
from scipy import linalg
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums
import networkx as nx
import math
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA
from random import shuffle
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.feature_selection import SelectFromModel
from scipy.sparse import csgraph
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.inspection import permutation_importance
from scipy.stats import wilcoxon
from sklearn.preprocessing import normalize
import scipy
from scipy.linalg import svd
import powerlaw
np.set_printoptions(threshold=sys.maxsize)
import scipy.special as sc
from scipy.special import comb


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


def get_split_reads_cutoff(g): # g is the number of reads in the specific sample
    # prior_a = 25 # prior probability
    prior_b = 50000000#63333330  # prior probability
    given_n = 42648185 # mean number of reads among samples
    # # n = 182534663 # max number of reads 
    # given_r = 2 # the cutoff with n reads  4
    prior_a = 20 # prior probability
    given_r = 2 # the cutoff with n reads  4

    for m in range(100):
        alpha = prior_a + m
        beta= prior_b + g - m
        p = 1
        for k in range(given_r + 1):
            choose = comb(given_n, k)
            beta1 = sc.beta(alpha + k, beta + given_n - k)
            beta2 = sc.beta(alpha, beta)
            # print (choose, beta1)
            p -= choose * beta1 / beta2
        if p > 0.9: #0.9
            break
    return m, p

def read_bkp(bkp_file, filtered_file):
    f = open(bkp_file)
    all_rows = csv.reader(f)

    out = open(filtered_file, 'w')
    csv_writer = csv.writer(out)
    final_row_num = 0

    for row in all_rows:
        if row[0][0] == "#":
            reads_num = int(row[0].split(";")[0].split(":")[1])
            min_split_num, p = get_split_reads_cutoff(reads_num)
            # print (bkp_file, reads_num, "cutoff", min_split_num)
        elif row[0] == "from_ref":
            pass
        else:
            # print (row)
            if len(row) < 13:
                continue
            eb = Acc_Bkp(row)
            if eb.cross_split_reads < min_split_num:
                continue
            final_row_num += 1
        csv_writer.writerow(row)

    f.close()
    out.close()
    if final_row_num == 0:
        os.system("rm %s"%(filtered_file))

def main(result_dir, sample_num):
    files = os.listdir(result_dir)
    # print (len(files))
    for acc_file in files:
        if not re.search("acc.csv", acc_file):
            continue
        if re.search("repeat.acc.csv", acc_file):
            continue
        # if hgt_result_dir + acc_file != "/mnt/d/breakpoints/script/analysis/hgt_results/CCIS98832363ST-4-0.acc.csv":
        #     continue
        # print (hgt_result_dir + acc_file)
        read_bkp(result_dir + acc_file, filter_hgt_result_dir + acc_file)
        sample_num += 1
        # break
    return sample_num


if __name__ == "__main__":
    hgt_result_dir = "/mnt/d/breakpoints/script/analysis/hgt_results/"
    tgs_dir = "/mnt/d/HGT/time_lines/SRP366030/"
    wenkui_dir = "/mnt/d/breakpoints/HGT/CRC/wenkui/"

    filter_hgt_result_dir = "/mnt/d/breakpoints/script/analysis/filter_hgt_results/"
    sample_num = 0
    sample_num = main(tgs_dir, sample_num)
    sample_num = main(wenkui_dir, sample_num)
    sample_num = main(hgt_result_dir, sample_num)


    print (sample_num)