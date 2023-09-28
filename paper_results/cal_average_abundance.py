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
# from math import comb #The comb function is new in Python 3.8
from scipy.special import comb
# from deep_learning import Classifier
# from KR_norm_juicer import KRnorm_sym
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from collections import defaultdict

cohort_abd = {"YuJ_2015":"2021-03-31.YuJ_2015.relative_abundance.xls",
"WirbelJ_2018":"2021-03-31.WirbelJ_2018.relative_abundance.xls",
"HanniganGD_2017":"2021-03-31.HanniganGD_2017.relative_abundance.xls",
"YachidaS_2019":"2021-10-14.YachidaS_2019.relative_abundance.xls",
"ThomasAM_2018a":"2021-03-31.ThomasAM_2018a.relative_abundance.xls",
"ThomasAM_2018b":"2021-03-31.ThomasAM_2018b.relative_abundance.xls",
"ZellerG_2014":"2021-03-31.ZellerG_2014.relative_abundance.xls",
"FengQ_2015":"2021-03-31.FengQ_2015.relative_abundance.xls",
"VogtmannE_2016":"2021-03-31.VogtmannE_2016.relative_abundance.xls",
"KarlssonFH_2013":"2021-10-14.KarlssonFH_2013.relative_abundance.xls",
"NielsenHB_2014":"2021-03-31.NielsenHB_2014.relative_abundance.xls",
"HallAB_2017":"2021-10-14.HallAB_2017.relative_abundance.xls",
"QinJ_2012":"2021-10-14.QinJ_2012.relative_abundance.xls",
"DavidLA_2015":"2021-03-31.DavidLA_2015.relative_abundance.xls",
"KieserS_2018":"2021-10-14.KieserS_2018.relative_abundance.xls"
}

level = 6

def each_cohort_get_taxa(abd_file_path):
    taxa_set = set()
    f = open(abd_file_path, 'r')
    i = 0
    for line in f:
        array = line.strip().split()
        if i > 0:
            species_name = array[0].split("|")[-level]
            taxa_set.add(species_name)
            # print (species_name)
        i += 1
    f.close()
    return taxa_set

def each_cohort(abd_file_path):
    taxa_set = set()
    f = open(abd_file_path, 'r')
    sample_abd = {}
    i = 0
    for line in f:
        array = line.strip().split()
        if i == 0:
            sample_list = array
            for sample in sample_list:
                sample_abd[sample] = defaultdict(float)
        else:
            species_name = array[0].split("|")[-level]
            genus_name = array[0].split("|")[-2]

            taxa_set.add(species_name)

            for j in range(len(sample_list)):
                sample = sample_list[j]
                abundance = float(array[j+1])
                sample_abd[sample][species_name] += abundance
        i += 1
    f.close()
    return sample_abd, taxa_set

    # elif genus_name in marker_genus:
    #     # print ("name", species_name, marker_genus)
    #     species_index = marker_genus[genus_name]
    #     for j in range(len(sample_list)):
    #         sample = sample_list[j]
    #         abundance = float(array[j+1])
    #         sample_abd[sample][species_index] += abundance


def taxa_average_abun():
    all_taxa_set = set()

    for cohort in cohort_abd:
        abd_file = cohort_abd[cohort]
        abd_file_path = "/mnt/d/breakpoints/script/analysis/use/" + abd_file
        taxa_set = each_cohort_get_taxa(abd_file_path)
        # print (taxa_set)
        all_taxa_set = all_taxa_set|taxa_set
    print (len(all_taxa_set))

    all_sample_abd = {}
    for taxon in all_taxa_set:
        all_sample_abd[taxon] = []

    for cohort in cohort_abd:
        abd_file = cohort_abd[cohort]
        abd_file_path = "/mnt/d/breakpoints/script/analysis/use/" + abd_file
        sample_abd, taxa_set = each_cohort(abd_file_path)
        for sample in sample_abd:
            for taxon in all_sample_abd:
                if taxon not in sample_abd[sample]:
                    all_sample_abd[taxon].append(0)
                else:
                    all_sample_abd[taxon].append(sample_abd[sample][taxon])
    
    ave_abun_dict = {}
    for taxon in all_sample_abd:
        ave_abun_dict[taxon] = np.mean(all_sample_abd[taxon])
        # print (taxon, np.mean(all_sample_abd[taxon]))
    sorted_dict = dict(sorted(ave_abun_dict.items(), key=lambda x: x[1], reverse = False))
    # for taxa in sorted_dict:
    #     print (taxa, sorted_dict[taxa])
    return sorted_dict

    # return all_sample_abd


