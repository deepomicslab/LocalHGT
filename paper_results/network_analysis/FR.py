import copy
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# check common
def check_common(list1, list2):
    diff = list(set(list1).difference(set(list2)))
    if len(diff) == 0:
        return True, []
    else:
        return False, diff

# compute fr
def fr_df(profile, distance_df):
    in_ref, diff = check_common(profile.columns, distance_df.columns)
    if not in_ref:
        print("These species cannot be found in GCN Reference: {}".format(diff))
        exit(0)
    d = copy.deepcopy(distance_df[list(profile.columns)].loc[list(profile.columns)])
    sp_list = list(d.columns)
    
    distance_matrix = d.values
    n_sp = len(sp_list)
    fr_matrix = np.zeros(shape=(n_sp, n_sp))
    max_fr = -np.Inf
    min_fr = np.Inf
    for i in range(n_sp):
        for j in range(i+1, n_sp):
            sp1 = sp_list[i]
            sp2 = sp_list[j]
            x = np.array(profile[sp1])
            y = np.array(profile[sp2])
            FR = np.dot(x, y)*(1-distance_matrix[i][j])
            fr_matrix[i][j] = FR
            # update max and min
            if FR > max_fr:
                max_fr = FR
            if FR < min_fr and FR > 0:
                min_fr = FR   
    log_max = np.log10(max_fr)
    log_min = np.log10(min_fr)
    fr_matrix = fr_normalize(fr_matrix, log_max, log_min)
    fr_matrix += fr_matrix.T
    row, col = np.diag_indices_from(fr_matrix)
    fr_matrix[row, col] = 1
    fr_df = pd.DataFrame(fr_matrix, columns=sp_list, index=sp_list)
    return fr_df, log_max, log_min

def fr_df_without_log(profile, distance_df):
    in_ref, diff = check_common(profile.columns, distance_df.columns)
    if not in_ref:
        print("These species cannot be found in GCN Reference: {}".format(diff))
        exit(0)

    d = copy.deepcopy(distance_df[list(profile.columns)].loc[list(profile.columns)])
    sp_list = list(d.columns)
    
    distance_matrix = d.values
    n_sp = len(sp_list)
    fr_matrix = np.zeros(shape=(n_sp, n_sp))
    max_fr = -np.Inf
    min_fr = np.Inf
    for i in range(n_sp):
        for j in range(i+1, n_sp):
            sp1 = sp_list[i]
            sp2 = sp_list[j]
            x = np.array(profile[sp1])
            y = np.array(profile[sp2])
            FR = np.dot(x, y)*(1-distance_matrix[i][j])
            fr_matrix[i][j] = FR
            # update max and min
            if FR > max_fr:
                max_fr = FR
            if FR < min_fr:
                min_fr = FR   
    log_max = max_fr
    log_min = min_fr
    fr_matrix = fr_normalize_without_log(fr_matrix, log_max, log_min)
    fr_matrix += fr_matrix.T
    row, col = np.diag_indices_from(fr_matrix)
    fr_matrix[row, col] = 1
    fr_df = pd.DataFrame(fr_matrix, columns=sp_list, index=sp_list)
    return fr_df, log_max, log_min

def fr_normalize_without_log(fr_matrix, log_max, log_min):
    range_interval = log_max - log_min
    #fr_matrix = np.array(list(map(np.log10, fr_matrix)))
    fr_matrix -= log_min
    fr_matrix /= range_interval
    fr_matrix[fr_matrix == -np.inf] = 0
    return fr_matrix

def ori_fr(profile, distance_df):
    in_ref, diff = check_common(profile.columns, distance_df.columns)
    if not in_ref:
        print("These species cannot be found in GCN Reference: {}".format(diff))
        exit(0)
    d = copy.deepcopy(distance_df[list(profile.columns)].loc[list(profile.columns)])
    sp_list = list(d.columns)
    
    distance_matrix = d.values
    n_sp = len(sp_list)
    fr_matrix = np.zeros(shape=(n_sp, n_sp))
    for i in range(n_sp):
        for j in range(i+1, n_sp):
            sp1 = sp_list[i]
            sp2 = sp_list[j]
            x = np.array(profile[sp1])
            y = np.array(profile[sp2])
            FR = np.dot(x, y)*(1-distance_matrix[i][j])
            fr_matrix[i][j] = FR
    fr_matrix += fr_matrix.T
    row, col = np.diag_indices_from(fr_matrix)
    fr_matrix[row, col] = 1
    fr_df = pd.DataFrame(fr_matrix, columns=sp_list, index=sp_list)
    return fr_df

# normalize by log10 and rescale to [0,1]
def fr_normalize(fr_matrix, log_max, log_min):
    range_interval = log_max - log_min
    fr_matrix = np.array(list(map(np.log10, fr_matrix)))
    fr_matrix -= log_min
    fr_matrix /= range_interval
    fr_matrix[fr_matrix == -np.inf] = 0
    return fr_matrix




