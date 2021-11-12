#!/usr/bin/env python3

"""
compare the output hgt and ground truth,
to evaluate our tool.
"""

import os
import numpy as np
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import random
import sys
from collections import Counter
from datetime import datetime
from simulation import Parameters

tolerate_dist = 50

def read_interval(interval_file, true_locus):
    gap = 100
    cover_flag = False
    f = open(interval_file)
    for line in f:
        line = line.strip()
        if line == '':
            continue
        chr_name = line.split(":")[0]
        start = int(line.split(":")[1].split("-")[0])
        end = int(line.split(":")[1].split("-")[1])
        if chr_name == true_locus[0] and int(true_locus[1]) > start + gap and int(true_locus[1]) < end - gap:
            cover_flag = True
        #     print (true_locus, array)
    f.close()
    return cover_flag

def check_if_bkp_in_extracted_ref(true, interval_file):
    true_bkp = read_true(true)
    i = 0
    for true_locus in true_bkp:
        if read_interval(interval_file, true_locus):
            i += 1
            # print ('Have ref:', true_locus)
        else:
            print ('Lack ref:', true_locus)
    ref_len = 0
    f = open(interval_file)
    for line in f:
        array = line.strip().split()
        start = int(line.split(":")[1].split("-")[0])
        end = int(line.split(":")[1].split("-")[1])
        ref_len += abs(end - start ) 
    f.close()
    print (len(true_bkp))
    return float(i)/len(true_bkp), ref_len

def read_true(true):
    true_bkp = []
    for line in open(true):
            array = line.strip().split()
            true_bkp.append([array[0], array[1], array[2], array[3]])
            true_bkp.append([array[0], array[1], array[2], array[4]])
    return true_bkp

def read_lemon(lemon):
    lemon_bkp = []
    past = ['', '', '', '']
    for line in open(lemon):
        array = line.strip().split(',')
        if array[0] == "from_ref": #skip the annotation line
            continue
        if '_'.join(array[:4]) == '_'.join(past):
            continue
        lemon_bkp.append(array[:4])
        past = array[:4]
    return lemon_bkp

def compare(true_bkp, our_bkp):
    right = 0
    error = 0
    correct_result_num = 0
    for true in true_bkp:
        identified = False
        for our in our_bkp:
            if true[0] == our[0] and true[2] == our[2] and abs(int(true[1])-int(our[1]))\
                < tolerate_dist and abs(int(true[3])-int(our[3])) < tolerate_dist:
                right += 1
                identified = True
                break
            elif true[0] == our[2] and true[2] == our[0] and abs(int(true[1])-int(our[3]))\
                < tolerate_dist and abs(int(true[3])-int(our[1])) < tolerate_dist:
                right += 1
                identified = True
                break
        if identified == False:
            print ("Missed bkp:", true)
    print ("-----------")  
    accuracy = right/len(true_bkp)

    #find false positive locus
    false_positive_locus = []
    for true in our_bkp:
        identified = False
        for our in true_bkp:
            if true[0] == our[0] and true[2] == our[2] and abs(int(true[1])-int(our[1]))\
                < tolerate_dist and abs(int(true[3])-int(our[3])) < tolerate_dist:
                    right += 1
                    identified = True
                    break
            elif true[0] == our[2] and true[2] == our[0] and abs(int(true[1])-int(our[3]))\
                < tolerate_dist and abs(int(true[3])-int(our[1])) < tolerate_dist:
                    right += 1
                    identified = True
                    break
        if not identified:
            false_positive_locus.append(true)
            print ("False bkp:", true)
    if len(our_bkp) > 0:
        FDR = len(false_positive_locus)/len(our_bkp)
    else:
        FDR = 0
    return round(accuracy,2), round(FDR,2)#, false_positive_locus

class Performance():
    def __init__(self, accuracy, FDR, tool_time, tool_mem):
        self.accuracy = accuracy
        self.FDR = FDR
        self.user_time = tool_time
        self.max_mem = tool_mem
        self.ref_accuracy = 0
        self.ref_len = 0

    def add_ref(self, ref_accuracy, ref_len):
        self.ref_accuracy = ref_accuracy
        self.ref_len = ref_len

class Sample():
    def __init__(self, ID, true_dir):
        self.ID = ID
        true_ID = self.get_true(ID)
        self.true_file = true_dir + '/' + true_ID + '.true.sv.txt'
        self.true_bkp = read_true(self.true_file)

    def get_true(self, ID):
        array = ID.split('_')
        if len(array) == 6:
            return ID # pure ID
        else:
            return '_'.join(array[:6])

    def eva_tool(self, tool_dir):
        acc_file = tool_dir + '/' + self.ID + '.acc.txt'
        bkp = read_lemon(acc_file)
        accuracy, FDR = compare(self.true_bkp, bkp)
        time_file = tool_dir + '/' + self.ID + '.time'
        tool_time = self.extract_time(time_file)
        tool_mem = self.extract_mem(time_file)
        pe = Performance(accuracy, FDR, tool_time, tool_mem)
        return pe

    def extract_time(self, time_file): #log file obtained by /usr/bin/time -v
         #if no time available
        for line in open(time_file):
            time_re = re.search('User time \(seconds\):(.*?)$', line)
            if time_re:
                user_time =  time_re.group(1).strip()

            time_sys = re.search('System time \(seconds\):(.*?)$', line)
            if time_sys:
                sys_time = time_sys.group(1).strip()
        # print (user_time, sys_time)
        all_time = float(user_time) + float(sys_time)
        final_time = round(all_time/3600, 2)
        return final_time

    def extract_mem(self, time_file):
        used_mem = 0 #if no time available
        for line in open(time_file):
            time_re = re.search('Maximum resident set size \(kbytes\):(.*?)$', line)
            if time_re:
                used_mem =  time_re.group(1).strip()
        final_mem = round(float(used_mem)/1000000, 2)
        return final_mem

    def eva_ref(self, tool_dir):
        interval_file = tool_dir + '/%s.interval.txt.bed'%(self.ID)
        print (interval_file)
        ref_accuracy, ref_len = check_if_bkp_in_extracted_ref(self.true_file, interval_file)
        return ref_accuracy, round(ref_len/1000000, 2) #M

class Figure():
    def __init__(self):
        self.data = []
        self.df = []
    
    def add_local_sample(self, pe, va): # any performance Object
        self.data.append([pe.user_time, pe.accuracy, pe.FDR, \
        pe.max_mem, pe.ref_accuracy, pe.ref_len, "LocalHGT", va])

    def add_lemon_sample(self, pe, va): # any performance Object
        self.data.append([pe.user_time, pe.accuracy, pe.FDR, \
        pe.max_mem, pe.ref_accuracy, pe.ref_len, "LEMON", va])

    def convert_df(self):
        self.df=pd.DataFrame(self.data,columns=['CPU time', 'Sensitivity','FDR', 'Max Memory (G)', \
        'Ref Accuracy',  'Extracted Ref (M)', 'Methods', 'Variation'])

    def plot(self):
        fi.convert_df()
        fig, axes = plt.subplots(3, 2, figsize=(15,10))
        sns.barplot(ax = axes[0][0], x='Variation',y='Sensitivity',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[0][1], x='Variation',y='FDR', hue= 'Methods',data=self.df) 
        sns.barplot(ax = axes[1][0], x='Variation',y='Max Memory (G)',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[1][1], x='Variation',y='CPU time', hue= 'Methods',data=self.df) 
        sns.barplot(ax = axes[2][0], x='Variation',y='Extracted Ref (M)',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[2][1], x='Variation',y='Ref Accuracy', hue= 'Methods',data=self.df)        
        #     plt.xticks(rotation=0)
        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        plt.savefig('/mnt/d/breakpoints/HGT/figures/HGT_comparison_%s.pdf'%(give_time))


if __name__ == "__main__":
    true_dir = "/mnt/d/breakpoints/HGT/uhgg_snp/"
    lemon_dir = "/mnt/d/breakpoints/HGT/lemon_snp/"
    local_dir = "/mnt/d/breakpoints/HGT/uhgg_snp_results/"
    fi = Figure()
    ba = Parameters()

    for snp_rate in ba.snp_level[1:-1]: # 0.01-0.09
        ba.change_snp_rate(snp_rate)
        for index in range(ba.iteration_times):
            ba.get_ID(index)    
            sa = Sample(ba.sample, true_dir)
            print (ba.sample)
            # lemon_pe = sa.eva_tool(lemon_dir)
            local_pe = sa.eva_tool(local_dir)
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
            local_pe.add_ref(ref_accuracy, ref_len)
            # print ("time:", lemon_pe.user_time, local_pe.user_time)
            # print ("sensitivity", lemon_pe.accuracy, local_pe.accuracy)
            # print ("mem",lemon_pe.max_mem, local_pe.max_mem)
            print ("ref", local_pe.ref_accuracy, local_pe.ref_len, "Mb")
            
            # fi.add_local_sample(local_pe, snp_rate)
            # fi.add_lemon_sample(lemon_pe, snp_rate)
            # break
        # break
    # fi.plot()
