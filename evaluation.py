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
ref_gap = 50

class Read_bed(object):
    def __init__(self):
        self.true_interval = {}
        self.gap = ref_gap
        self.ref_len = 0

    def read(self, interval_file):
        f = open(interval_file)
        for line in f:
            line = line.strip()
            if line == '':
                continue
            chr_name = line.split(":")[0]
            start = int(line.split(":")[1].split("-")[0])
            end = int(line.split(":")[1].split("-")[1])   
            self.ref_len += (end - start)
            self.add(chr_name, start, end)  
        f.close()   

    def add(self, chr_name, start, end):
        if chr_name in self.true_interval:
            self.true_interval[chr_name].append([start, end])
        else:
            self.true_interval[chr_name] = [[start, end]]

    def search(self, true_locus):
        if true_locus[0] not in self.true_interval:
            return False
        else:
            for true in self.true_interval[true_locus[0]]:
                start = true[0]
                end = true[1]
                if int(true_locus[1]) > start + self.gap and int(true_locus[1]) < end - self.gap:
                    return True
            return False

def check_if_bkp_in_extracted_ref(true, interval_file):
    bed_obj = Read_bed()
    bed_obj.read(interval_file)

    all_pos = read_all_frag(true)
    i = 0
    for single_locus in all_pos:
        if bed_obj.search(single_locus):
            i += 1
            # print ('Have ref:', true_locus)
        else:
            print ('Lack ref:', single_locus)
    return round(float(i)/len(all_pos),2), bed_obj.ref_len

def read_all_frag(true):
    all_pos = []
    for line in open(true):
            array = line.strip().split()
            all_pos.append([array[0], array[1]])
            all_pos.append([array[2], array[3]])
            all_pos.append([array[2], array[4]])

    return all_pos

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

def read_localHGT(lemon):
    lemon_bkp = []
    past = ['', '', '', '']
    for line in open(lemon):
        if line[0] == "#":
            continue
        array = line.strip().split(',')
        if array[0] == "from_ref": #skip the annotation line
            continue
        from_ref = array[0]
        from_pos = int(array[1])
        to_ref = array[4]
        to_pos = int(array[5])
        lemon_bkp.append([from_ref, from_pos, to_ref, to_pos])
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
    accuracy = right/len(true_bkp)
    recall = accuracy

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
            # print ("False bkp:", true)
    if len(our_bkp) > 0:
        FDR = len(false_positive_locus)/len(our_bkp)
    else:
        FDR = 0
    precision = 1-FDR
    if precision > 0 and recall > 0:
        F1_score = 2/((1/precision) + (1/recall)) 
    else:
        F1_score = 0
    return round(accuracy,2), round(FDR,2), round(F1_score,2)#, false_positive_locus

class Performance():
    def __init__(self, accuracy, FDR, tool_time, tool_mem, F1_score, complexity):
        self.accuracy = accuracy
        self.FDR = FDR
        self.user_time = tool_time
        self.max_mem = tool_mem
        self.ref_accuracy = 0
        self.ref_len = 0
        self.F1_score = F1_score
        self.complexity = complexity

    def add_ref(self, ref_accuracy, ref_len):
        self.ref_accuracy = ref_accuracy
        self.ref_len = ref_len

class Sample():
    def __init__(self, ID, true_dir):
        self.ID = ID
        true_ID = self.get_true(ID)
        self.true_file = true_dir + '/' + true_ID + '.true.sv.txt'
        self.true_bkp = read_true(self.true_file)
        self.complexity = ''

    def get_true(self, ID):
        array = ID.split('_')
        if len(array) == 6:
            return ID # pure ID
        else:
            return '_'.join(array[:6])

    def eva_tool(self, tool_dir, tool):
        acc_file = tool_dir + '/' + self.ID + '.acc.csv'
        if tool == "lemon":
            bkp = read_lemon(acc_file)
        else:
            bkp = read_localHGT(acc_file)
        accuracy, FDR, F1_score = compare(self.true_bkp, bkp)
        time_file = tool_dir + '/' + self.ID + '.time'
        tool_time = self.extract_time(time_file)
        tool_mem = self.extract_mem(time_file)
        pe = Performance(accuracy, FDR, tool_time, tool_mem, F1_score, self.complexity)
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
        # print (interval_file)
        ref_accuracy, ref_len = check_if_bkp_in_extracted_ref(self.true_file, interval_file)
        return ref_accuracy, round(ref_len/1000000, 2) #M

    def change_ID(self, new_id):
        self.ID = new_id

class Figure():
    def __init__(self):
        self.data = []
        self.df = []
        self.variation = "Depth"
    
    def add_local_sample(self, pe, va): # any performance Object
        self.data.append([pe.user_time, pe.accuracy, pe.FDR, \
        pe.max_mem, pe.ref_accuracy, pe.ref_len, "LocalHGT", va, pe.F1_score, pe.complexity])

    def add_lemon_sample(self, pe, va): # any performance Object
        self.data.append([pe.user_time, pe.accuracy, pe.FDR, \
        pe.max_mem, pe.ref_accuracy, pe.ref_len, "LEMON", va, pe.F1_score, pe.complexity])

    def convert_df(self):
        self.df=pd.DataFrame(self.data,columns=['CPU time', 'Recall','FDR', 'Peak RAM', \
        'Ref Accuracy',  'Extracted Ref (M)', 'Methods', self.variation, "F1", "Complexity"])

    def plot(self):
        self.convert_df()
        # fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        # sns.boxplot(ax = axes[0], x=self.variation,y='Recall',hue= 'Methods',data=self.df)
        # sns.barplot(ax = axes[1], x=self.variation,y='FDR', hue= 'Methods',data=self.df) 
        # axes[1].set_ylim(0,0.05)   
        ax = sns.barplot(x=self.variation, y="F1 score",hue= 'Methods',data=self.df)   
            # plt.xticks(rotation=0)
        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        plt.savefig('/mnt/d/breakpoints/HGT/figures/HGT_comparison_%s.pdf'%(give_time))
        self.df.to_csv('analysis/depth_comparison.csv', sep=',')

    def plot_all(self):
        self.convert_df()
        fig, axes = plt.subplots(3, 2, figsize=(15,10))
        sns.barplot(ax = axes[0][0], x=self.variation,y='Recall',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[0][1], x=self.variation,y='FDR', hue= 'Methods',data=self.df) 
        axes[0,1].set_ylim(0,0.05)
        sns.barplot(ax = axes[1][0], x=self.variation,y='Peak RAM',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[1][1], x=self.variation,y='CPU time', hue= 'Methods',data=self.df) 
        sns.barplot(ax = axes[2][0], x=self.variation,y='Extracted Ref (M)',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[2][1], x=self.variation,y='Ref Accuracy', hue= 'Methods',data=self.df)        
        #     plt.xticks(rotation=0)
        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        plt.savefig('/mnt/d/breakpoints/HGT/figures/HGT_comparison_%s.pdf'%(give_time))
        
    

    def plot_cami(self):
        self.convert_df()
        print (self.df)
        fig, axes = plt.subplots(3, 3, figsize=(15,10))
        sns.barplot(ax = axes[0][0], x=self.variation, y='CPU time', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'low']).set_title('Low')  
        sns.barplot(ax = axes[0][1], x=self.variation, y='CPU time', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'medium']).set_title('Medium') 
        sns.barplot(ax = axes[0][2], x=self.variation, y='CPU time', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'high']).set_title('High')

        sns.barplot(ax = axes[1][0], x=self.variation, y='Peak RAM', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'low'])  
        sns.barplot(ax = axes[1][1], x=self.variation, y='Peak RAM', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'medium']) 
        sns.barplot(ax = axes[1][2], x=self.variation, y='Peak RAM', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'high']) 

        sns.barplot(ax = axes[2][0], x=self.variation, y='Recall', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'low'])  
        sns.barplot(ax = axes[2][1], x=self.variation, y='Recall', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'medium']) 
        sns.barplot(ax = axes[2][2], x=self.variation, y='Recall', hue= 'Methods',data=self.df.loc[self.df['Complexity'] == 'high']) 

  
        # sns.lineplot(ax = axes[1], x=self.variation, y='CPU time', style="Complexity", \
        # hue= 'Methods',data=self.df, markers=True, dashes=False) 
        # sns.lineplot(ax = axes[2], x=self.variation, y='Peak RAM', style="Complexity", \
        # hue= 'Methods',data=self.df, markers=True, dashes=False) 
        
        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        plt.savefig('/mnt/d/breakpoints/HGT/figures/HGT_comparison_%s.pdf'%(give_time))
        self.df.to_csv('analysis/cami_comparison.csv', sep=',')

def cami():
    fi = Figure()
    ba = Parameters()
    ba.get_dir(true_dir)
    fi.variation = "Mutation Rate"
    for snp_rate in [0.01, 0.02, 0.03, 0.04, 0.05]:
    # for snp_rate in [0.05]:
        ba.change_snp_rate(snp_rate)
        index = 0
        ba.get_ID(index)
        sa = Sample(ba.sample, true_dir)
        for level in ba.complexity_level:
            cami_ID = ba.sample + '_' + level
            sa.change_ID(cami_ID)
            sa.complexity = level
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
            local_pe = sa.eva_tool(local_dir) 
            local_pe.add_ref(ref_accuracy, ref_len)
            
            fi.add_local_sample(local_pe, snp_rate)
            lemon_pe = sa.eva_tool(lemon_dir)
            fi.add_lemon_sample(lemon_pe, snp_rate)
            print ("#",cami_ID, ref_accuracy, ref_len, "Mb", local_pe.accuracy, local_pe.F1_score, local_pe.complexity)
            print ("############ref" ,ba.sample, ref_accuracy, ref_len, "Mb", local_pe.accuracy ,lemon_pe.accuracy)
    fi.plot_cami()
            
def snp():
    fi = Figure()
    ba = Parameters()

    for snp_rate in ba.snp_level: # 0.01-0.09
        # if snp_rate == 0.07:
        #     continue
    # for snp_rate in [0.05]:
        ba.change_snp_rate(snp_rate)
        for index in range(ba.iteration_times):
            ba.get_ID(index)    
            sa = Sample(ba.sample, true_dir)
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
            
            local_pe = sa.eva_tool(local_dir)
            print ("############ref" ,ba.sample, ref_accuracy, ref_len, "Mb", local_pe.accuracy, local_pe.FDR)
            # print ("############ref" ,ba.sample, ref_accuracy, ref_len, "Mb", local_pe.accuracy)

def depth():
    fi = Figure()
    ba = Parameters()
    true_dir = "/mnt/d/breakpoints/HGT/uhgg_depth/"
    lemon_dir = "/mnt/d/breakpoints/HGT/lemon_depth/"
    local_dir = "/mnt/d/breakpoints/HGT/uhgg_depth_results/"

    for depth in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        ba.change_depth(depth)
        for index in range(10):
            ba.get_ID(index)    
            sa = Sample(ba.sample, true_dir)
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
            
            # print ("############ref" ,ba.sample, ref_accuracy, ref_len, "Mb")
            local_pe = sa.eva_tool(local_dir)
            local_pe.add_ref(ref_accuracy, ref_len)
            print ("--------next lemon-------")
            lemon_pe = sa.eva_tool(lemon_dir)
            print ("############ref" ,ba.sample, ref_accuracy, ref_len, "Mb", local_pe.accuracy,\
             local_pe.FDR, local_pe.F1_score ,lemon_pe.accuracy, lemon_pe.FDR, lemon_pe.F1_score)
            fi.add_local_sample(local_pe, depth)
            fi.add_lemon_sample(lemon_pe, depth)
    fi.plot()

def pure_snp():
    fi = Figure()
    ba = Parameters()
    true_dir = "/mnt/d/breakpoints/HGT/uhgg_snp/"
    lemon_dir = "/mnt/d/breakpoints/HGT/lemon_snp_pure/"
    local_dir = "/mnt/d/breakpoints/HGT/uhgg_snp_pure/"

    for snp_rate in ba.snp_level:
        ba.change_snp_rate(snp_rate)
        for index in range(10):
            ba.get_ID(index)    
            sa = Sample(ba.sample, true_dir)
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
        
            local_pe = sa.eva_tool(local_dir, 'localHGT')
            local_pe.add_ref(ref_accuracy, ref_len)
            lemon_pe = sa.eva_tool(lemon_dir, "lemon")
            fi.add_local_sample(local_pe, snp_rate)
            fi.add_lemon_sample(lemon_pe, snp_rate)
    fi.variation = "snp"
    fi.convert_df()
    fi.df.to_csv('/mnt/c/Users/swang66/Documents/For_methods/pure_snp_comparison.csv', sep=',')

def pure_length():
    fi = Figure()
    ba = Parameters()
    true_dir = "/mnt/d/breakpoints/HGT/uhgg_length/"
    lemon_dir = "/mnt/d/breakpoints/HGT/lemon_length_results/"
    local_dir = "/mnt/d/breakpoints/HGT/uhgg_length_results/"

    # for snp_rate in ba.snp_level:
    #     ba.change_snp_rate(snp_rate)
    for read_length in [75]:
        ba.reads_len = read_length
        for index in range(10):
            ba.get_ID(index)    
            sa = Sample(ba.sample, true_dir)
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
        
            local_pe = sa.eva_tool(local_dir, 'localHGT')
            local_pe.add_ref(ref_accuracy, ref_len)
            lemon_pe = sa.eva_tool(lemon_dir, "lemon")
            fi.add_local_sample(local_pe, read_length)
            fi.add_lemon_sample(lemon_pe, read_length)
            print ("#ref" ,ba.sample, ref_accuracy, ref_len, "Mb", local_pe.accuracy, local_pe.FDR)
    fi.variation = "length"
    fi.convert_df()
    fi.df.to_csv('/mnt/c/Users/swang66/Documents/For_methods/pure_length_comparison.csv', sep=',')

def pure_donor():
    fi = Figure()
    ba = Parameters()
    true_dir = "/mnt/d/breakpoints/HGT/uhgg_length/"
    lemon_dir = "/mnt/d/breakpoints/HGT/lemon_length_results/"
    local_dir = "/mnt/d/breakpoints/HGT/uhgg_length_results/"

    # for snp_rate in ba.snp_level:
    #     ba.change_snp_rate(snp_rate)
    for read_length in [75]:
        ba.reads_len = read_length
        for index in range(10):
            ba.get_ID(index)    
            sa = Sample(ba.sample, true_dir)
            ref_accuracy, ref_len = sa.eva_ref(local_dir)
        
            local_pe = sa.eva_tool(local_dir, 'localHGT')
            local_pe.add_ref(ref_accuracy, ref_len)
            lemon_pe = sa.eva_tool(lemon_dir, "lemon")
            fi.add_local_sample(local_pe, read_length)
            fi.add_lemon_sample(lemon_pe, read_length)
            print ("#ref" ,ba.sample, ref_accuracy, ref_len, "Mb", local_pe.accuracy, local_pe.FDR)
    fi.variation = "length"
    fi.convert_df()
    fi.df.to_csv('/mnt/c/Users/swang66/Documents/For_methods/pure_length_comparison.csv', sep=',')

if __name__ == "__main__":
    true_dir = "/mnt/d/breakpoints/HGT/uhgg_snp/"
    lemon_dir = "/mnt/d/breakpoints/HGT/lemon_snp/"
    local_dir = "/mnt/d/breakpoints/HGT/uhgg_snp_results_paper/"
    print ("evaluation")
    # cami()
    # snp()
    # pure_snp()
    # depth()
    # pure_length()



