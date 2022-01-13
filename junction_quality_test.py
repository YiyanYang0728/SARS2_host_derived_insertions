############################# Purpose ##############################
# This script is designed to calculate the quality around juction 
# region in host-cov chimeric reads by doing random permutation test
# to get a quality distribution of 1000 randomly regions with equal 
# length in the whole reads. 
####################################################################

# Input: junction regions (start, end, length); 
# single read quality string (from fastq file)

# Process: calculate the mean quality of junction reads;
# sample 100 same-length regions and their avg quality score from the reads** with mapped region**;
# get p-value (times<=the junction quality score/1000 times)

# Output: junction info
# the quality distribution of 1000 random regions (plot)
# the pvalue

import sys
import numpy as np
from Bio import SeqIO
import random

def merge_region(regions):
    pre_region = regions[0]
    merged_region = [pre_region]
    for cur_region in regions:
        if pre_region!=cur_region:
            if cur_region[0]<pre_region[1]:
                pre_region = [pre_region[0], cur_region[1]]
                merged_region.pop()
                merged_region.append(pre_region)
            else:
                pre_region = cur_region
                merged_region.append(pre_region)
    return merged_region

def get_avg_qual(whole_qual_list, start, end):
    # start, end: 0-based
    return np.mean(whole_qual_list[start:end])

def permutate_region(whole_qual_list, sample_len, valid_region, times=10000):
    bg_avg_qual_list = []
    valid_region = [[i, j-sample_len+1] for i,j in valid_region if j-sample_len+1 > i]
    count = 0
    while 1:
        # randomly choose a start site from whole read
        if count >= times:
            break
        rand_start = random.choice(range(0, len(whole_qual_list)-sample_len+1))
        valid_sample = 0
        for r in valid_region:
            if rand_start <= r[1] and rand_start >= r[0]:
                valid_sample = 1
                break
        if valid_sample:
            rand_end = rand_start + sample_len
            avg_qual = get_avg_qual(whole_qual_list, rand_start, rand_end)
            bg_avg_qual_list.append(avg_qual)
            count += 1
    return bg_avg_qual_list

def get_perm_pval(score, bg_score_list):
    hit_len = len([i for i in bg_score_list if i <= score])
    return hit_len/len(bg_score_list)

# only handle one read per time
# 1. read junction ranges (might be neg or pos)
junction_region = sys.argv[1]
ref_junction_region = sys.argv[2]
# start_end_tuple = [(idx[0], idx[1]) for idx in eval(junction_region) if isinstance(idx, tuple)][0] # we assume only one elment in the tuple
start_end_tuples = [list(map(int, i.split(","))) for i in junction_region.split(";")]
ref_start_end_tuples = [i for i in ref_junction_region.split(";")]

# 2. get mapping region
mapping_region = sys.argv[3]
valid_region = [list(map(int, i.split("-"))) for i in mapping_region.split(",")]
valid_region = merge_region(valid_region)
# print(valid_region)

# 3. read avg quality score
fq = sys.argv[4]  #"/data/yangy34/projects/coronavirus_insertion/vis_read_quality/5385292c-e25b-4dda-9481-b7b798b4d7e3.read.fq"
for record in SeqIO.parse(fq, "fastq"):
    whole_qual_list = record.letter_annotations["phred_quality"]

# for storing results
p_value_list = []
avg_qual_list = []
region_start_list = []
region_end_list = []
region_len_list = []

for index, start_end_tuple in enumerate(start_end_tuples):
    ref_start_end_tuple = ref_start_end_tuples[index]
    start, end = start_end_tuple
    region_start = max(min(start, end)-20, 0) # add up 20 bp
    region_end = min(max(start, end)+20+1, len(whole_qual_list)) # add down 20 bp
    region_len = region_end - region_start
    avg_qual = get_avg_qual(whole_qual_list, region_start, region_end)
    region_start_list.append(str(region_start))
    region_end_list.append(str(region_end-1))
    region_len_list.append(str(region_len))
    avg_qual_list.append(str(avg_qual))
    
    # 3. sampling 1000 times to get quality score distribution
    sample_len = region_len
    bg_avg_qual_list = permutate_region(whole_qual_list, sample_len, valid_region, times=10000)
    # print(len(bg_avg_qual_list))

    # 4. calculate p-value
    p_value = get_perm_pval(avg_qual, bg_avg_qual_list)
    p_value_list.append(str(p_value))

print("\t".join([ref_start_end_tuple, junction_region, ";".join(region_start_list), ";".join(region_end_list), ";".join(region_len_list), ";".join(avg_qual_list), ";".join(p_value_list)]))
