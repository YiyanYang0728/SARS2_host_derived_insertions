# This script is used to get the reads 
# that can be fully matched to an reference 
# with or without the insertion of interest. 

# all inputs are bed files with bowtie2 --xeq options

# Input: 2 bed files mapped to w_insert and wo_insert ref 
# within a specific region
# Output: Reads that belongs to w_insert or wo_insert ref

# Idea: a read can be considered to be exclusively mapped to a reference, 
# 1) if it is only mapped to one reference and no mismatch/indel in given region
# 2) if it is mapped 2 references, but mapped to one reference with the full match or with no mismatch/indel in given region
# Notice: Reads 1) mapped to both references with full match and 2) mapped to both references all with mismatch/indel in given regions
# will be discarded.

import sys, os, re

def variant_in_region(start_pos, cigar, insert_type="w_insert"):
    # for wo_insert, varaint cannot be >= 7121 and <= 7122
    if insert_type == "wo_insert":
        given_region = list(range(7117, 7122+1))
    # for w_insert, varaint cannot be >= 7121 and <= 7149
    elif insert_type == "w_insert":
        given_region = list(range(7121, 7149+1))

    comp = re.compile('[0-9]+[SHMNDIPX=]')
    str_list = re.findall(comp, cigar)

    curr_pos = start_pos
    for indx, item in enumerate(str_list):
        span, letter = int(item[:-1]), item[-1]
        if letter in "=DX":
            curr_pos += span
            letter_start = start_pos
            letter_end = curr_pos-1
            letter_region = list(range(letter_start, letter_end+1))
            start_pos = curr_pos
            # print(letter, span, letter_start, letter_end)
        # D: deletion; I : insertion; X: mismatch
        if insert_type=="wo_insert" and letter in "DXI" and len(set.intersection(set(letter_region), set(given_region)))!=0:
            return 1
        if insert_type=="w_insert" and letter in "DXI" and len(set.intersection(set(letter_region), set(given_region)))!=0:
            return 1
    return 0

w_insert_bed = sys.argv[1]
wo_insert_bed = sys.argv[2]

reads = []

valid_reads_w_insert = {}
with open(w_insert_bed) as f1:
    for line in f1:
        cols = line.strip().split("\t")
        start_pos = int(cols[1]) + 1
        cigar = cols[6]
        seq = cols[7]

        flag = variant_in_region(start_pos, cigar, insert_type="w_insert")
        if flag==0: # no varaint found in read region
            key = "|".join([cols[3], cols[4]])
            if key not in reads:
                reads.append(key)
            if key in valid_reads_w_insert:
                print("repeat read!", key)
            else:
                valid_reads_w_insert[key] = (cols[3], cols[4], cols[5], cols[6])

valid_reads_wo_insert = {}
with open(wo_insert_bed) as f2:
    for line in f2:
        cols = line.strip().split("\t")
        start_pos = int(cols[1]) + 1
        cigar = cols[6]
        seq = cols[7]

        flag = variant_in_region(start_pos, cigar, insert_type="wo_insert")
        if flag==0: # no varaint found in read region
            key = "|".join([cols[3], cols[4]])
            if key not in reads:
                reads.append(key)
            if key in valid_reads_wo_insert:
                print("repeat read!", key)
            else:
                valid_reads_wo_insert[key] = (cols[3], cols[4], cols[5], cols[6])

for key in reads:
    if key in valid_reads_wo_insert and key not in valid_reads_w_insert:
        print("\t".join(valid_reads_wo_insert[key])+"\two_insert")
    elif key in valid_reads_w_insert and key not in valid_reads_wo_insert:
        print("\t".join(valid_reads_w_insert[key])+"\tw_insert")
