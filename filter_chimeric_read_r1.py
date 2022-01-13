#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import re
import sys

def read_line(hit):
    pos, match = hit
    pos_start, pos_end =pos.split("-")
    m = re.search("^(.*):(.*):(\d*)-(\d*)\(([+-])\)",match)
    hit_name = m.group(1)
    hit_length = int(m.group(2))
    hit_start = int(m.group(3))
    hit_end = int(m.group(4))
    hit_strand = m.group(5)
    return (int(pos_start),int(pos_end),hit_name,hit_start,hit_end,hit_length,hit_strand)

def get_junction_site(name, start, end, strand, order):
    if order == "former":
        if strand == "+":
            return name+":"+str(end)
        else:
            return name+":"+str(start)
    elif order == "latter":
        if strand == "+":
            return name+":"+str(start)
        else:
            return name+":"+str(end)

if __name__ == '__main__':
    infile = sys.argv[1]
    with open(infile) as csvfile:
        reader = csv.reader(csvfile,delimiter="\t")
        for row in reader:
            # read 1st column: readId:length
            read,length = row[0].split(":")
            # read 2nd & 3rd columns and combine accordingly
            hits = zip(row[1].split(","),row[2].split(","))
            hits = list(hits)

            junction_array = []
            ref_junction_array = []

            flag = False
            for index in range(1, len(hits)):
                prelin = hits[index-1] # get previous line            
                previous_line = read_line(prelin)
                pre_pos_start,pre_pos_end,pre_name,pre_start,pre_end,pre_length,pre_strand = previous_line
                
                line = hits[index] # get current line
                current_line = read_line(line)
                pos_start,pos_end,hit_name,hit_start,hit_end,hit_length,hit_strand = current_line
                
                # print(current_line, previous_line)
                
                if previous_line[2] != current_line[2] and (previous_line[2]=="NC_045512.2" or current_line[2]=="NC_045512.2" ): # junction 
                    if pre_pos_end + 15 > pos_start:  # junction gap should be less than 15
                        if ((pre_strand == "+" and pre_end <= pre_length -50 ) or (pre_strand == "-" and pre_start >=50 )) and ((hit_strand == "-" and hit_end <= hit_length -50 ) or (hit_strand == "+" and hit_start >=50 )): # not end 
                            flag = True
                            junction_array.append((pre_pos_end,pos_start))
                            jun_site_1 = get_junction_site(pre_name, pre_start, pre_end, pre_strand, "former")
                            jun_site_2 = get_junction_site(hit_name, hit_start, hit_end, hit_strand, "latter")
                            ref_junction_array.append((jun_site_1, jun_site_2))

            if flag:
                print("\t".join(row)+"\t"+";".join([a+","+b for a,b in ref_junction_array]) +"\t" +";".join([str(a)+","+str(b) for a,b in junction_array]))
            
        # if it is junction
            # previous hit is not end
            # current hit is not start
            # no huge space on junction
            # change flag to true
            # push the junction into read_array

        # if the coverage is enough 
        # push read_array to res_array
