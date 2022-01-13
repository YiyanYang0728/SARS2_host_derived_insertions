import sys, re
from collections import defaultdict

# Usage:
# python get_chimeric_pattern.py res_qual_spl_filtered.tsv > read_chimeric_pattern.tsv

infile = sys.argv[1]
pattern_dict = defaultdict(list)

def get_ref_strand(x):
    ref = x.split(":")[0]
    strandness = re.search(r"\((.*?)\)", x).group(1)
    if ref == "NC_045512.2":
        ref = "s"
    else:
        ref = "h"
    return ref+strandness

with open(infile) as f:
    for line in f:
        p_list = []
        cols = line.strip().split("\t")
        readId = cols[0]
        sample = cols[-1]
        mapping = cols[8].split(";")
        for index, map in enumerate(mapping):
            pattern = ""
            value = "\t".join([readId]+[i.split(";")[index] for i in cols[1:-1]]+[sample])
            for m in map.split(","):
                p = get_ref_strand(m)
                pattern += p
            pattern_dict[pattern].append(value)
            
for key in pattern_dict:
    if key in ["s+h+", "s+h-", "s-h+", "s-h-", "h+s+", "h+s-", "h-s+", "h-s-"]:
        # print(key, len(pattern_dict[key]))
        for ele in pattern_dict[key]:
            print(key+"\t"+ele)
