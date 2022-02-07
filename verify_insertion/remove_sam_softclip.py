# This script is to remove softclipping in sam file 
# by providing ivar ref primer bed file.
# If a softclipping region falls into the bed regions,
# it will be removed, while other softclipping will be remained

import sys, os, re

sam_file = sys.argv[1]
bed_file = sys.argv[2] # first 3 cols: refname, start_pos, end_pos
out_file = sys.argv[3]

def parse_bed_file(bed_file):
    prime_regions = []
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            prime_regions.append((int(cols[1]), int(cols[2])))
    return prime_regions

def softclip_in_prime_region(softclip_region, prime_regions):
    s1, e1 = softclip_region[0], softclip_region[1]-1
    for region in prime_regions:
        s2, e2 = region
        if s1 <= s2 and e1 <= e2 and s2 <= e1:
            return s2, e1+1
        elif s2 <= s1 and e2 <= e1 and e2 >= s1:
            return s1, e2+1
        elif s2 <= s1 and e1 <= e2:
            return s1, e1+1
        elif s1 <= s2 and e2 <= e1:
            return s2, e2+1
    return None, None

def remove_softclip(cigar, seq, qual, ref_pos, prime_regions):
    comp = re.compile('[0-9]+[SHMNDIPX]')
    str_list = re.findall(comp, cigar)
    new_cigar = ""
    new_seq = seq
    new_qual = qual
    softclip_regions = []
    cur_pos = int(ref_pos)
    for indx, item in enumerate(str_list):
        span, letter = int(item[:-1]), item[-1]
        if letter in "MND":
            cur_pos += span
        if letter == "S":
            if indx == 0: # S is at the beginning
                softclip_start = cur_pos - span
                softclip_end = cur_pos
                softclip_region = softclip_start, softclip_end
                new_s, new_e = softclip_in_prime_region(softclip_region, prime_regions)
                # print(new_s, new_e)
                if new_s!=None and new_e!=None:
                    new_seq = new_seq[span:]
                    new_qual = new_qual[span:]
                else:
                    new_cigar += str(span) + letter
            if indx == len(str_list)-1: # S is at the end
                softclip_start = cur_pos
                softclip_end = cur_pos + span
                softclip_region = softclip_start, softclip_end
                new_s, new_e = softclip_in_prime_region(softclip_region, prime_regions)
                # print(new_s, new_e)
                if new_s!=None and new_e!=None:
                    new_seq = new_seq[:-span]
                    new_qual = new_qual[:-span]
                else:
                    new_cigar += str(span) + letter

        else:
            new_cigar += str(span) + letter
    return new_cigar, new_seq, new_qual

prime_regions = parse_bed_file(bed_file)

g = open(out_file, "w")
with open(sam_file) as f:
    for line in f:
        if line.startswith("@"):
            g.write(line)
        else:
            cols = line.strip().split("\t")
            ref_pos = cols[3]
            cigar = cols[5]
            seq = cols[9]
            qual = cols[10]
            # print(ref_pos, cigar, seq)
            new_cigar, new_seq, new_qual = remove_softclip(cigar, seq, qual, ref_pos, prime_regions)
            # print(new_cigar, new_seq)
            cols[5] = new_cigar
            cols[9] = new_seq
            cols[10] = new_qual
            # print("="*20)
            g.write("\t".join(cols)+"\n")
g.close()
