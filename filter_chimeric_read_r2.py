# This script is designed to filter out reads with low-quality junctions
# if a read has multiple junctions, we only keep the junction with high quality (p-val>0.05).
import sys, csv
sys.path.append("/data/Irp-jiang/share/covid_data/direct_RNAseq/analysis_v2/Homo_sapiens")
sys.path.append("/data/Irp-jiang/share/covid_data/direct_RNAseq/analysis_v2/Chlorocebus_sabaeus")

PVAL_THERSHOLD = 0.2
infile = sys.argv[1]

# get ref info
junction_ref = {} # key: read+junction -> ref
with open("res_filtered.tsv") as f:
    for line in f:
        cols = line.strip().split("\t")
        read,length = cols[0].split(":")
        ref_list = cols[2].split(",")
        for index, jun in enumerate(cols[4].split(";")):
            junction_ref[read+":"+str(jun)] = ",".join(ref_list[index:index+2])

with open(infile) as csvfile:
    reader = csv.reader(csvfile,delimiter="\t")
    for row in reader:
        keep_hits = []
        read = row[0]
        sample = row[-1]
        # read 2nd & 3rd columns and combine accordingly
        hits = zip(row[1].split(";"),row[2].split(";"), row[3].split(";"), row[4].split(";"), row[5].split(";"), row[6].split(";"), row[7].split(";"))
        hits = list(hits)
        for hit in hits:
            hit = list(hit)
            jun = hit[1]
            if float(hit[-1]) >= PVAL_THERSHOLD:
                ref = junction_ref[read+":"+str(jun)]
                hit.append(ref)
                keep_hits.append(hit)

        if len(keep_hits) == 2:
            outrow = [i + ";" + j for i, j in zip(keep_hits[0], keep_hits[1])]
            print("\t".join([read]+outrow+[sample]))
        elif len(keep_hits) == 1:
            # print(keep_hits[0])
            print("\t".join([read]+list(keep_hits[0])+[sample]))
