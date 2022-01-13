#!/bin/bash

#####################################################
######### Pipeline to detect chimeric reads #########
#####################################################

mode=$1
WD=./${mode} # for Chlorocebus_sabaeus please change this variable accordingly
cd ${WD}

# 1. map fastq file to reference and convert the alignment to paf format
srr=`cat sra.list` # sra.list stores each srr sample name per line
for i in $srr; do ../fastq2paf.sh $i; done

# 2. extract host-virus chimeric reads from all samples and merged into one file
> res.paf
for i in $srr; do ../extract.sh $i >> res.paf; done

# 3. get junction start and end positions for each chimeric read and summarize them into a tab-delimited file
cat res.paf | cut -f 1-9 |sort -k1,1 -k3,3n -k4,4n | awk '{print $1":"$2"\t"$3"-"$4"\t"$6":"$7":"$8"-"$9"("$5")"}' | datamash -g 1 collapse 2,3 > res.tsv

# 4. filter out chimeric reads whose 1) junction length >= 15bp
# 2) junction occurs within the last 50 nucleotides of the first gene sequence or within the first 50 nucleotides of the second gene sequence
python ../filter_chimeric_read_r1.py res.tsv > res_filtered.tsv
# mapping sample information to each read
awk -F"\t" 'NR==FNR{a[$1]=$2;next}NR!=FNR{split($1,b,":");print a[b[1]]"\t"b[1]"\t"$(NF-1)"\t"$NF"\t"$2}' read_spl.mapping res_filtered.tsv > res_spl_read_junction.tsv

# 5. split fastq file into multiple single-read fastq files for following quality score checking
cut -f1,2 res_spl_read_junction.tsv | parallel -k --colsep '\t' "../read_fq.sh {}"

# 6. do permutation test to check if the quality score within 20 bp of either side of the junction region
# is higher than the 20th percentile quality score for each read
cat res_spl_read_junction.tsv | parallel -k --colsep '\t' "../single_junc_qual_test.sh {}" > res_qual.tsv
# mapping sample information to each junction
awk -F"\t" 'NR==FNR{a[$1]=$2;next}{print $0"\t"a[$1]}' read_spl.mapping res_qual.tsv > res_qual_spl.

# 7. filter out chimeric reads that don't meet the above demand
python ../filter_chimeric_read_r2.py res_qual_spl.tsv > res_qual_spl_filtered.tsv

# 8. annotate junctions with chimeric pattern (host-virus or virus-host)
python ../get_chimeric_pattern.py res_qual_spl_filtered.tsv | sed "s/\-/\+/g" | sed -e "s/s+h+/sh/g" -e "s/h+s+/hs/g"  > read_chimeric_pattern.tsv
