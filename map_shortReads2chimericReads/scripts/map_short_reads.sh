#!/bin/bash
cores=${1:-32}
for srr in `cat data/sra.list`
do
    scripts/map_short_reads_per_sample.sh ${cores} $srr
done

awk -F"\t" 'BEGIN{OFS="\t";}NR==FNR{a[$1"\t"$2"\t"$3]=$4"\t"$5;next}NR!=FNR{split($4,b,",");\
if((b[1]<b[2]) && (($2"\t"b[1]"\t"b[2]) in a)){key=($2"\t"b[1]"\t"b[2]);print $11,$2,$4,$10,a[key]};\
if((b[1]>=b[2]) && (($2"\t"b[2]"\t"b[1]) in a)){key=($2"\t"b[2]"\t"b[1]);print $11,$2,$4,$10,a[key]}}' \
<(cat *.supported.read.tsv) data/read_chimeric_pattern.tsv > results/junction.supported.read.tmp

echo -e "BioSample\tDirectRNAseqRun\tShortReadSeqRun\tReadName\tJunctionRegion\tReferenceMapping\tSpanningJunctionReadCount\tSpanningJunctionReadPairCount" > results/final.supported.read.tsv
join -t $'\t' -a 1 -1 1 -2 1 -o 1.3,1.1,1.2,2.2,2.3,2.4,2.5,2.6 <(sort -k1,1 ../DRS2NGS2sample.tsv) <(sort -k1,1 junction.supported.read.tmp) | sort -t$'\t' -k1,1 -k7,7nr -k8,8nr   >> results/final.supported.read.tsv
rm results/junction.supported.read.tmp
