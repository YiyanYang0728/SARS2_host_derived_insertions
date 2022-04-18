#!/bin/bash

srr=$1 # sample name
cores=${2:-32}

# 1. parse chimera to reference
cd ../data
DRS_sample=`grep ${srr} DRS2NGS | cut -f1`
chimeric_read_names=`grep ${DRS_sample} read_chimeric_pattern.tsv | cut -f2| sort | uniq`
read_dir=/gpfs/gsfs12/users/Irp-jiang/share/covid_data/direct_RNAseq/analysis_v2/Chlorocebus_sabaeus/read_fq
ref=${srr}_DRS_chimera
> ${ref}.fasta
for read_name in ${chimeric_read_names}
do
    awk -v n="${read_name}" 'NR==1{print ">"n}NR==2{print $0}' ${read_dir}/${read_name}.read.fq >> ${ref}.fasta
done

# 2. use bwa to map reads and convert them into paf format
module load bwa
bwa index ${ref}.fasta
# Notice: the follwoing two steps will take a long time, use HPC to run it if possible #
bwa mem -x ont2d -t ${cores} ${ref}.fasta ${srr}.f1.fq.gz ${srr}.r2.fq.gz > ${ref}.bwa2.sam
paftools.js sam2paf ${ref}.bwa2.sam> ${ref}.bwa2.paf

# 3. find short reads supportting DRS chimeric reads
flank=10
# get chimeric junction bed file
awk -F"\t" -v drs="${DRS_sample}" '$NF==drs{split($4,a,",");if(a[1]>a[2]){print $2"\t"a[2]"\t"a[1]"\toverlap"}else{print $2"\t"a[1]"\t"a[2]"\tgap"}}' read_chimeric_pattern.tsv | sortBed -i > ${ref}.bed
# get single read bed file
awk -F"\t" '{print $6"\t"$8"\t"$9"\t"$1}' ${srr}_DRS_chimera.bwa2.paf | sortBed -i > ${srr}_read.bed

# 3.1 single read spanning the junction sites
bedtools closest -a ${ref}.bed -b ${srr}_read.bed | \
awk -v flank="${flank}" '( $4=="gap" && ($6<=($2-flank) && $7>=($3+flank)) ) || ($4=="overlap" && ($6<=($2-flank) && $7>=($2+flank) || $6<=($3-flank) && $7>=($3+flank)) ) {print $0}' > ${srr}.1.supported.reads
datamash -s -g 1,2,3 count 8  < ${srr}.1.supported.reads > ${srr}.1.supported.count

# 3.2 the insert of a pair of reads spanning the junction sites
# get a pair of reads bed file
cut -f1,4 ${srr}_read.bed | cut -f1 -d"/" | sort | uniq -c | awk '$1==2{print $2"\t"$3}' > ${srr}.1.tmp
awk -F"\t" 'NR==FNR{a[$2"/1\t"$1];next}($1"\t"$6 in a){print $0}' ${srr}.1.tmp ${srr}_DRS_chimera.bwa2.paf | sed "s#/1##g" > ${srr}.1.1.tmp
awk -F"\t" 'NR==FNR{a[$2"/2\t"$1];next}($1"\t"$6 in a){print $0}' ${srr}.1.tmp ${srr}_DRS_chimera.bwa2.paf | sed "s#/2##g" > ${srr}.1.2.tmp
cut -f1,6 ${srr}.1.1.tmp | awk -F"\t" 'NR==FNR{a[$1"\t"$2];next}$1"\t"$6 in a{print $0}' - ${srr}.1.2.tmp | cut -f1,6 > ${srr}.2.tmp 
paste <(awk -F"\t" 'NR==FNR{a[$1"\t"$2];next}$1"\t"$6 in a{print $0}' ${srr}.2.tmp ${srr}.1.1.tmp | sort -k1,1 -k6,6) <(awk -F"\t" 'NR==FNR{a[$1"\t"$2];next}$1"\t"$6 in a{print $0}' ${srr}.2.tmp ${srr}.1.2.tmp | sort -k1,1 -k6,6) > ${srr}.3.tmp
# $4!=$7 different strands; get leftmost and rightmost region
awk -F"\t" 'BEGIN{OFS="\t"}{print $6,$8,$9,$5,$26,$27,$23,$1}' ${srr}.3.tmp | awk 'BEGIN{OFS="\t"}$4!=$7{if($6>=$3){right_lim=$6}else{right_lim=$3};if($5<=$2){left_lim=$5}else{left_lim=$2};print $1,left_lim,right_lim,$8}' | sortBed -i > ${srr}_pairread.bed
# a pair of reads should span the junction sites plus some flanks
bedtools closest -a ${ref}.bed -b ${srr}_pairread.bed | \
awk -v flank="${flank}" '($2!=-1) && (( $4=="gap" && ($6<=($2-flank) && $7>=($3+flank)) ) || ($4=="overlap" && ($6<=($2-flank) && $7>=($2+flank) || $6<=($3-flank) && $7>=($3+flank)) ))' > ${srr}.2.supported.reads
datamash -s -g 1,2,3 count 8  < ${srr}.2.supported.reads > ${srr}.2.supported.count

# 4. get final table
awk -F"\t" 'NR==FNR{a[$1"\t"$2"\t"$3]=$4;next}NR!=FNR{if($1"\t"$2"\t"$3 in a){key=$1"\t"$2"\t"$3;print key"\t"a[key]}else{key=$1"\t"$2"\t"$3;print key"\t0"}}' ${srr}.1.supported.count ${ref}.bed |\
awk -F"\t" 'NR==FNR{a[$1"\t"$2"\t"$3]=$4;next}NR!=FNR{if($1"\t"$2"\t"$3 in a){key=$1"\t"$2"\t"$3;print key"\t"$4"\t"a[key]}else{key=$1"\t"$2"\t"$3;print key"\t"$4"\t0"}}' ${srr}.2.supported.count - > ${srr}.supported.read.tsv
