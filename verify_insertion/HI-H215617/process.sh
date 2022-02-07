#!/bin/bash
module load bowtie/2-2.4.4 bedtools/2.30.0 bedops/2.4.40 samtools/1.14

bamfile=H215617_trimmed_filtered.sorted.bam

samtools collate -u -O $bamfile | samtools fastq -@ 32 -1 paired1.fq -2 paired2.fq -0 /dev/null -s /dev/null -n

fq1=paired1.fq
fq2=paired2.fq

ref_list=("w_insert" "wo_insert")
for ref in ${ref_list[@]}; do
    bowtie2-build ${ref}.fasta ${ref}
    bowtie2 --xeq -p 32 -a -x ${ref} -1 paired1.fq -2 paired2.fq -S ${ref}.bowtie2.xeq.sam
    samtools view -@ 32 -F 4 -h ${ref}.bowtie2.xeq.sam | sam2bed > ${ref}.bowtie2.xeq.bed
    closestBed -b hit_${ref}.bed -a <(awk 'length($12)>=10{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$12}' ${ref}.bowtie2.xeq.bed) -d |awk '$12==0{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > ${ref}.extract.bowtie2.xeq.tsv
done

python ../insertion_match_reads.py w_insert.extract.bowtie2.xeq.tsv wo_insert.extract.bowtie2.xeq.tsv > assign_label.reads

# get reads exclusively aligned to the consensus genome with the insertion
awk '$NF=="w_insert"' assign_label.reads > w_insert.reads

# get reads exclusively aligned to the consensus genome without the insertion
awk '$NF=="wo_insert"' assign_label.reads > wo_insert.reads
