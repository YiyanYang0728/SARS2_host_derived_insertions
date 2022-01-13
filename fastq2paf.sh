#!/bin/bash

#########################################
###### convert fastq to paf format ######
#########################################
srr=$1  # srr sample name

# 1. preprocess fastq file
NanoFilt -q 10 --headcrop 50 ${srr}.fastq | gzip > ${srr}.q10hc50.fq.gz

# 2. align reads to host genome, transcriptome and SARS-CoV-2 genome
cores=28 # CPUS
genome_file=genome.fa
transcriptome_file=transcriptome.fa
sars_file=sars.fa

minimap2 -ax splice -p 0.2 -k 14 -t ${cores} ${genome_file} ${srr}.q10hc50.fq.gz > ${srr}.genome.sam
minimap2 -ax map-ont --secondary=no -t ${cores} ${transcriptome_file} ${srr}.q10hc50.fq.gz > ${srr}.transcriptome.sam
minimap2 -ax map-ont -p 0.2 -t ${cores} ${sars_file} ${srr}.q10hc50.fq.gz > ${srr}.sars.sam

# 3. convert file from sam format to paf format
paftools.js sam2paf ${srr}.genome.sam > ${srr}.genome.paf
paftools.js sam2paf ${srr}.transcriptome.sam > ${srr}.transcriptome.paf
paftools.js sam2paf ${srr}.sars.sam > ${srr}.sars.paf
