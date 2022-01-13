#!/bin/bash
srr=$1
i=$2
outdir=read_fq
mkdir -p ${outdir}
module load seqtk
echo $i|seqtk subseq ${srr}.q10hc50.fq.gz - | sed '2s/U/T/g' > ${outdir}/${i}.read.fq

