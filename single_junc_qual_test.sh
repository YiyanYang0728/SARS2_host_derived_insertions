#!/bin/bash
srr=$1
i=$2
ref_junction=$3
junc_region=$4
map_region=$5
outdir=read_fq
mkdir -p ${outdir}
info=`python ../junction_quality_test.py "${junc_region}" "${ref_junction}" "${map_region}" ${outdir}/${i}.read.fq`
echo -e "${i}\t${info}"
