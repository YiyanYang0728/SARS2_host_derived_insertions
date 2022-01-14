#!/bin/bash

####################################################################
########################## Data instruction ########################
# The raw direct RNA-sequencing reads for Homo sapiens samples 
# (SRR13089348, SRR13069479, SRR13089340, SRR13089336, SRR13089338, 
# and SRR13089334) were downloaded from Sequence Read Archive (SRA) 
# from the National Center for Biotechnology Information (NCBI) 
# database under accession number PRJNA675370. The raw direct 
# RNA-sequencing reads for Chlorocebus sabaeus samples (ERR4808851, 
# SRS6201528, SRS6201528, SRS6313628, SRS6313628, SRS6430809, 
# SRS6442399, SRS6442400, SRS6442401, SRS7237280, SRS7741765, 
# SRS7741763, SRS7741767, SRS7945957, SRS7945957, SRS7945956, 
# and SRS7945956) were downloaded from SRA from the NCBI database 
# under accession numbers PRJNA608224, PRJEB39337, PRJNA623285, 
# PRJNA658490, PRJNA675370, and PRJNA688696. The Chlorocebus sabaeus 
# sample VeroInf24h.fastq is downloaded from https://osf.io/8f6n9/ on 
# the Open Science Framework. The Chlorocebus sabaeus samples 
# (CRR197323, CRR197322, CRR197321, CRR197320, CRR127517, and 
# CRR127516) are downloaded from GSA:CRA002508 from the BIG Data 
# Center (https://bigd.big.ac.cn/) under accession number PRJCA002477. 
####################################################################

mode=$1
WD=./${mode}
cd ${WD}

if [ "${mode}" == "Homo_sapiens" ]
then
    cat sra.list | parallel fasterq-dump {}
fi

if [ "${mode}" == "Chlorocebus_sabaeus" ]
then
    grep -P "^ERR|^SRR" sra.list | parallel fasterq-dump {}
    # fastq files VeroInf24h.all.fastq, CRR197323.fastq, CRR197322.fastq, 
    # CRR197321.fastq, CRR197320.fastq, CRR127517.fastq, and CRR127516.fastq 
    # need to be downloaded from https://osf.io/8f6n9/ and PRJCA002477 at https://bigd.big.ac.cn/
fi

if [ "${mode}" != "Homo_sapiens" ] && [ "${mode}" != "Chlorocebus_sabaeus" ]
then
    echo "Incorrect host name! Must be 'Homo_sapiens' or 'Chlorocebus_sabaeus'. Please try again. Exit."
    exit 1
fi
