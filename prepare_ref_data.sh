#!/bin/bash
mode=$1
WD=./${mode}
cd ${WD}

if [ "${mode}" == "Homo_sapiens" ]
then
    cp ../sars.fa .
    wget ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

    gzip -d *.gz
    mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
    cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > transcriptome.fa
    rm Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa
fi

if [ "${mode}" == "Chlorocebus_sabaeus" ]
then
    cp ../sars.fa .
    wget ftp://ftp.ensembl.org/pub/release-105/fasta/chlorocebus_sabaeus/dna/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-105/fasta/chlorocebus_sabaeus/cdna/Chlorocebus_sabaeus.ChlSab1.1.cdna.all.fa.gz +
    wget ftp://ftp.ensembl.org/pub/release-105/fasta/chlorocebus_sabaeus/ncrna/Chlorocebus_sabaeus.ChlSab1.1.ncrna.fa.gz
    
    gzip -d *.gz
    mv Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa genome.fa
    cat Chlorocebus_sabaeus.ChlSab1.1.cdna.all.fa Chlorocebus_sabaeus.ChlSab1.1.ncrna.fa > transcriptome.fa
    rm Chlorocebus_sabaeus.ChlSab1.1.cdna.all.fa Chlorocebus_sabaeus.ChlSab1.1.ncrna.fa
fi

if [ "${mode}" != "Homo_sapiens" ] && [ "${mode}" != "Chlorocebus_sabaeus" ]
then
    echo "Incorrect host name! Must be 'Homo_sapiens' or 'Chlorocebus_sabaeus'. Please try again. Exit."
    exit 1
fi
