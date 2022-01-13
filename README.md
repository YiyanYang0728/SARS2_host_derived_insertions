# SARS2_host_derived_insertions
This repository contains scripts and data for the direct RNA-seq analysis described in paper "Putative host-derived insertions in the genomes of circulating SARS-CoV-2 variants" by Yiyang Yang, Keith Dufault-Thompson, Rafaela Salgado Fontenele, Xiaofang Jiang.

## Prerequisites
- Python >=3.8
- minimap2
- paftools
- NanoFilt 
- datamash

## Usage
```
# git repo
git clone https://github.com/YiyanYang0728/SARS2_host_derived_insertions.git
cd SARS2_host_derived_insertions

# preparing reference data into host folder
./prepare_ref_data.sh Homo_sapiens

# for Chlorocebus_sabaeus, please run:
# ./prepare_ref_data.sh Chlorocebus_sabaeus

# Run on Linux
# Default host is Homo sapiens
./Pipeline.sh Homo_sapiens

# for Chlorocebus_sabaeus, please run:
# ./Pipeline.sh Chlorocebus_sabaeus

```

## Citation
Putative host-derived insertions in the genome of circulating SARS-CoV-2 variants
Yiyang Yang, Keith Dufault-Thompson, Rafaela Salgado Fontenele, Xiaofang Jiang
bioRxiv 2022.01.04.474799; doi: https://doi.org/10.1101/2022.01.04.474799
