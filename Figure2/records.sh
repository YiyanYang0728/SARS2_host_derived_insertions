# get host gene annotation
host=Homo_sapiens # Chlorocebus_sabaeus
keyword=ENST # ENSCSAT
grep -oP "${keyword}[0-9]+.[0-9]+" <(cut -f3 ../results/${host}/read_chimeric_pattern.tsv) > 1.tmp
paste 1.tmp <(awk -F"\t" '{print $NF}' ../results/${host}/read_chimeric_pattern.tsv) > 2.tmp
awk 'NR==FNR{a[$1]=$0;next}$1 in a{print a[$1]"\t"$2}' ${host}.transcriptome.tsv 2.tmp > 3.tmp

echo -e "Transcript\tGene\tType1\tType2\tGeneSymbol\tSample\tCellline" > ${host}_chimeric_read_host_genes.tsv.1
awk -F"\t" 'NR==FNR{a[$1]=$2;next}$NF in a{print $0"\t"a[$NF]}' sample2celllin.tsv 3.tmp >> ${host}_chimeric_read_host_genes.tsv.1

rm *.tmp
