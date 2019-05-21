
make

./CRISPinatoR index data/hg19_chr/22.fa

samtools faidx data/hg19_chr/22.fa

./CRISPinatoR library -r data/hg19_chr/22.fa -d data/martquery_GRCh37_0819193157_65_tr_type_added.chr22.txt -e data/ese.txt -o hg19_sgRNA_list.txt
