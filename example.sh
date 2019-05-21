
make

./CRISPinatoR index example_data/22.fa

samtools faidx example_data/22.fa

./CRISPinatoR library -r example_data/22.fa -d example_data/martquery_GRCh37_0819193157_65_tr_type_added.chr22.txt -e example_data/ese.txt -o hg19_sgRNA_list.txt

