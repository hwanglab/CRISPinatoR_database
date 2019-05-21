# CRISPinatoR_database
This repository contains database generating source codes for **[CRISPinatoR](www.crispinator.com)**


This source codes includes following functions:
* Load Ensembl data including gene, transcript, and exon annotation
* Search sgRNA candidates from exonic regions defined by Ensembl data
* Annotate distance from the sgRNA candidates to closest ESE (Exon Splicing Enhancer)
* Off-target search for the sgRNA candidates
* Write the sgRNA candidates into a tab-delimetered text file.


To generate the database, below data is required. Example data can be downloaded **[here](https://drive.google.com/open?id=1_o31dnXGel2W-SdCw-XW3DtBaZdCJdVJ)**.
* Ensembl data downloaded from [here](http://useast.ensembl.org/biomart/martview) (Preprocessing is required). A whole data file for hg19 can be downloaded **[here](https://drive.google.com/open?id=13Yqca46UlRM7QEVLvvy7fLEEdRVBmdaB)**
* Reference FASTA file. It should be indexed by bwa and samtools faidx command.
* List of known ESE as a plain text file. For each line, the first tab-delimetered word will be regarded as an ESE. 


To test this source codes, please follow instructions below.
1. clone this repository.
```
git clone 
```


2. Download the example data above, then extract the zip file to a directory.
```
wget https://drive.google.com/open?id=1_o31dnXGel2W-SdCw-XW3DtBaZdCJdVJ
unzip example_data.zip
```


3.Compile source codes.
```
make
```


4.Then index the reference file.
```
./CRISPinatoR index example_data/22.fa
samtools faidx example_data/22.fa
```


5.Now you are ready to run the CRISPinatoR library construction. This command takes appoximately 1 hour to finish.
```
./CRISPinatoR library -r example_data/22.fa -d example_data/martquery_GRCh37_0819193157_65_tr_type_added.chr22.txt -e example_data/ese.txt -o hg19_sgRNA_list.txt
```
