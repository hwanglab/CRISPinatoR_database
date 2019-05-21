# CRISPinatoR_database
This repository contains database generating source codes for [CRISPinatoR](www.crispinator.com)

It includes following functions:
* Load Ensembl data including gene, transcript, and exon annotation
* Search sgRNA candidates from exonic regions defined by Ensembl data
* Annotate distance from the sgRNA candidates to closest ESE (Exon Splicing Enhancer)
* Off-target search for the sgRNA candidates

To generate the database, below data is required:
* Ensembl data from [here](http://useast.ensembl.org/biomart/martview). List of required field can be found at example data below.
* Reference FASTA file. It should be indexed by bwa and samtools faidx command.
* List of known ESE as a plain text file. For each line, the first tab-delimetered word will be regarded as an ESE. 

Example data can be downloaded [here](https://drive.google.com/open?id=1_o31dnXGel2W-SdCw-XW3DtBaZdCJdVJ)

To test this source codes, clone this repository and download the example data above, then extract the zip file to a directory.
