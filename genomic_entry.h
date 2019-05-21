/*
 * genomic_entry.h
 *
 *  Created on: 2016. 6. 24.
 *      Author: mummy
 */

#ifndef GENOMIC_ENTRY_H_
#define GENOMIC_ENTRY_H_

#include <stdio.h>
#include "utils.h"

typedef struct genome 		genome;
typedef struct chromosome 	chr;
typedef struct exon		exon;
typedef struct transcript	transcript;
typedef struct gene		gene;
typedef struct domain		domain;

#define GENE	1
#define TRSC	2
#define EXON	3
#define HASH_SIZE		2048		// Define number of hashes.
#define MAX_NAME_LEN		1024
#define MAX_ATTR_LEN		1024
#define FORWARD		1
#define REVERSE 	-1
#define ADD_3BP_BUFFER 	1
#define NO_BUFFER	0

// linked list struct with entry
typedef struct lentry lentry;
struct lentry {
	lentry* next;	// next list
	void* entry;	// point any entry
};

// Represent a whole genome
struct genome {
	chr** 		chrs;
	unsigned short	num_chrs;
	FILE*		fasta;
	char*		fasta_path;
	lentry**	gene_head;	// store all genes, be hashed for faster lookup
	lentry**	trsc_head;
	lentry**	exon_head;

	char*		prefix_gene;
	char*		prefix_trsc;
	char*		prefix_exon;
	char*		prefix_protein;
	int		plen_gene;
	int		plen_trsc;
	int		plen_exon;
	int		plen_protein;
};

struct chromosome {
	char* 		name;
	unsigned short 	idx;
	lentry*		genes;
	unsigned int	num_genes;
	unsigned int	length;			// store .fai information
	unsigned int	offset;
	unsigned short	linebases;
	unsigned short	linelength;
};

struct exon {
	unsigned int	ensembl_id;
	chr* 		chr;
	short		strand;
	unsigned int 	start;
	unsigned int 	end;
	short		start_phase;
	short		end_phase;
	gene*		gene;
	lentry*		first_trsc;
	//short		is_coding;
	short		is_functional;
};

struct domain {
	int interProID;
	char* description;
	int start;
	int end;
};

struct transcript {
	unsigned int 	ensembl_id;
	unsigned int	ensembl_protein_id;
	chr*			chr;
	short			strand;
	unsigned int 	tr_start;		// include 5' & 3' UTR
	unsigned int 	tr_end;
	unsigned int	cdna_start;
	unsigned int	cdna_end;
	unsigned int	trsc_length;
	lentry*		first_exon;
	lentry*		is_coding_exon;
	gene*		gene;
	char*		trsc_type;
	int 	num_exons;
	int	num_coding_exons;
	unsigned int	aa_length;
	short		functional_trsc;
	short		has_CDS;		// has normal ATG
	int		first_cdna_exon;
	lentry*		domains;
};

struct gene {
	unsigned int	ensembl_id;
	char* 		symbol;
	chr* 		chr;
	unsigned int 	start;			// gene's position in genome.
	unsigned int 	end;
	short			strand;
	lentry*			first_trsc;
	lentry*			first_exon;
	int	num_exons;
	int	num_transcripts;
	int	n_functional_trsc;
};

genome* load_genome(char*);
int load_ensembl_data(genome*, char*);
void free_genome(genome*);
chr* find_chr_name(genome*, char*);
chr* find_chr_idx(genome*, int);
char* get_subseq_by_name(genome*, char*, unsigned int, unsigned int);
char* get_subseq_by_idx(genome*, int, unsigned int, unsigned int);
gene* find_gene_by_name(genome*, char*);
int is_coding_exon(transcript*, int);
char* concat_exon_seqs(genome*, transcript*, int);
unsigned int convert_genomic_coord_to_trsc_coord(transcript*, unsigned int);
gene* find_gene_by_id(genome* gen, unsigned int);
transcript* find_trsc_by_id(genome*, unsigned int);
transcript* find_trsc_by_id_in_gene(gene*, unsigned int);
exon* find_exon_by_id(genome*, unsigned int);
int load_ensembl_domain(genome*, char*);

#endif /* GENOMIC_ENTRY_H_ */
