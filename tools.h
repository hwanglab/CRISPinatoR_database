/*
 *
 * tools.c     CRISPinatoR
 *
 *  Date: 2016. 6. 27.
 *  Author: Yunku Yeu
 */

#ifndef TOOLS_
#define TOOLS_

#include "genomic_entry.h"

#define ERROR			-1
#define FIRST_EXON		0
#define	LAST_EXON		1
#define SINGLE_EXON_GENE	2
#define NONCODING_EXON		4
#define NORMAL_EXON		8

double calculate_fs_offtarget_score(genome*, crispr_candidate*, bntseq_t*, FILE*);
int mark_and_list_kmers(const char*, char*, char*, kmers*, int, int);
int calculate_aa_len(genome*, transcript*, int);
char* convert_dna_to_aa_seq(const char*, unsigned int);
void print_seq(char*, unsigned int, unsigned int, int, FILE*, int, int, int);
unsigned int get_stop_codon_position_from_seq(transcript*, char*, unsigned int);
void list_known_trascripts(genome*, gene*);
void print_transcript(genome*, transcript*, int, FILE*, int, int, int);
int create_skipped_transcript(genome*, transcript*, transcript*, int, FILE*);
void free_skipped_transcript(transcript*);

/**BWA related functions**/
void seq_reverse(int, ubyte_t*, int);
void init_seq(bwa_seq_t*, int);
bwt_t* load_bwt(char*);
void do_alignment(genome*, bwt_t*, bntseq_t*, crispr_candidate*, bwa_seq_t*, int, gap_opt_t*);

#endif

