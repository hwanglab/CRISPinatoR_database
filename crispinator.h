
#ifndef CRISPINATOR_
#define CRISPINATOR_

#include "genomic_entry.h"
#include "bntseq.h"
#include "bwtaln.h"
#include "motif.h"

#define ESE_K	6
#define ESS_K	6
#define BP_PER_LINE	120
#define SGRNA_LEN	23
#define MAX_OFF_TARGET	1000

#define NORMAL	0
#define TRUNC_START	1
#define TRUNC_STOP	2
#define SKIPPED		4
#define FIND_ATG	8
#define PRINT_AA	16

// START of editing type flag
#define EDIT_START_TRSC		1
#define EDIT_END_TRSC		2	
#define EDIT_MIDDLE_TRSC	4	

#define EDIT_NONCODING		8
#define EDIT_5UTR		16
#define EDIT_START_CODON	32
#define EDIT_CODING		64
#define EDIT_STOP_CODON		128
#define EDIT_3UTR		256
#define EDIT_SPLICE_SITE	512
// END of editing type flag

#define W_ALPHA	1
#define W_BETHA	0.1

#define W_CODING	1
#define W_NONCODING	-1

#define TOP_N		100

#define W_POS		0.30
#define W_OFF		0.25
#define W_TRSC		0.25
#define W_ESE		0.2

#define PADDED_5	4
#define PADDED_3	3

#define FS_INS 0
#define FS_DEL 1

#define PAM_NGG	0
#define PAM_CCN 1

// struct representing a CRISPR/Cas9 sgRNA candidate
typedef struct crispr_candidate crispr_candidate;
struct crispr_candidate {
	unsigned int pos_in_chr;
	unsigned int pos_in_gene;
	unsigned int pam_pos;
	unsigned int cut_site;
	int dist5_ese, dist3_ese;
	int* n_in_gene;
	int* n_out_gene;
	int* pos_in_trsc;
	int n_ese;
	int n_aff_trsc;
	int n_aff_functional_trsc;
	double s_pos;
	double s_off;
	double s_trsc;
	double s_ese;
	double s_sum;
	bwa_seq_t *aln;
	gene *target_gene;
	transcript *target_lct;		// target transcript (longest coding transcript) containing target exon
	exon *target_exon;
	char* rna_seq;
	char* rna_seq_padded;
	int pam_direction;
};

int run_interactive_offtarget(int argc, char* argv[]);
int crispinator_write_library(int argc, char* argv[]);
int crispinator_annotate_base_editing(int argc, char* argv[]);

#endif

