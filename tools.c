/*
 *
 * tools.c     CRISPinatoR
 *
 *  Date: 2016. 6. 27.
 *  Author: Yunku Yeu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "crispinator.h"
#include "genomic_entry.h"
#include "kseq.h"	
#include "bwtaln.h"
#include "bntseq.h"
#include "bwase.h"
#include "motif.h"
#include "tools.h"
#include <math.h>

// bwa struct
extern unsigned char nst_nt4_table[256];

// Codon table, will be used for printing AA sequence
// See calculate_aa_len() function to see how this is used.
// Nucleotide representation: T(00) < C(01) < A(10) < G(11) in binary
// Accoding to the rule above, TTT will be represented 000000 in binary, then 0th element of this array.
// ATG is 100011 in binary --> 35th element of this array.
char codon_tab[64][4] = {"Phe", "Phe", "Leu", "Lue", "Ser", "Ser", "Ser", "Ser",
			"Tyr", "Tyr", "***", "***", "Cys", "Cys", "***", "Trp",
			"Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro",
			"His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg",
			"Ile", "Ile", "Ile", "Met", "Thr", "Thr", "Thr", "Thr",
			"Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg",
			"Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala",
			"Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly"
};

// Short version of the codon table
char AA_tab[64] = {'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S',
		   'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
		   'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P',
		   'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
		   'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T',
		   'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
		   'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A',
		   'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
};

// Requried values for MIT off-target score calculation
double fs_m[20] = {0.000, 0.000, 0.014, 0.000, 0.000,
		   0.395, 0.317, 0.000, 0.389, 0.079, 
		   0.445, 0.508, 0.613, 0.851, 0.732, 
		   0.828, 0.615, 0.804, 0.685, 0.583};


// Calculate MIT off-target score. Unfortunately, the shut down the web page. 
// Parameters
// *********************************************************************************
//   genome* gen: 	Genome struct
//   crispr_candidate* cand: A sgRNA candidate struct
//   bntseq_t* bns: 	BWA related data structure that is required to retrieve alignment position
//   FILE* stream: 	Log stream. NULL can be given (No output)
// *********************************************************************************
// Return value: Off-target score for a sgRNA candidate
// *********************************************************************************
double calculate_fs_offtarget_score(genome* gen, crispr_candidate* cand, bntseq_t *bns, FILE* stream) {
	int i, j, k, n_mm, pos, seqid, *diff = NULL, aln_pos;//, noticed = 0, num_to_be_printed = 10;
	double overall_score, single_score, term1, term2, term3;
	char *q_seq, *r_seq;
	bwa_seq_t *p;

	r_seq = NULL;
	overall_score = 0.0;
	p = cand->aln;

	//if(stream != NULL)
	//	fprintf(stream, "Cand name = %s (Chromosome %s, pos: %u)\n", cand->aln->name, (cand->target_gene != NULL) ? cand->target_gene->chr->name : "NA", cand->pos_in_chr);

	q_seq = (char*)malloc(sizeof(char) * (p->full_len + 1));
	for(i = 0; i < p->full_len; i++) {
		q_seq[i] = "ACGT"[p->seq[i]%4];
	}
	q_seq[p->full_len] = '\0';

	// Convert sgRNA sequence into NGG form
	if(q_seq[0] == 'C') {
		q_seq = create_rc_free_exist(q_seq);
	}

	// BWA store main alignment at the front of their mapping list
	bns_cnt_ambi(bns, p->pos, p->len, &seqid);
	aln_pos = (int)(p->pos - bns->anns[seqid].offset + 1);
	if(stream != NULL) {
		//fprintf(stream, "Main alignment: %s:%d\t", bns->anns[seqid].name, aln_pos);
		//fprintf(stream, "Another aln: %d alignments (first %d alignments are printed.)\n", p->n_multi, num_to_be_printed);
		fprintf(stream, "#OT= %d\n", p->n_multi);
	}

	//printf("Main alignment: %s:%d\n", bns->anns[seqid].name, aln_pos);
	//printf("Another aln: %d alignments (first 100 alignments are printed.)\n", p->n_multi);

	//if(p->n_multi == 0)
	//	fprintf(stderr, "ERROR: No alignment\n");

	// Traverse all alignments
	i = 0;
	while(i < p->n_multi) {
		bwt_multi1_t *q = p->multi + i;	
		bns_cnt_ambi(bns, q->pos, p->len, &seqid);
		aln_pos = (int)(q->pos - bns->anns[seqid].offset + 1);
		//printf("Aln %d: %d-\n", i, aln_pos);
		
		// Get reference sequence for each mapping
		r_seq = get_subseq_by_idx(gen, find_chr_name(gen, bns->anns[seqid].name)->idx, (unsigned int)aln_pos, (unsigned int)(aln_pos + SGRNA_LEN - 1));

		if(stream != NULL) {
			fprintf(stream, "OT %d\n", i);
		}
		
		if(r_seq == NULL) {
			//if(stream != NULL)
			//	fprintf(stream, "Candidate %s matched to %s:%d\n", cand->aln->name, bns->anns[seqid].name, aln_pos);
			for(j = i + 1; j < p->n_multi; ++j)
				p->multi[j - 1] = p->multi[j]; 
			p->n_multi--;
			if(stream != NULL)
				fprintf(stream, "Wrong position - filtered\n");

			continue;
		}
		// Because bwa converts N in the reference genome into random nucleotide when it is indexed,
		// our 23-bp sgRNA sequence might match to the NNN...NNN area by random chance (= 1/(4^23) ). 
		// Moreover, we accept 3 mismatches thus the chance is much greater than the random.
		else if(r_seq[0] == 'N' || r_seq[0] == 'n') {
			//if((stream != NULL) && (!noticed)) {
			//	fprintf(stream, "Candidate %s matched to NNN..NNN area\n", cand->aln->name);
			//	fprintf(stdout, "Candidate %s matched to NNN..NNN area\n", cand->aln->name);
			//	noticed = 1;
			//}

			for(j = i + 1; j < p->n_multi; ++j)
				p->multi[j - 1] = p->multi[j]; 
			p->n_multi--;

			free_cnt(r_seq);
			if(stream != NULL)
				fprintf(stream, "N sequences - filtered\n");
			continue;
		}

		make_str_hc(r_seq);

		if(q->strand)
			r_seq = create_rc_free_exist(r_seq);
		// r_seq and q_seq are either FORWARD strand now
		
		//if(stream != NULL && i < num_to_be_printed) {
		//	fprintf(stream, "  Alignment %d\t%s:%d (%d mm)\n", i, bns->anns[seqid].name, aln_pos, q->mm);
		//	fprintf(stream, "    Reference  : %s\n", r_seq);
		//	fprintf(stream, "    Query(cand): %s\n", q_seq);
		//}

		if(q->mm > 0)
			diff = (int*)malloc_cnt(sizeof(int) * q->mm);
		if(stream != NULL)
			fprintf(stream, "#MM= %d\nMM positions= ", q->mm);
		pos = n_mm = 0;

		term1 = 1.0;
		for(j = 0; j < strlen(q_seq) - 3; j++) {
			if(r_seq[j] != q_seq[j]) {
				if(n_mm > q->mm) {
					//fprintf(stderr, "Query and Ref are diff:\n  Query: %s\n  Refer: %s\n", r_seq, q_seq);
					break;
				}
				// diff record mismatch position
				if(q->mm > 0)
					diff[n_mm] = pos;

				// Accodring to the mismatch position, different panalty was applied
				// Term 1 summarizes the penalties.
				term1 *= (1.0 - fs_m[pos]);
				if(stream != NULL)
					fprintf(stream, "%d ", pos);
				n_mm++;
				//if(stream != NULL && i < num_to_be_printed) {
				//	fprintf(stream, "      Diff at %d, e[M] = %f\n", pos, fs_m[pos]);
				//}
			}
			pos++;
		}
			
		if(stream != NULL)
			fprintf(stream, "\n");

		if(n_mm > 1) {
			// Term 2 represents mean pairwise distance between mismatches
			term2 = 0.0;
			for(j = 0; j < n_mm; j++) {
				for(k = j + 1; k < n_mm; k++) {
					term2 += abs(diff[j] - diff[k]);
				}
			}
			term2 = (term2 / n_mm) - 1.0;   // mean pairwise distance
			if(stream != NULL)
				fprintf(stream, "MPD= %f\n", term2);
			term2 = 1.0 / (((19.0 - term2) / 19.0 * 4.0) + 1);
		}
		else 
			term2 = 1.0;

		// Term 3 represents number of mismatches
		if(n_mm > 1)
			term3 = 1.0 / (double)(n_mm * n_mm);
		else
			term3 = 1.0;

		// Three terms are multiplied => score for A OFF-TARGET. 
		single_score = term1 * term2 * term3 * 100.0;
		//if(stream != NULL && i < num_to_be_printed) {
		if(stream != NULL) {
			fprintf(stream, "SCORE= %f * %f * %f = %f\n", term1, term2, term3, single_score);
		}

		// Off-target score was summed.
		overall_score += single_score;

		free_cnt(r_seq);
		if(q->mm > 0)
			free_cnt(diff);

		i++;
	}

	// Convert sum of single off target scores == MIT score from Feng Zhang lab.
	overall_score = 100.0 / (100.0 + overall_score) * 100.0;
	if(stream != NULL)
		fprintf(stream, "OVERALL_SCORE= %f\n", overall_score);
	
	free_cnt(q_seq);

	return overall_score;
}


// return type of exon: FIRST_EXON or LAST_EXON or SINGLE_EXON_GENE or NONCODING_EXON or NORMAL_EXON
int get_exon_type(transcript* trsc, exon* target_exon) {
	int start, end, exon_size;
	exon* cur_exon;
	lentry *cursor;

	start = trsc->cdna_start;
	end   = trsc->cdna_end;
	FOR_LENTRY(cursor, trsc->first_exon) {
		cur_exon = (exon*)cursor->entry;
		exon_size = cur_exon->end - cur_exon->start + 1;

		if(cur_exon == target_exon) {
			if(exon_size < start) {			// EXON contains 5'UTR only
				return NONCODING_EXON;
			}
			else if((start <= 0) && (end <= 0)) {	// EXON contains 3'UTR only
				return NONCODING_EXON;
			}
			// here, either / both of start and end are positive
			else if((start > 0) && (exon_size >= start)) {	// coding sequence starts in this exon
				if(exon_size >= end) 	// if coding sequence ends in the same exon --> single exon gene
					return SINGLE_EXON_GENE;
				else
					return FIRST_EXON;
			}
			else if(exon_size >= end) {
				return LAST_EXON;
			}
			return NORMAL_EXON;
		}
		start -= exon_size;
		end   -= exon_size;
	}

	return ERROR;
}

// return exon type, in string format
const char* get_exon_type_str (int type) {
	if(type == NORMAL_EXON)
		return "Normal exon";
	else if(type == FIRST_EXON)
		return "First coding exon";
	else if (type == LAST_EXON) 
		return "Last coding exon";
	else if (type == SINGLE_EXON_GENE)
		return "Single exon gene";
	else if (type == NONCODING_EXON)
		return "Noncoding exon";
	return "Uncategorized exon";
}

// This function find k-mers in the input 'seq', and writes the position of k-mers to the list_seq
// It assumes that 
//   1) seq, mark_seq, list_seq, and mers should not be NULL
//   2) seq is already converted into upper case
//   3) mark_seq is a copy of seq
//
// Parameters ***************************************************************
//   const char* seq: 		sgRNA sequence to be annotated with ESE/ESS. Should be upper case string.
//   char* mark_seq:  		Copy of the sgRNA sequence. If ESE/ESS found, it will be changed lower case
//   char* list_seq:		string contains locations of found ESE/ESS
//   kmers* mers:		k-mer struct representing ESE or ESS
//   int padded_len_start, int padded_len_end:	
//   				Length of padded sequence on both ends   				   
// **************************************************************************
// Return: Number of k-mers found in the input 'seq'
// **************************************************************************
int mark_and_list_kmers(const char* seq, char* mark_seq, char *list_seq, kmers *mers, int padded_len_start, int padded_len_end) {
	int i, j, k, len = strlen(seq), cnt = 0, start, end;
	char buf[BUF_LEN];

	// For each character of 'seq'
	for(i = 0; i < len; i++) {
		// Compare all K-mers in 'mers'
		for(j = 0; j < mers->length; j++) {
			// size of each mer is stored in k_size member variable and it can be vary.
			// So determine searching range for each k-mer
			start = (padded_len_start >= mers->k_size[j]) ? padded_len_start - mers->k_size[j] : -1;
			end = (mers->k_size[j] <= padded_len_end) ? len - padded_len_end : len - mers->k_size[j] + 1;
			if((start < i) && (i < end)) {
				// Find k-mer in the 'seq'
				if(strncmp(seq + i, mers->mers[j], mers->k_size[j]) == 0) {
					cnt++;
					// append location of k-mer to the list_seq string
					if(list_seq != NULL) {
						if(strlen(list_seq) > 0)
							strcat(list_seq, " / ");
					}
					
					sprintf(buf, "%d", j);
					if(list_seq != NULL)
						strcat(list_seq, buf);
					
					// mark k-mer as lower case
					for(k = 0; k < mers->k_size[j]; k++) {	if(mark_seq[i + k] >= 'A' && mark_seq[i + k] <= 'Z') 
						if(mark_seq[i + k] >= 'A' && mark_seq[i + k] <= 'Z') 
							mark_seq[i + k] -= 'A' - 'a';
					}

					break;
				}
			}
		}
	}

	return cnt;
}

/* BWA related functions */
void seq_reverse(int len, ubyte_t *seq, int is_comp)
{
	int i;
	if (is_comp) {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			if (tmp < 4) tmp = 3 - tmp;
			seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			seq[len-1-i] = seq[i]; seq[i] = tmp;
		}
	}
}

void init_seq(bwa_seq_t *p, int length) {
	p->clip_len = p->len = p->full_len = length;
	p->seq  = (ubyte_t*)malloc_cnt(p->full_len * 1);
	p->qual = (ubyte_t*)malloc_cnt(p->full_len * 1);
	p->rseq = (ubyte_t*)malloc_cnt(p->full_len * 1);
	p->md = NULL;
	p->cigar = NULL;
	p->multi = NULL;
	return;
}

bwt_t* load_bwt(char *path) {
	bwt_t* bwt;
	char *str = (char*)malloc_cnt((strlen(path) + 10) * sizeof(char));
	strncpy(str, path, strlen(path) + 1);
	strcat(str, ".bwt");  
	bwt = bwt_restore_bwt(str);

	strncpy(str, path, strlen(path) + 1);
	strcat(str, ".sa");  
	bwt_restore_sa(str, bwt);
	free_cnt(str);

	return bwt;
}

// Off-target search function.
void do_alignment(genome *gen, bwt_t *bwt, bntseq_t *bns, crispr_candidate *candidates, bwa_seq_t *seqs, int n_seq, gap_opt_t *opt) {
	int i, j, best, cnt, n_multi, seqid;
	bwa_seq_t *p;

	// do alignment
	bwa_cal_sa_reg_gap(0, bwt, n_seq, seqs, opt);

	// this loop includes three funtions in bwa. bwa_aln2seq_core(), bwa_cal_pac_pos(), bwa_print_sam1()
	for(i = 0; i < n_seq; i++) {
		int z = 0, rest, strand;
		//bwtint_t l;
		p = &seqs[i];
		if(!p->n_aln) { 	// n_aln is number of alignment. It cannot be null, because the inputs (sgRNAs sequences) were extracted from reference genome. 
			p->n_multi = 0;
			p->type = BWA_TYPE_NO_MATCH;
			continue;
		}
		p->ref_shift = 0;	// indel is out of scope
		//for(l = p->aln[0].k; l <= p->aln[0].l; ++l) {
		//}
		//p->sa = p->aln[0].k + (bwtint_t)((p->aln[0].l - p->aln[0].k + 1) * drand48());  // Basically, bwa chooses a main alignment position randomly among best matches
		best = p->aln[0].score;	// For us, it is not important
		for (j = cnt = 0; j < p->n_aln; ++j) {
			const bwt_aln1_t *q = p->aln + j;
			if (q->score > best) break; 
			cnt += q->l - q->k + 1;	
		}
		
		p->c1 = cnt;
		for (; j < p->n_aln; ++j) cnt += p->aln[j].l - p->aln[j].k + 1;
		p->c2 = cnt - p->c1;
		p->type = p->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;

		rest = MAX_OFF_TARGET;
		if(p->multi != NULL) {
			free(p->multi);
		}
		p->multi = malloc_cnt(rest * sizeof(bwt_multi1_t));
		for (j = 0; (j < p->n_aln) && (z < rest); ++j) { // j starts from 0, it means the best alignment is also stored in p->multi
			const bwt_aln1_t *q = p->aln + j;
			bwtint_t l;
			for (l = q->k; (l <= q->l) && (z < rest); ++l) {
				p->multi[z].pos = l;
				//p->multi[z].gap = q->n_gapo + q->n_gape;
				//p->multi[z].ref_shift = (int)q->n_del - (int)q->n_ins;
				p->multi[z].gap = 0;
				p->multi[z].ref_shift = 0;
				p->multi[z++].mm = q->n_mm;
				//printf("%d / %d\n", z, rest);
			}
		}
		p->n_multi = z;

		// convert SA interval to position
		//fprintf(result, "SA calculation\n");
		//p->pos = bwa_sa2pos(bns, bwt, p->sa, p->len + p->ref_shift, &strand);	
		//p->strand = strand;
		//if (p->pos == (bwtint_t)-1) 
		//	p->type = BWA_TYPE_NO_MATCH;
		p->pos = 0;
		p->sa = 0;
		for (j = n_multi = 0; j < p->n_multi; ++j) {
			bwt_multi1_t *q = p->multi + j;
			q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand);	
			q->strand = strand;
			bns_cnt_ambi(bns, q->pos, p->len, &seqid);
			if (q->pos != (bwtint_t)-1) {
				if( ( (int)(q->pos - bns->anns[seqid].offset + 1) == candidates[i].pos_in_chr) &&
						(find_chr_name(gen, bns->anns[seqid].name) == candidates[i].target_gene->chr) ){
					p->pos = q->pos;
					p->strand = q->strand;
				}
				else {
					p->multi[n_multi++] = *q;
				}
			}
		}
		p->n_multi = n_multi;
	}
}

// end of bwa functions


// This function calculate length of Amino Acid for the input transcript
// Parameters ***************************************************************
//   genome* gen: 		genome struct
//   transcript* trsc: 		Input transcript
//   int flag:			If FIND_ATG is set with this flag, this function trys to find the first ATG in its sequence.
//    				  Otherwise, trsc->cdna_start and phase info will be used to find the location of the ATG.
// **************************************************************************
// Return: Length of amino acid sequence
// **************************************************************************
int calculate_aa_len(genome* gen, transcript* trsc, int flag) {
	// 3bp buffer is added to transcript due to phase. 
	// Some partial transcripts may have phase information in the first exon.
	unsigned int cdna_start = trsc->cdna_start + 3, i, j, len;//cdna_end = trsc->cdna_end;
	int code, term = 0;
	char *gen_seq;
	int aa_cnt = 0;

	// Generate DNA sequence by concatenating exons
	gen_seq = concat_exon_seqs(gen, trsc, ADD_3BP_BUFFER);
	len = strlen(gen_seq);
	//printf("Genomic sequence length = %d\n", len - 3);

	if(((exon*)(trsc->first_exon->entry))->start_phase > 0 )
		cdna_start -= ((exon*)(trsc->first_exon->entry))->start_phase;

	// make CDS coordinates to 0-base
	cdna_start--;
	//cdna_end--;

	if(flag & FIND_ATG) {
		// find a start codon. phase information is ignored
		while(cdna_start < strlen(gen_seq) - 3) {
			if(((gen_seq[cdna_start]   == 'a') || (gen_seq[cdna_start]   == 'A')) &&
				((gen_seq[cdna_start+1] == 't') || (gen_seq[cdna_start+1] == 'T')) &&
				((gen_seq[cdna_start+2] == 'g') || (gen_seq[cdna_start+2] == 'G')))
				break;

			cdna_start++;	
		}
	}
	
	// fill amino acid sequences
	i = 0;
	while( i < len ) {
		if(term)
			break;

		else if((!term) && (i >= cdna_start) && (i+2 < len)) {
			code = 0;
			//printf ("AA %d: DNA = ", aa_cnt + 1);
			for(j = 0; j < 3; j++) {
				code = code << 2;
				//printf("%c", gen_seq[i+j]);
				if((gen_seq[i+j] == 'T') || (gen_seq[i+j] == 't')) code = code | 0;
				else if((gen_seq[i+j] == 'C') || (gen_seq[i+j] == 'c')) code = code | 1;
				else if((gen_seq[i+j] == 'A') || (gen_seq[i+j] == 'a')) code = code | 2;
				else if((gen_seq[i+j] == 'G') || (gen_seq[i+j] == 'g')) code = code | 3;
				else {
					code = 999;   // when N is in the codon 
					break;
				}
			}
			//printf(" => ");

			if(code < 64) {
				//printf("%s\n", codon_tab[code]);
				if(strcmp(codon_tab[code], "***") == 0)
					term = 1;
				else
					aa_cnt++;	// stop codons aren't included AA count
			}	

			i = i + 3;
		}
		else {
			i++;
		}
	}

	free_cnt(gen_seq);

	//getchar();

	return aa_cnt;
}

// Variation of calculate_aa_len.
// This returns AA seq of given DNA sequence. Returned string should be freed by caller
char* convert_dna_to_aa_seq(const char* gen_seq, unsigned int cdna_start) {
	int i, j, cur, code, max_len = strlen(gen_seq) / 3, len = strlen(gen_seq);
	char* aa_seq = malloc_cnt(sizeof(chr) * (max_len + 1));

	cdna_start--;

	aa_seq[max_len] = '\0';

	for(i = 0; i < max_len; i++) {
		aa_seq[i] = ' ';
	}

	// fill amino acid sequences
	cur = 0;
	for(i = cdna_start; i < len; ) {
		if(i + 2 < len) {
			code = 0;
			for(j = 0; j < 3; j++) {
				code = code << 2;
				if((gen_seq[i+j] == 'T') || (gen_seq[i+j] == 't')) code = code | 0;
				else if((gen_seq[i+j] == 'C') || (gen_seq[i+j] == 'c')) code = code | 1;
				else if((gen_seq[i+j] == 'A') || (gen_seq[i+j] == 'a')) code = code | 2;
				else if((gen_seq[i+j] == 'G') || (gen_seq[i+j] == 'g')) code = code | 3;
				else {
					code = 999;   // when N is in the codon 
					break;
				}
			}
			//printf("seq = %c%c%c, code = %d, ", gen_seq[i], gen_seq[i+1], gen_seq[i+2], code);
			if(code < 64) {
				aa_seq[cur++] = AA_tab[code];
				if(AA_tab[code] == '*')
					break;
			}
			else {
				aa_seq[cur++]   = ' ';
			}
			//printf("aa seq = %c\n", aa_seq[cur]);

			i = i + 3;
		}
		else {
			break;
		}
	}
	aa_seq[cur] = '\0';

	return aa_seq;
}

// Function for printing DNA & AA sequence to the output stream
// Parameters
// *********************************************************************************
//   char* gen_seq:		DNA sequence to be printed
//   unsigned int cdna_start, cdna_end:	 start and end of coding region.
//   int flag:			If PRINT_AA is set, DNA and AA sequence will be printed.
//   				Otherwise, only DNA sequence will be printed.
// 				If FIND_ATG is set with this flag, this function trys to find the first ATG in its gen_seq.
//    				Otherwise, cdna_start and phase info will be used to find the location of the ATG.
//   FILE* stream:		Output stream
//   int hl_start, hl_end: 	start and end for highlighted (by '@') region
//   int tap: 			Size of empty space before each line
// *********************************************************************************
// Return value: void
// *********************************************************************************

void print_seq(char* gen_seq, unsigned int cdna_start, unsigned int cdna_end, int flag, FILE* stream, int hl_start, int hl_end, int tap) {
	unsigned int cur, start, end, i, j, len = strlen(gen_seq);
	int code, term = 0;
	char *upper, *aa_seq = NULL;
	
	//printf("cdna: %u-%u\n", cdna_start, cdna_end);
	//printf("FLAG = %d\n", flag);

	// make coordinates to 0-base
	// when print AA, these coordinates represents ranges of cdna.
	// When print gene sequence only, these coordinates mean start / end of sequence.
	cdna_start--;
	cdna_end--;


	if(flag & PRINT_AA) {
		cur = 1;
		aa_seq = (char*)malloc_cnt(sizeof(char) * (len + 1));
		aa_seq[len] = '\0';

		for(i = 0; i < len; i++) {
			aa_seq[i] = ' ';
		}

		if(flag & FIND_ATG) {
			// find a start codon. phase information is ignored
			while(cdna_start < strlen(gen_seq) - 3) {
				if(((gen_seq[cdna_start]   == 'a') || (gen_seq[cdna_start]   == 'A')) &&
						((gen_seq[cdna_start+1] == 't') || (gen_seq[cdna_start+1] == 'T')) &&
						((gen_seq[cdna_start+2] == 'g') || (gen_seq[cdna_start+2] == 'G')))
					break;

				cdna_start++;	
			}
		}

		// fill amino acid sequences
		for(i = 0; i < len; ) {
			if((!term) && (i >= cdna_start) && (i+2 < len)) {
				code = 0;
				for(j = 0; j < 3; j++) {
					code = code << 2;
					if((gen_seq[i+j] == 'T') || (gen_seq[i+j] == 't')) code = code | 0;
					else if((gen_seq[i+j] == 'C') || (gen_seq[i+j] == 'c')) code = code | 1;
					else if((gen_seq[i+j] == 'A') || (gen_seq[i+j] == 'a')) code = code | 2;
					else if((gen_seq[i+j] == 'G') || (gen_seq[i+j] == 'g')) code = code | 3;
					else {
						code = 999;   // when N is in the codon 
						break;
					}
				}
				//printf("seq = %c%c%c, code = %d, ", gen_seq[i], gen_seq[i+1], gen_seq[i+2], code);
				if(code < 64) {
					aa_seq[i]   = codon_tab[code][0];
					aa_seq[i+1] = codon_tab[code][1];
					aa_seq[i+2] = codon_tab[code][2];
					if(strcmp(codon_tab[code], "***") == 0)
						term = 1;
				}
				else {
					aa_seq[i]   = ' ';
					aa_seq[i+1] = ' ';
					aa_seq[i+2] = ' ';
				}
				//printf("aa seq = %c%c%c\n", aa_seq[i], aa_seq[i+1], aa_seq[i+2]);

				i = i + 3;
			}
			else {
				aa_seq[i]   = ' ';
				i++;
			}
		}
	}
	else {
		cur = cdna_start;
		len = cdna_end;
	}

	// upper line == header
	upper  = (char*)malloc_cnt(sizeof(char) * (len + 1));
	upper[len]  = '\0';

	for(i = cur; i <= len; i++) {
		upper[i-1]  = ((i % 10) == 0) ? '#' : '=';
	}

	if((hl_start > 0) && (hl_end > 0)) {
		for(i = hl_start; i <= hl_end; i++)
			upper[i-1] = '@';
	}


	// print out three strings
	while(cur <= len) {
		start = cur;	
		end = cur + BP_PER_LINE;

		if((start < cdna_start + 1) && (cdna_start <= end + 1)) {
			end += (cdna_start + 1 - start) % 3;
		}	
		if(end > len)
			end = len;

	 	for(i = 0; i < tap; i++)
			fputc(' ', stream);
		fprintf(stream, "%6u ", start);
		for(i = start - 1; i < end; i++) {
			fputc(upper[i], stream);
		}
		fprintf(stream, " %-6u\n", end);

		// print genomic sequences	
	 	for(i = 0; i < tap; i++)
			fputc(' ', stream);
		fprintf(stream, "       ");
		for(i = start - 1; i < end; i++) {
			fputc(gen_seq[i], stream);
		}
		fputc('\n', stream);

		if(flag & PRINT_AA) {
			// print amino acid sequences
			for(i = 0; i < tap; i++)
				fputc(' ', stream);
			fprintf(stream, "       ");
			for(i = start - 1; i < end; i++) {
				fputc(aa_seq[i], stream);
			}
			fputc('\n', stream);
		}

		fputc('\n', stream);

		cur = end + 1;
	}
	fputc('\n', stream);

	free_cnt(upper);
	free_cnt(aa_seq);	

}


// Find location of the first appeared termination codon in the input sequence. 
// Parameters
// *********************************************************************************
//   transcript *trsc: 	Transcript struct. Used for get basic information
//   char *gen_seq:	transcript (cDNA) sequence. Maybe modified by INDEL
//   unsigned int cdna_start: starting point for PTC search
// *********************************************************************************
// Return value: Coordinate of the first Stop codon
// *********************************************************************************
unsigned int get_stop_codon_position_from_seq(transcript* trsc, char* gen_seq, unsigned int cdna_start) {
	unsigned int i, j, len = strlen(gen_seq), term = 0;
	int code = 0;
	
	// make coordinates to 0-base
	cdna_start--;

	// Check the start codon
	if(((gen_seq[cdna_start]   == 'a') || (gen_seq[cdna_start]   == 'A')) &&
	   ((gen_seq[cdna_start+1] == 't') || (gen_seq[cdna_start+1] == 'T')) && 
	   ((gen_seq[cdna_start+2] == 'g') || (gen_seq[cdna_start+2] == 'G'))) {
		// Okay
	}
	else {
		fprintf(stderr, "[ERROR] Wrong cDNA_start: Sequence at the cdna_start is not \"ATG\"\n");
		fprintf(stderr, "[ERROR] Transcript: ENST%011d\n", trsc->ensembl_id);
		fprintf(stderr, "[ERROR] Input sequence: %s\n", gen_seq);
		fprintf(stderr, "[ERROR] Input cdna_start: %u\n", cdna_start);
		return 0;
	}

	// find matched amino acid sequences
	for(i = 0; i < len; ) {
		if((!term) && (i >= cdna_start) && (i+2 < len)) {
			code = 0;
			for(j = 0; j < 3; j++) {
				code = code << 2;
				if((gen_seq[i+j] == 'T') || (gen_seq[i+j] == 't')) code = code | 0;
				else if((gen_seq[i+j] == 'C') || (gen_seq[i+j] == 'c')) code = code | 1;
				else if((gen_seq[i+j] == 'A') || (gen_seq[i+j] == 'a')) code = code | 2;
				else if((gen_seq[i+j] == 'G') || (gen_seq[i+j] == 'g')) code = code | 3;
				else {
					code = 999;   // when N is in the codon 
					break;
				}
			}
			//printf("seq = %c%c%c, code = %d, ", gen_seq[i], gen_seq[i+1], gen_seq[i+2], code);
			if(code < 64) {
//				aa_seq[i]   = codon_tab[code][0];
//				aa_seq[i+1] = codon_tab[code][1];
//				aa_seq[i+2] = codon_tab[code][2];
				if(strcmp(codon_tab[code], "***") == 0) {
					term = i;
				}
			}
			else {
//				aa_seq[i]   = ' ';
//				aa_seq[i+1] = ' ';
//				aa_seq[i+2] = ' ';
			}
			//printf("aa seq = %c%c%c\n", aa_seq[i], aa_seq[i+1], aa_seq[i+2]);

			i = i + 3;
		}
		else {
//			aa_seq[i]   = ' ';
			i++;
		}
	}

	return (term + 1); // Return 1-base DNA coordinate of stop codon
}


// Generate DNA sequece of given transcript (trsc), then print its sequence with coded AA sequence
void print_transcript(genome* gen, transcript* trsc, int flag, FILE* stream, int hl_start, int hl_end, int tap) {
	unsigned int cdna_start = trsc->cdna_start + 3;	// 3bp of buffer will be added to the concatenated exon seq.
	char *merged_seq = concat_exon_seqs(gen, trsc, ADD_3BP_BUFFER);

	if(((exon*)(trsc->first_exon->entry))->start_phase > 0 )
		cdna_start -= ((exon*)(trsc->first_exon->entry))->start_phase;

	print_seq(merged_seq, cdna_start, trsc->cdna_end, flag, stream, hl_start, hl_end, tap);		

	free_cnt(merged_seq);
	return;
}


// Create a new transcript struct from original, without the skipped_exon_idx th exon.
int create_skipped_transcript(genome* gen, transcript* original, transcript* skipped, int skipped_exon_idx, FILE* stream) {
	int flag = 0, i;
	lentry *cursor_orig, *cursor_exon, *cursor_coding, *cursor_orig_coding;
	exon *cur_exon;
	unsigned int exon_start, exon_end, exon_length;	

	// exon list should be copyed, not referenced, because one of them will be removed
	memcpy(skipped, original, sizeof(transcript));
	skipped->first_exon = NULL;
	skipped->is_coding_exon = NULL;
	cursor_orig = original->first_exon;
	cursor_orig_coding = original->is_coding_exon;

	cursor_exon = NULL;
	cursor_coding = NULL;
	cur_exon = (exon*)cursor_orig->entry;
	exon_start = 1;
	exon_end = exon_start + (cur_exon->end - cur_exon->start);
	exon_length = cur_exon->end - cur_exon->start + 1;

	for(i = 0; i < original->num_exons; i++) {
		if(skipped_exon_idx == i) {
			if(stream != NULL) {
				fprintf(stream, "Exon %s%011d (%d th exon) has ESE(s). ", gen->prefix_exon, cur_exon->ensembl_id, skipped_exon_idx + 1);
			}
			// Adjust coding region information

			// cur transcript doesn't have any coding region --> nothing to do
			if((skipped->cdna_start == 0) && (skipped->cdna_end == 0)) {
				// To find ORF, uncomment below lines
				//skipped->cdna_start = 1;
				//skipped->trsc_length -= exon_length;
				//skipped->cdna_end = skipped->trsc_length;			
			}
			//[cur exon]  <-----cDNA----->
			//target exon is on 5'UTR
			else if((exon_start < skipped->cdna_start) && (exon_end < skipped->cdna_start)) {
				// adjust position and length of cDNA in transcript
				skipped->cdna_start -= exon_length;
				skipped->cdna_end -= exon_length;
				skipped->trsc_length -= exon_length;
			}
			//[cur exon]
			//      <-----cDNA----->
			// The start codon has been damaged / removed
			else if((exon_start <= skipped->cdna_start) && (skipped->cdna_start <= exon_end)) {
				lentry* inner_cursor;
				exon* inner_exon;

				// adjust cDNA start postition
				inner_cursor = original->first_exon;
				skipped->cdna_start = 1;
				while(inner_cursor != NULL) {
					inner_exon = (exon*)inner_cursor->entry;
					if(inner_exon != cur_exon)
						break;
					skipped->cdna_start += inner_exon->end - inner_exon->start + 1;
					inner_cursor = inner_cursor->next;
				}

				flag = flag | FIND_ATG;				// start codon has damaged
				skipped->trsc_length -= exon_length;	
			}
			//         [cur exon]
			//      <-----cdna----->
			else if((skipped->cdna_start <= exon_start) && (exon_end <= skipped->cdna_end)) {
				if(exon_start - skipped->cdna_start <= 2) {
					flag = flag | FIND_ATG;			// check if start codon has damaged
					skipped->cdna_start = exon_start;
				}
				skipped->trsc_length -= exon_length;
			}
			//                 [cur exon]
			//      <-----cdna----->
			else if((exon_start <= skipped->cdna_end) && (skipped->cdna_end < exon_end)) {
				if(exon_start - skipped->cdna_start <= 2) {
					flag = flag | FIND_ATG;			// check if start codon has damaged
					skipped->cdna_start = exon_start;
				}

				skipped->trsc_length -= exon_length;
			}
			//      <-----cdna----->     [cur exon]
			// target exon is on 3'UTR
			else if((skipped->cdna_end < exon_start) && (skipped->cdna_end < exon_end)) {
				// cdna_start is unchanged
				skipped->trsc_length -= exon_length;
			}
		
			// go to the next exon (target exon is skipped)
			cursor_orig = cursor_orig->next;
			cursor_orig_coding = cursor_orig_coding->next;
			i++;
		}

		if(cursor_orig != NULL) {
			// Create next entries and move the cursors
			if(cursor_exon == NULL) {
				skipped->first_exon = (lentry*)malloc_cnt(sizeof(lentry));
				skipped->is_coding_exon = (lentry*)malloc_cnt(sizeof(lentry));
				cursor_exon = skipped->first_exon;
				cursor_coding = skipped->is_coding_exon;
			}
			else {
				cursor_exon->next = (lentry*)malloc_cnt(sizeof(lentry));
				cursor_exon = cursor_exon->next;
				cursor_coding->next = (lentry*)malloc_cnt(sizeof(lentry));
				cursor_coding = cursor_coding->next;
			}
			cursor_exon->entry = cursor_orig->entry;
			cursor_exon->next = NULL;
			cursor_coding->entry = cursor_orig_coding->entry;
			cursor_coding->next = NULL;

			if(cursor_orig->next != NULL) {
				cursor_orig = cursor_orig->next;
				cur_exon = (exon*)cursor_orig->entry;
				exon_start = exon_end + 1;
				exon_end = exon_start + (cur_exon->end - cur_exon->start);
				exon_length = cur_exon->end - cur_exon->start + 1;
				cursor_orig_coding = cursor_orig_coding->next;
			}
		}
	}

	skipped->num_exons--;
	skipped->aa_length = calculate_aa_len(gen, skipped, flag);
	
	return flag;
}

// Because the skipped transcript create own lentry struct for exon, 
// the lengty struct should be freed. 
// For other, only addresses are copied.
void free_skipped_transcript(transcript* skipped) {
	lentry* cursor_exon;

	while(skipped->first_exon != NULL) {
		cursor_exon = skipped->first_exon;
		skipped->first_exon = skipped->first_exon->next;
		free_cnt(cursor_exon);
	}
	free_cnt(skipped);
}

