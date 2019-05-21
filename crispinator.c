/*
 * crispinator.c     CRISPinatoR
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
#include <math.h>
#include "tools.h"

/*************************************************************************************
 * A main function that write sgRNAs for a gene.
 * Parameters
 * ***********************************************************************************
 * genome *gen:		genome struct
 * kmers *ese:		kmer struct storing ESE
 * bwt_t *bwt, bntseq_t *bns, gap_opt_t *opt:	bwt structs required for off-target searching
 * gene *target_gene:	gene struct of the target gene. All sgRNA for this gene will be found
 * FILE *result:	FILE struct for output file
 * FILE *log:		FILE struct for log file. It can be NULL
 * double cut_offs:	Threshold for off-target score. sgRNAs having higher score than this will be written
 * int cut_5:		Maximum upstream distance from cut site to any ESE
 * int cut_3:		Maximum downstream distance from cut site to any ESE
 * int cut_adjust:	Assumed cut site from PAM sequence. If -3, 'A' of ABCNGG will be assumed as a cut site
 * int verbose:
 * ***********************************************************************************
 * Return : void
 * ***********************************************************************************/
void write_sgrna_gene(genome *gen, kmers *ese, bwt_t *bwt, bntseq_t *bns, gap_opt_t *opt, gene* target_gene, FILE *result, FILE* log, double cut_offs, int cut_5, int cut_3, int cut_adjust, int verbose) {
	char *exon_seq = NULL, *ese_seq = NULL, *gen_seq = NULL, *rna_only = NULL, *trsc_seq, *rna_padded = NULL;
	int i, j, k, l, is_coding, longest_aa_len, longest_dna_len, larger_padded = (PADDED_5 >= PADDED_3) ? PADDED_5 : PADDED_3, target_exon_idx_total, target_exon_idx_coding, dist_ptc_lee_skip;
	unsigned int start, end, len, pam_pos, cut_site, pos_cut_cdna, pos_cut_cds, pos_last_ee, pos_ptc_fs1, pos_ptc_fs2, pos_ptc_skip, pos_last_ee_when_skip;
	double cut_site_ratio;
	lentry *cursor, *cursor_exon, *cur_coding;
	exon* cur_exon, *prev_exon, *target_exon;
	transcript** original = NULL, *lc_trsc = NULL, *skipped;
	crispr_candidate* candidates = NULL;
	prev_exon = cur_exon = NULL;

	// variables for alignment
	bwa_seq_t *seqs, *p;
	int n_cand, idx_cand, seqid;
	

	// ************************************
	// STEP 1: Find sgRNA candidates
	// ************************************
	if(verbose)
		printf("    Gene info: %s:%u-%u\n", target_gene->chr->name, target_gene->start, target_gene->end);

	if(log != NULL) {
		fprintf(log, "Gene %s (%s%011d)\n", target_gene->symbol, gen->prefix_gene, target_gene->ensembl_id);
		fprintf(log, "Gene info: %s:%u-%u\n", target_gene->chr->name, target_gene->start, target_gene->end);
	}

	// copy pointers to transcripts.
	// AA lengths are counted and AA lengths for NMD transcripts are set to 0
	// Also pick the longest coding transcript of the target gene
	if(target_gene->num_transcripts > 0) {
		original = (transcript**)malloc_cnt(sizeof(transcript*) * target_gene->num_transcripts);
		longest_aa_len = longest_dna_len = 0;
		j = 0;
		FOR_LENTRY(cursor, target_gene->first_trsc) {
			// Store transcripts of the target gene, and calculate AA length.
			original[j] = (transcript*)cursor->entry;
			if( (!strcmp(original[j]->trsc_type, "nonsense_mediated_decay"))
			 || (!strcmp(original[j]->trsc_type, "non_stop_decay")) 
			 || (!strcmp(original[j]->trsc_type, "polymorphic_pseudogene")) ){
				original[j]->aa_length = 0;
			}
			else if(original[j]->ensembl_protein_id > 0) {
				if(strcmp(original[j]->trsc_type, "protein_coding")) {
					printf("Unexpected transcript type: %s, but processed anyway.\n", original[j]->trsc_type);
				}

				original[j]->aa_length = calculate_aa_len(gen, original[j], 0);
				target_gene->n_functional_trsc++;
			}
			else 
				original[j]->aa_length = 0;

			// The longest coding transcript should coding transcript, and should have proper start and stop codon
			if((original[j]->aa_length > 0) && (original[j]->has_CDS)) {
				if(original[j]->aa_length > longest_aa_len) {
					longest_aa_len = original[j]->aa_length;
					longest_dna_len = (original[j]->cdna_end - original[j]->cdna_start + 1);
					lc_trsc = original[j];
				}
				else if(( original[j]->aa_length == longest_aa_len ) &&
					((original[j]->cdna_end - original[j]->cdna_start + 1) > longest_dna_len) ) {
					longest_aa_len = original[j]->aa_length;
					longest_dna_len = (original[j]->cdna_end - original[j]->cdna_start + 1);
					lc_trsc = original[j];
				}
			}

			j++;
		}
	}
	else {
		if(log != NULL)
			fprintf(log, "[END] target gene does not have any transcript.\n"); 
		return;
	}

	
	if(lc_trsc == NULL) {
		if(log != NULL)
			fprintf(log, "[END] target gene does not have proper coding transcript.\n");
		return;
	}

	gen_seq  = get_subseq_by_idx(gen, target_gene->chr->idx, target_gene->start, target_gene->end);
	if(gen_seq == NULL) {
		if(log != NULL)
			fprintf(log, "Genomic range is out of chromosome\n");		
		printf("returned\n");
		return;
	}

	len = strlen(gen_seq);
	make_str_hc(gen_seq);
	if(target_gene->strand == REVERSE) {
		if(verbose)
			printf("    Converting seq. to reverse complement..\n");
		gen_seq = create_rc_free_exist(gen_seq);
	}

	if(verbose)
		printf("    Counting candidate..\n");

	// Ccount number of candidate sgRNAs
	// Number of input sequences should be pre-counted for bwa.
	n_cand = 0;
	for(i = SGRNA_LEN + larger_padded - 1; i < len - larger_padded; i++) {
		if(((gen_seq[i-1] == 'G') && (gen_seq[i] == 'G')) ||
		   ((gen_seq[i-1] == 'C') && (gen_seq[i] == 'C'))) {
			// pam_pos is N's genomic coordinate of the NGG 
			// start & end are genomic coordinates of search range for this sgRNA
			if(target_gene->strand == FORWARD) {
				if(gen_seq[i] == 'G') {	// Gene is forward, NGG gRNA (sense)
					pam_pos = target_gene->start + i - 2;
					cut_site = pam_pos + cut_adjust;
					start = cut_site - cut_5;
					end   = cut_site + cut_3 - 1;
				}
				else {	// Forward gene, CCN gRNA (antisense)
					pam_pos = target_gene->start + i + 1;
					cut_site = pam_pos - cut_adjust;
					start = cut_site - cut_3 + 1;
					end   = cut_site + cut_5;
				}
			}
			else {
				if(gen_seq[i] == 'G') {	// Reverse gene, NGG gRNA (sense)
					pam_pos = target_gene->end - i + 2;
					cut_site = pam_pos - cut_adjust;
					start = cut_site - cut_3 + 1;
					end   = cut_site + cut_5;
				}
				else {	// Reverse gene, CCN gRNA (antisense)
					pam_pos = target_gene->end - i - 1;
					cut_site = pam_pos + cut_adjust;
					start = cut_site - cut_5;
					end   = cut_site + cut_3 - 1;
				}
			}

			is_coding = 0;
			cur_coding = lc_trsc->is_coding_exon;
			// Traverse all exons of the lc_trsc (longest coding transcript) of the gene
			FOR_LENTRY(cursor_exon, lc_trsc->first_exon) {
				cur_exon = (exon*)cursor_exon->entry;
				// the exon is a coding exon
				if(cur_coding->entry != NULL) {
					// the search range should be inside of the coding exon, 
					//   to avoid editing exon boundary
					if(((cur_exon->start <= start) && (start <= cur_exon->end)) &&
					   ((cur_exon->start <= end) && (end <= cur_exon->end))) {
						is_coding = 1;
						break;
					}
				}
				cur_coding = cur_coding->next;
			}

			// Convert genomic coordinate of cut_site as transcript coordiante
			pos_cut_cdna = convert_genomic_coord_to_trsc_coord(lc_trsc, cut_site);

			if(pos_cut_cdna > 0) { 	// Convert DNA-level CDS postion to AA level. The function returned 1-based coordiate. So the 1st-3rd nt should be matched to 1st AA.
				if(pos_cut_cdna >= lc_trsc->cdna_start) 
					pos_cut_cds = (pos_cut_cdna - lc_trsc->cdna_start) / 3 + 1;
				else
					pos_cut_cds = 0;	// PAM is on UTR
			}
			else {
				pos_cut_cds = 0;	// PAM is outside of transcript
			}

			// if the sgRNA targets at least one coding exon, count it as a candidate
			if((is_coding) && (pos_cut_cds > 0)) {
				n_cand++;
			}
		}
	}

	if(n_cand == 0) {
		free_cnt(gen_seq);
		free_cnt(original);
		if(log != NULL)
			fprintf(log, "[END] No candidate.\n");
		return;
	}

	if(verbose)
		printf("    Generating %d candidates..\n", n_cand);

	if(log != NULL)
		fprintf(log, " Generating %d candidates..\n", n_cand);

	
	// ************************************
	// STEP 2: Makes read structs for bwa. 
	// In step 2~4, many codes are from bwa source code
	// ************************************

	// Create candidate struct, containing information about sgRNA and alignments
	candidates = (crispr_candidate*)malloc_cnt(sizeof(crispr_candidate) * n_cand);
	for(i = 0; i < n_cand; i++) {
		candidates[i].target_gene = target_gene;
		candidates[i].n_in_gene = (int*)malloc_cnt(sizeof(int) * (opt->max_diff + 1));
		for(j = 0; j <= opt->max_diff; j++) {
			candidates[i].n_in_gene[j] = 0;
		}
		candidates[i].pos_in_chr = 0;
		candidates[i].pam_pos = 0;
		candidates[i].n_ese = 0;
		candidates[i].n_aff_trsc = 0;
		candidates[i].n_aff_functional_trsc = 0;
		candidates[i].s_pos = 0.0;
		candidates[i].s_off = 0.0;
		candidates[i].s_trsc = 0.0;
		candidates[i].s_ese = 0.0;
		candidates[i].s_sum = 0.0;
		candidates[i].dist5_ese = -1;
		candidates[i].dist3_ese = -1;
		candidates[i].rna_seq = (char*)malloc_cnt(sizeof(char) * (SGRNA_LEN + 1));
		candidates[i].rna_seq_padded = (char*)malloc_cnt(sizeof(char) * (SGRNA_LEN + PADDED_5 + PADDED_3 + 1));
		candidates[i].rna_seq[0] = '\0';
		candidates[i].rna_seq_padded[0] = '\0';
	}

	// Temperory string buffers
	ese_seq = (char*)malloc(sizeof(char) * 102400);
	rna_only = (char*)malloc(sizeof(char) * SGRNA_LEN);
	rna_padded = (char*)malloc(sizeof(char) * (SGRNA_LEN + PADDED_5 + PADDED_3 + 1));
	ese_seq[0] = '\0';

	idx_cand = 0;
	// Prepare bwa input sequences
	seqs = (bwa_seq_t*)malloc_cnt(n_cand * sizeof(bwa_seq_t));

	// Search gene sequence again, to create the sgRNA.
	for(i = SGRNA_LEN + larger_padded - 1; i < len - larger_padded; i++) {
		if(((gen_seq[i-1] == 'G') && (gen_seq[i] == 'G')) ||
		   ((gen_seq[i-1] == 'C') && (gen_seq[i] == 'C'))) {
			if(n_cand == idx_cand)
				break;
			// pam_pos is chrosomal position of 'N'GG 
			// start & end are genomic coordinates of search range for this sgRNA
			if(target_gene->strand == FORWARD) {
				if(gen_seq[i] == 'G') {  // ####################ngG on the reference genome. i points the uppercase G
					// store informations for each candidate
					candidates[idx_cand].pam_pos = target_gene->start + i - 2;
					candidates[idx_cand].cut_site = candidates[idx_cand].pam_pos + cut_adjust;
					candidates[idx_cand].pos_in_chr = target_gene->start + i - SGRNA_LEN + 1;	// Genomic coordinate of PAM
					candidates[idx_cand].pos_in_gene = i - SGRNA_LEN + 1;	// Coordinate in gene sequence. Will be used to substr() or strcpy()
					start = candidates[idx_cand].cut_site - cut_5;
					end   = candidates[idx_cand].cut_site + cut_3 - 1;
					candidates[idx_cand].pam_direction = PAM_NGG;	// Direction of PAM (sense / antisense)

					if(log != NULL) {
						fprintf(log, "Case FORWARD ###NGG. i = %d, pam_pos = %u, pos_in_chr = %u, pos_in_gene = %u", i, candidates[idx_cand].pam_pos, candidates[idx_cand].pos_in_chr, candidates[idx_cand].pos_in_gene);
						fprintf(log, ", start-end = %u-%u\n", start, end);
					}
				}
				else {			//  cCn#################### on the reference genome. i points the uppercase C
					candidates[idx_cand].pam_pos = target_gene->start + i + 1;
					candidates[idx_cand].cut_site = candidates[idx_cand].pam_pos - cut_adjust;
					candidates[idx_cand].pos_in_chr = target_gene->start + i - 1;
					candidates[idx_cand].pos_in_gene = i - 1;
					start = candidates[idx_cand].cut_site - cut_3 + 1;
					end   = candidates[idx_cand].cut_site + cut_5;
					candidates[idx_cand].pam_direction = PAM_CCN;
					if(log != NULL) {
						fprintf(log, "Case FORWARD CCN###. i = %d, pam_pos = %u, pos_in_chr = %u, pos_in_gene = %u", i, candidates[idx_cand].pam_pos, candidates[idx_cand].pos_in_chr, candidates[idx_cand].pos_in_gene);
						fprintf(log, ", start-end = %u-%u\n", start, end);
					}
				}
			}
			else {
				if(gen_seq[i] == 'G') {	// Ggn#################### (reverse, view on reference genome). i points the uppercase G
					candidates[idx_cand].pam_pos = target_gene->end - i + 2;
					candidates[idx_cand].cut_site = candidates[idx_cand].pam_pos - cut_adjust;
					candidates[idx_cand].pos_in_chr = target_gene->end - i;
					candidates[idx_cand].pos_in_gene = i - SGRNA_LEN + 1;
					start = candidates[idx_cand].cut_site - cut_3 + 1;	// Genomic coordinates
					end   = candidates[idx_cand].cut_site + cut_5;
					candidates[idx_cand].pam_direction = PAM_NGG;
					if(log != NULL) {
						fprintf(log, "Case REVERSE ###NGG. i = %d, pam_pos = %u, pos_in_chr = %u, pos_in_gene = %u", i, candidates[idx_cand].pam_pos, candidates[idx_cand].pos_in_chr, candidates[idx_cand].pos_in_gene);
						fprintf(log, ", start-end = %u-%u\n", start, end);
					}
				}
				else {			// ####################nCc (reverse, view on reference genome). i points the uppercase C.
					candidates[idx_cand].pam_pos = target_gene->end - i - 1;
					candidates[idx_cand].cut_site = candidates[idx_cand].pam_pos + cut_adjust;
					candidates[idx_cand].pos_in_chr = target_gene->end - i - SGRNA_LEN + 2;
					candidates[idx_cand].pos_in_gene = i - 1;
					start = candidates[idx_cand].cut_site - cut_5;
					end   = candidates[idx_cand].cut_site + cut_3 - 1;
					candidates[idx_cand].pam_direction = PAM_CCN;
					if(log != NULL) {
						fprintf(log, "Case REVERSE CCN###. i = %d, pam_pos = %u, pos_in_chr = %u, pos_in_gene = %u", i, candidates[idx_cand].pam_pos, candidates[idx_cand].pos_in_chr, candidates[idx_cand].pos_in_gene);
						fprintf(log, "start-end = %u-%u\n", start, end);
					}
				}
			}

			// Find a target exon on the longest coding transcript
			target_exon = NULL;
			cur_coding = lc_trsc->is_coding_exon;
			FOR_LENTRY(cursor_exon, lc_trsc->first_exon) {
				cur_exon = (exon*)cursor_exon->entry;
				// the search range should be inside of the coding exon, 
				//   to avoid editing exon boundary
				if(((cur_exon->start <= start) && (start <= cur_exon->end)) &&
				   ((cur_exon->start <= end) && (end <= cur_exon->end))) {
					// if cur_exon is a coding exon		
					if(cur_coding->entry != NULL) {
						target_exon = cur_exon;
						break;
					}
				}
				cur_coding = cur_coding->next;
			}	
	
			// Convert genomic coordinate of cut_site as transcript coordiante
			pos_cut_cdna = convert_genomic_coord_to_trsc_coord(lc_trsc, candidates[idx_cand].cut_site);

			if(pos_cut_cdna > 0) { 	// Convert DNA-level CDS postion to AA level. The function returned 1-based coordiate. So the 1st-3rd nt should be matched to 1st AA.
				if(pos_cut_cdna >= lc_trsc->cdna_start) 
					pos_cut_cds = (pos_cut_cdna - lc_trsc->cdna_start) / 3 + 1;
				else
					pos_cut_cds = 0;	// PAM is on UTR
			}
			else {
				pos_cut_cds = 0;	// PAM is outside of transcript
			}
			
			// the sgRNA can be a candidate
			if((pos_cut_cds > 0) && (target_exon != NULL)) {
				// prepare sgRNA seq as a string.
				if(candidates[idx_cand].pam_direction == PAM_NGG) {
					strncpy(candidates[idx_cand].rna_seq, (gen_seq + candidates[idx_cand].pos_in_gene), SGRNA_LEN);	
					strncpy(candidates[idx_cand].rna_seq_padded, (gen_seq + candidates[idx_cand].pos_in_gene - PADDED_5), SGRNA_LEN + PADDED_5 + PADDED_3);
				}
				else {
					strncpy_rc(candidates[idx_cand].rna_seq, (gen_seq + candidates[idx_cand].pos_in_gene), SGRNA_LEN);
					strncpy_rc(candidates[idx_cand].rna_seq_padded, (gen_seq + candidates[idx_cand].pos_in_gene - PADDED_3), SGRNA_LEN + PADDED_5 + PADDED_3);
				}

				// gRNA sequences are prepared (with/without padded sequence)
				candidates[idx_cand].rna_seq[SGRNA_LEN] = '\0';
				candidates[idx_cand].rna_seq_padded[SGRNA_LEN + PADDED_5 + PADDED_3] = '\0';
				if(log != NULL) {
					fprintf(log, "sgRNA seq = %s\n", candidates[idx_cand].rna_seq);
					fprintf(log, "sgRNA seq with padding = %s\n", candidates[idx_cand].rna_seq_padded);
				}

				// extract exon sequence and mark ESE on the exon (keep previous exon and skip this if gRNA is on the same exon)
				if( (prev_exon == NULL) || 
				    ((prev_exon != NULL) && (prev_exon != target_exon)) ) {
					if(exon_seq != NULL)
						free_cnt(exon_seq);

					// get exon sequence
					exon_seq = get_subseq_by_idx(gen, target_exon->chr->idx, target_exon->start, target_exon->end);
					make_str_hc(exon_seq);
					if(target_gene->strand == REVERSE)
						exon_seq = create_rc_free_exist(exon_seq);

					for(j = 0; j < strlen(exon_seq); ++j)
						ese_seq[j] = ' ';
					ese_seq[strlen(exon_seq)] = '\0';

					// mark ESE as * on the exon sequence
					for(k = 0; k < ese->length; k++) {
						for(j = 0; j < strlen(exon_seq) - ese->k_size[k] + 1; j++) {
							if(strncmp(exon_seq + j, ese->mers[k], ese->k_size[k]) == 0) {
								for(l = 0; l < ese->k_size[k]; l++)
									ese_seq[j + l] = '*';
							}
						}
					}
				}

				// Position of sgRNA on the exon
				unsigned int pos_in_exon = (target_gene->strand > 0) ? candidates[idx_cand].cut_site - target_exon->start : target_exon->end - candidates[idx_cand].cut_site;

				if(log != NULL) {
					for(j = 0; j < pos_in_exon + 10; j++) 
						fprintf(log, " ");
					fprintf(log, "*\n");
					fprintf(log, "Exon seq: %s (%d)\n", exon_seq, (int)strlen(exon_seq));
					fprintf(log, "ESE seq : %s (%d)\n", ese_seq,  (int)strlen(ese_seq));
				}

				// check distance from PAM to closest ESE on 5' & 3' sides, using the ese_seq
				j = 0;
				while(pos_in_exon >= j) {
					if(ese_seq[pos_in_exon - j] == '*') {	// ESE was marked as *
						if(candidates[idx_cand].pam_direction == PAM_NGG)
							candidates[idx_cand].dist5_ese = j;
						else
							candidates[idx_cand].dist3_ese = j;
						break;
					}
					j++;
				}
				j = 0;
				while(pos_in_exon + j < strlen(ese_seq)) {
					if(ese_seq[pos_in_exon + j] == '*') {
						if(candidates[idx_cand].pam_direction == PAM_NGG)
							candidates[idx_cand].dist3_ese = j;
						else
							candidates[idx_cand].dist5_ese = j;
						break;
					}
					j++;
				}

				if(log != NULL)
					fprintf(log, "Pos. in exon = %u, dist to 5/3 = %d/%d\n", pos_in_exon, candidates[idx_cand].dist5_ese, candidates[idx_cand].dist3_ese);

				// if there is no ESE in boundary cut. To uncomment it, n_cand should be modified....
				//if( ((candidates[idx_cand].dist5_ese < 0) && (candidates[idx_cand].dist3_ese < 0)) ||
				//    ((candidates[idx_cand].dist5_ese > cut_5) && (candidates[idx_cand]dist3_ese > cut_3))) {
				//		candidates[idx_cand].dist5_ese = -1;
				//		candidates[idx_cand].dist3_ese = -1;
				//		continue;
				//}

				// create read sequence for off-target search using bwa
				// source code below is from bwa
				p = &seqs[idx_cand];
				candidates[idx_cand].aln = p;
				init_seq(p, SGRNA_LEN);

				for (j = 0; j != p->full_len; ++j) {
					p->seq[j] = nst_nt4_table[(int)candidates[idx_cand].rna_seq[j]];
					p->qual[j] = 'h';
				}

				memcpy(p->rseq, p->seq, p->len);
				seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
				seq_reverse(p->len, p->rseq, 0);
				p->name = (char*)malloc_cnt(128 * sizeof(char));
				sprintf(p->name, "@candidate_%d (pos of PAM on chr: %u)", idx_cand + 1, candidates[idx_cand].pam_pos);
				// bwa code ends
				
				candidates[idx_cand].target_exon = target_exon;
				candidates[idx_cand].target_lct = lc_trsc;

				idx_cand++;

				prev_exon = target_exon;
			}
		}
	}

	if(verbose)
		printf("    Candidates were prepared.\n");

	if(idx_cand != n_cand) {
		fprintf(stderr, "idx_cand = %d, n_cand = %d\n", idx_cand, n_cand);
	}

	// ************************************
	// STEP 3: Do alignment
	// ************************************
	if(verbose)
		printf("    Start alignment (%d candidates)\n", idx_cand);
	do_alignment(gen, bwt, bns, candidates, seqs, idx_cand, opt);
	if(verbose)
		printf("    Finished\n");

	if(log != NULL)
		fprintf(log, "Alignment result: \n");

	// ************************************
	// STEP 4: Print details for each candidates
	// ************************************
	for(i = 0; i < idx_cand; i++) {
		p = candidates[i].aln;

		// convert each alignment position to readable number
		for(j = 0; j < p->n_multi; j++) {
			bwt_multi1_t *q = p->multi + j;	
			bns_cnt_ambi(bns, q->pos, p->len, &seqid);

			if((q->mm == 0) && ((unsigned int)(q->pos - bns->anns[seqid].offset + 1) == candidates[i].pos_in_chr)) {
				for(k = j + 1; k < p->n_multi; ++k)
					p->multi[k - 1] = p->multi[k]; 
				p->n_multi--;
				break;
			}
		}

		// Calculate off-target score
		candidates[i].s_off = calculate_fs_offtarget_score(gen, &candidates[i], bns, NULL);

		/**** header, FYI ****
		fprintf(result, "Gene.ID\tSymbol\tChr\tG.start\tG.end\t");					// 0 1 2 3 4
		fprintf(result, "Num.trscs\tLCT.ID\tT.exon.ID\tSym\tTrsc.info\t");				// 5 6 7 8 9
		fprintf(result, "Num.CT.with.TE\tNum.total.CT\tNum.NCT.with.TE\tNum.total.NCT\tgRNA.seq\t");	// 10 11 12 13 14
		fprintf(result, "gRNA.PAM.seq\tgRNA.PAM.pad.seq\tSense\tPos.PAM\t");		// 15 16 18 19
		fprintf(result, "Pos.CUT.CDS(PAM%+d)\tLTC.CDS.len\tPos.PAM.CDS.Ratio\tNo.T.exon.LCT\tNum.CE.in.LCT\t", cut_adjust);	// 20 21 22 23 24
		fprintf(result, "Exon.Loc\tPos.LEE\tPos.PTC.fs1\tPos.PTC.fs2\tPos.PTC.skip\t");			// 25 26 27 28 29
		fprintf(result, "Pos.LEE.skip\tDist.PTC.LEE.skip\tDist.ESE.5\tDist.ESE.3\tDist.ESE.min\t");	// 30 31 32 33 34
		fprintf(result, "OffT.score\tOffT.list\tOffT.detail.list\n");					// 35 36 37
		*/

		if(candidates[i].s_off >= cut_offs) {
			int sym, affected_coding_trsc, affected_noncoding_trsc;

			cur_exon = candidates[i].target_exon;
			lc_trsc = candidates[i].target_lct;

			if((cur_exon->start_phase >= 0) && (cur_exon->start_phase == cur_exon->end_phase)) {
				sym = 1;
			}
			else if((cur_exon->start_phase == -1) && (cur_exon->end_phase == 0)) {
				sym = 1;
			}
			else if((cur_exon->start_phase == 0) && (cur_exon->end_phase == -1)) {
				sym = 1;
			}
			else {
				sym = 0;
			}
			// print 0-4
			//	          0        1   2   3   4        
			fprintf(result, "%s%011d\t%s\t%s\t%u\t%u\t", 
					gen->prefix_gene, target_gene->ensembl_id, 	// 0
				        target_gene->symbol, target_gene->chr->name, 	// 1 2
					target_gene->start, target_gene->end); 		// 3 4

			//		  5   6        7       8
			fprintf(result, "%d\t%s%011d\t%s%011d\t%s\t",
					target_gene->num_transcripts, gen->prefix_trsc, lc_trsc->ensembl_id,		  // 5 6 
					gen->prefix_exon, cur_exon->ensembl_id, (sym == 1) ? "Symmetric" : "Asymmetric"); // 7 8

			// print ordinal number of the target exon and count total number of affected coding / noncoding transcripts
			affected_noncoding_trsc = affected_coding_trsc = 0;

			// start-end region, potentially affedted by the candidate i
			j = (target_gene->strand == FORWARD) ? 1 : -1;
			if(candidates[i].pam_direction == PAM_CCN)
				j *= -1;

			if(j > 0) {
				start = candidates[i].cut_site - cut_5 - 1;
				end   = candidates[i].cut_site + cut_3 - 1;
			}
			else {
				start = candidates[i].cut_site - cut_3 + 1;
				end   = candidates[i].cut_site + cut_5 + 1;
			}

			// traverse all transcripts and write potential effect of sgRNA for each transcript
			target_exon_idx_total = target_exon_idx_coding = 0;
			for(j = 0; j < target_gene->num_transcripts; j++) {
				k = 0;
				l = 0;
				// is_coding_exon->entry points an exon if it is coding exon.
				// if it is not a coding exon, it points NULL
				cur_coding = original[j]->is_coding_exon;
				FOR_LENTRY(cursor_exon, original[j]->first_exon) {
					cur_exon = cursor_exon->entry;
					// If this exon is a coding exon
					if(cur_coding->entry != NULL) {
						k++;
					}
					l++;
					// k and l are 1 for the fisrt exon
					if(((cur_exon->start <= start) && (start <= cur_exon->end)) &&
					   ((cur_exon->start <= end)   && (end <= cur_exon->end))) {
						if((original[j] == lc_trsc) && (cur_exon == candidates[i].target_exon)) {
							target_exon_idx_total = l;
							target_exon_idx_coding = k;
						}

						if(!strcmp(original[j]->trsc_type, "nonsense_mediated_decay")) {
							fprintf(result, "-/%d[%d/%d](%s%011d-%s%011d:0aa) ", original[j]->num_coding_exons, l, original[j]->num_exons, gen->prefix_exon, cur_exon->ensembl_id, gen->prefix_trsc, original[j]->ensembl_id);
							affected_noncoding_trsc++;
						}
						else if(original[j]->ensembl_protein_id > 0) {
							fprintf(result, "%d/%d[%d/%d](%s%011d-%s%011d:%daa) ", k, original[j]->num_coding_exons, l, original[j]->num_exons, gen->prefix_exon, cur_exon->ensembl_id, gen->prefix_trsc, original[j]->ensembl_id, original[j]->aa_length);
							affected_coding_trsc++;
						}
						else {
							fprintf(result, "-/%d[%d/%d](%s%011d-%s%011d:0aa) ", original[j]->num_coding_exons, l, original[j]->num_exons, gen->prefix_exon, cur_exon->ensembl_id, gen->prefix_trsc, original[j]->ensembl_id);
							affected_noncoding_trsc++;
						}
						break;
					}
					cur_coding = cur_coding->next;
				}
			}

			if((target_exon_idx_total <= 0) || (target_exon_idx_coding <= 0)) {
				//fprintf(stderr, "[ERROR] Cannot find the ordinal number of target exon in the LCT!!\n");
				//fprintf(stderr, "[ERROR] %s%011d - %s%011d\n", gen->prefix_trsc, lc_trsc->ensembl_id, gen->prefix_exon, candidates[i].target_exon->ensembl_id);
				fprintf(result, "Edit exon boundary\n");
				continue;
			}

			strncpy(rna_only, candidates[i].rna_seq, SGRNA_LEN - 3);
			rna_only[SGRNA_LEN - 3] = '\0';

			//		   10  11  12  13  	
			fprintf(result, "\t%d\t%d\t%d\t%d\t", 
					affected_coding_trsc, target_gene->n_functional_trsc, //10, 11
					affected_noncoding_trsc, target_gene->num_transcripts - target_gene->n_functional_trsc); // 12, 13
			//		 14  15  16  17  
			fprintf(result, "%s\t%s\t%s\t%s\t", 
					rna_only, candidates[i].rna_seq, candidates[i].rna_seq_padded,		// 14 15 16
					//rna_only, candidates[i].rna_seq, candidates[i].rna_seq_padded, rna_padded,		// 14 15 16 17
					(candidates[i].pam_direction == PAM_NGG) ? "Sense" : "Antisense"); 	// 18

			// Calculate Cut.site.CDS (Pos.PAM + adjust_pos) (on AA sequence)
			pos_cut_cdna = convert_genomic_coord_to_trsc_coord(lc_trsc, candidates[i].cut_site);
			if(pos_cut_cdna > 0) { 	// Convert DNA-level CDS postion to AA level. The function returned 1-based coordiate. So the 1st-3rd nt should be matched to 1st AA.
				if(pos_cut_cdna >= lc_trsc->cdna_start) {
					if(pos_cut_cdna <= lc_trsc->cdna_end)
						pos_cut_cds = (pos_cut_cdna - lc_trsc->cdna_start) / 3 + 1;
					else
						pos_cut_cds = 0;	// cut site is on 3'UTR
				}
				else {
					pos_cut_cds = 0;	// cut site is on 5'UTR
				}
			}
			else {
				pos_cut_cds = 0;	// PAM is outside of transcript
			}

			cut_site_ratio = (double)pos_cut_cds / lc_trsc->aa_length;
			fprintf(result, "%u\t%u\t%u\t%1.4f\t%d\t%d\t%s\t", 
					candidates[i].pam_pos, pos_cut_cds, lc_trsc->aa_length, cut_site_ratio,	// 19 20 21 22
					target_exon_idx_coding, lc_trsc->num_coding_exons,	// 23 24
					(target_exon_idx_coding == 1) ? "First" : 
								        ((target_exon_idx_coding == lc_trsc->num_coding_exons) ? 
										"Last" : "Middle"));	// 25

			if(lc_trsc->num_coding_exons > 1) {
				// Traverse all coding exons
				FOR_LENTRY(cursor_exon, lc_trsc->is_coding_exon) {
					// If this exon is a coding exon, remember it
					if(cursor_exon->entry != NULL) {
						cur_exon = cursor_exon->entry;
					}
				}
				// cur_exon points the last coding exon. Calculate Pos.Last.EE.CDS (on AA sequence)
				if(target_gene->strand == FORWARD)
					pos_last_ee = convert_genomic_coord_to_trsc_coord(lc_trsc, cur_exon->start);
				else
					pos_last_ee = convert_genomic_coord_to_trsc_coord(lc_trsc, cur_exon->end);

				//if(pos_last_ee > 0) {
				//	pos_last_ee = (pos_last_ee - lc_trsc->cdna_start) / 3 + 1;
				//}
			}
			else {
				pos_last_ee = 0;	// Single exon genes.
			}

			// Predict altered AA sequence 
			if(pos_cut_cdna > lc_trsc->cdna_start + 2) {	// Check if PAM is before the coding region or on the start codon
				trsc_seq = concat_exon_seqs(gen, lc_trsc, NO_BUFFER);

				// Generate cDNA sequence with 1 nt frameshift at PAM
				for(j = pos_cut_cdna; j < strlen(trsc_seq); j++) {
					trsc_seq[j] = trsc_seq[j + 1];
				}
				
				// Predict location of PTC
				pos_ptc_fs1 = get_stop_codon_position_from_seq(lc_trsc, trsc_seq, lc_trsc->cdna_start);

				// Generate cDNA sequence with 2 nt frameshift at PAM
				for(j = pos_cut_cdna; j < strlen(trsc_seq); j++) { // 1 more framshift
					trsc_seq[j] = trsc_seq[j + 1];
				}
				// Predict location of PTC
				pos_ptc_fs2 = get_stop_codon_position_from_seq(lc_trsc, trsc_seq, lc_trsc->cdna_start);
				free_cnt(trsc_seq);
			}
			else {
				pos_ptc_fs1 = pos_ptc_fs2 = 0;
			}

			// Generate target exon skipped transcript
			// If the target exon is the first or the last coding exon, or the on the translation start site,
			if((convert_genomic_coord_to_trsc_coord(lc_trsc, candidates[i].target_exon->start) <= 3) ||
				(target_exon_idx_coding == 1) || 
				(target_exon_idx_coding == lc_trsc->num_coding_exons)) {
				pos_ptc_skip = 0;	// 	No exon skipping. 0 means NA
				pos_last_ee_when_skip = 0;
			}
			else {
				skipped = (transcript*)malloc_cnt(sizeof(transcript));
				create_skipped_transcript(gen, lc_trsc, skipped, target_exon_idx_total - 1, NULL);
				if(log != NULL) {
				}
				trsc_seq = concat_exon_seqs(gen, skipped, NO_BUFFER);
				pos_ptc_skip = get_stop_codon_position_from_seq(lc_trsc, trsc_seq, lc_trsc->cdna_start);
				free_cnt(trsc_seq);
			
				//if(skipped->num_exons > 1) {
				if(skipped->num_coding_exons > 1) {
					// Traverse all coding exons
					FOR_LENTRY(cursor_exon, skipped->is_coding_exon) {
						// If this exon is a coding exon, remember it
						if(cursor_exon->entry != NULL) {
							cur_exon = cursor_exon->entry;
						}
					}
					// cur_exon points the last exon. Calculate Pos.Last.EE (on cDNA sequence)
					// -- > cur_exon points the last coding exon. Calculate Pos.Last.EE.CDS (on AA sequence)
					if(target_gene->strand == FORWARD)
						pos_last_ee_when_skip = convert_genomic_coord_to_trsc_coord(skipped, cur_exon->start);
					else
						pos_last_ee_when_skip = convert_genomic_coord_to_trsc_coord(skipped, cur_exon->end);

					//if(pos_last_ee_when_skip > 0) {
					//	pos_last_ee_when_skip = (pos_last_ee_when_skip - skipped->cdna_start) / 3 + 1;
					//}
				}
				else {
					pos_last_ee_when_skip = 0;	// Single exon genes.
				}

				free_skipped_transcript(skipped);
			}

			dist_ptc_lee_skip = (int)pos_last_ee_when_skip - (int)pos_ptc_skip;
			//		 26  27  28  29  30  31
			fprintf(result, "%u\t%u\t%u\t%u\t%u\t%d\t",  
					pos_last_ee, pos_ptc_fs1, pos_ptc_fs2, pos_ptc_skip, 
					pos_last_ee_when_skip, dist_ptc_lee_skip);
		

			// No ESE on both side
			if((candidates[i].dist5_ese < 0) && (candidates[i].dist3_ese < 0)) {
				fprintf(result, "-\t-\t-\t");	// 32 33 34.  Print nothing
			}
			else if(candidates[i].dist5_ese < 0) {	// No ESE on 5' side
				fprintf(result, "-\t%d\t%d\t", candidates[i].dist3_ese, candidates[i].dist3_ese);	// 32 33 34. Print 3' side
			}
			else if(candidates[i].dist3_ese < 0) {	// No ESE on 3' side
				fprintf(result, "%d\t-\t%d\t", candidates[i].dist5_ese, candidates[i].dist5_ese);	// 32 33 34. Print 5' side
			}
			else {	// Else. Print closer ESE
				if(candidates[i].dist5_ese < candidates[i].dist3_ese) {
					fprintf(result, "%d\t%d\t%d\t", candidates[i].dist5_ese, 
									candidates[i].dist3_ese, 
									candidates[i].dist5_ese);	// 32 33 34
				}
				else {
					fprintf(result, "%d\t%d\t%d\t", candidates[i].dist5_ese,
									candidates[i].dist3_ese, 
									candidates[i].dist3_ese);	// 32 33 34
				}
			}

			//		 35
			fprintf(result, "%2.2f\t", candidates[i].s_off);

			for(j = 0; j < p->n_multi; j++) {
				bwt_multi1_t *q = p->multi + j;	
				if(q->mm <= opt->max_diff)
					candidates[i].n_in_gene[q->mm]++;
				else
					fprintf(stderr, "%d mm is found\n", q->mm);
			}
			//	36
			for(j = 0; j <= opt->max_diff; j++) {
				if(j > 0) 
					fprintf(result, " / ");
				fprintf(result, "%d-mm: %d", j, candidates[i].n_in_gene[j]);
			}
			fprintf(result, "\t");

			//	37
			for(j = 0; j < p->n_multi; j++) {
				bwt_multi1_t *q = p->multi + j;	
				bns_cnt_ambi(bns, q->pos, p->len, &seqid);
				fprintf(result, "%s:%d(%dmm) ", bns->anns[seqid].name, (int)(q->pos - bns->anns[seqid].offset + 1), q->mm);
			}

			fprintf(result, "\n");
		}
	}
	
	//if(verbose)
	//printf("   %s%011d is processed.\n", gen->prefix_gene, target_gene->ensembl_id);

	// free read structs
	for(i = 0; i < n_cand; i++) {
		free_cnt(candidates[i].n_in_gene);
		free_cnt(candidates[i].rna_seq);
		free_cnt(candidates[i].rna_seq_padded);
		p = &seqs[i];

		if(p->aln)	free_cnt(p->aln);
		if(p->md)	free_cnt(p->md); 
		if(p->multi)	free_cnt(p->multi);
		if(p->cigar)	free_cnt(p->cigar);
		free_cnt(p->name);
		free_cnt(p->seq);
		free_cnt(p->rseq);
		free_cnt(p->qual);
	}

	free_cnt(candidates);
	free_cnt(seqs);
	free_cnt(original);
	free_cnt(gen_seq);
	if(exon_seq != NULL)
		free_cnt(exon_seq);
	free_cnt(ese_seq);
	free_cnt(rna_only);
	free_cnt(rna_padded);

	return;
}


// Wrapper of write_sgrna_gene()
// Run write_sgrna_gene() for hash buckets from start to end
void run_write_library_batch (genome* gen, kmers* ese, char* output_path, double cut_offs, int cut_5, int cut_3, int cut_adjust, int start, int end, char* log_path, int verbose) {
	int hash_i, hash_start, hash_end;
	lentry *head_cursor;
	gene *target_gene = NULL;
	FILE *result, *log = NULL;

	// variables for alignment
	bwt_t *bwt;
	bntseq_t *bns;
	gap_opt_t *opt;
	opt = gap_init_opt();
	opt->max_diff = 3;
	opt->seed_len = SGRNA_LEN;
	opt->max_seed_diff = 3;
	opt->max_entries = 200000000; // 100 times to default
	opt->max_gapo = 0; 
	opt->max_gape = 0;
	opt->fnr = 0.0;

	printf("cut_5 / cut_3 = %d / %d\n", cut_5, cut_3);

	result = fopen(output_path, "w");
	if(result == NULL) {
		fprintf(stderr, "Cannot open %s\n", output_path);
		return;
	}

	if(log_path != NULL) {
		log = fopen(log_path, "w");
		if(log == NULL) {
			fprintf(stderr, "Cannot open log file.\n");
			return;
		}
	}

	printf("Load index..\n");
	bwt = load_bwt(gen->fasta_path);
	bns = bns_restore(gen->fasta_path);
	srand48(bns->seed);

	printf("Begin calculation..\n");
	// Write header
	// Column  0: Gene.ID: 		Ensembl Gene ID (ENSG000..)
	// Column  1: Symbol: 		Gene symbol
	// Column  2: Chr: 		Chromosome name
	// Column  3: G.start:	 	Gene start
	// Column  4: G.end: 		Gene end
	// Column  5: Num.trscs		Number of transcripts in the gene
	// Column  6: LCT.ID		Ensembl Transcript ID of the longest coding transcript (LCT) of the gene
	// Column  7: T.exon.ID		Ensembl Exon ID of the target exon (in the LCT)
	// Column  8: Sym		Asym (Asymmetric) or Sym (Symmetric)
	// Column  9: Trsc.info		Information of target site in all transcripts i.e) 1/1[3/3](ENSE00002326974-ENST00000392830:385aa) 1/1[3/3](ENSE00002326974-ENST00000240050:385aa)
	// Column 10: Num.CT.with.T	Number of coding transcripts containing target site
	// Column 11: Num.total.CT	Total number of coding transcripts of the gene
	// Column 12: Num.NCT.with.T	Number of noncoding transcripts containing target site
	// Column 13: Num.total.NCT	Total number of noncoding transcripts of the gene
	// Column 14: gRNA.seq		Guid RNA sequence (20bp)
	// Column 15: gRNA.PAM.seq	Guid RNA sequence with target PAM sequence (23bp)
	// Column 16: gRNA.PAM.pad.seq	Guid RNA sequence with target PAM sequence and padding sequence. (23 + PADDED_5 + padded_3 bp. In default 23 + 4 + 3 = 30bp)
	// Column 18: Sense		Sense or Antisense
	// Column 19: Pos.PAM		Position of PAM as genomic coordiate
	// Column 20: Pos.CUT.CDS	Position of PAM on CDS of the LCT (AA level). The CDS does not include UTRs
	// 				Adjusted by cut_adjust parameter. 
	// Column 21: LTC.CDS.len	CDS length of the LCT (AA sequence length)
	// Column 22: Pos.PAM.CDS.Ratio	Pos.CUT.CDS / LTC.CDS.Len
	// Column 23: No.T.exon.LCT	Ordinal number of the target exon among coding exons of the LCT
	// Column 24: Num.CE.in.LCT	Number of Coding exons of the LCT
	// Column 25: Exon.Loc		The target exon is the first / middle / last coding exon
	// Column 26: Pos.LEE		Position of the last coding Exon-Exon junction on cDNA (of the LCT)
	// Column 27: Pos.PTC.fs1	Position of PTC if there is 1bp frame shifting on cDNA
	// Column 28: Pos.PTC.fs2	Position of PTC if there is 2bp frame shifting on cDNA
	// Column 29: Pos.PTC.skip	Position of PTC if the target exon is skipped on cDNA
	// Column 30: Pos.LEE.skip	Position of the last coding exon-exon junction on exon-skipped cDNA
	// Column 31: Dist.PTC.LEE.skip	Distance from PTC from LEE on exon-skipped cDNA
	// Column 32: Dist.ESE.5	Distance to closest ESE at upstream
	// Column 33: Dist.ESE.3	Distance to closest ESE at downstream
	// Column 34: Dist.ESE.min	Distance to closest ESE at both side (min of 31 and 32)
	// Column 35: OffT.score	Off target score
	// Column 36: OffT.list		List of off targets
	// Column 37: OffT.detail.list	A detailed list of off targets
	fprintf(result, "Gene.ID\tSymbol\tChr\tG.start\tG.end\t");					// 0 1 2 3 4
	fprintf(result, "Num.trscs\tLCT.ID\tT.exon.ID\tSym\tTrsc.info\t");				// 5 6 7 8 9
	fprintf(result, "Num.CT.with.TE\tNum.total.CT\tNum.NCT.with.TE\tNum.total.NCT\tgRNA.seq\t");	// 10 11 12 13 14
	fprintf(result, "gRNA.PAM.seq\tgRNA.PAM.pad.seq\tSense\tPos.PAM\t");		// 15 16 18 19
	fprintf(result, "Pos.CUT.CDS(PAM%+d)\tLTC.CDS.len\tPos.CUT.CDS.Ratio\tNo.T.exon.LCT\tNum.CE.in.LCT\t", cut_adjust);	// 20 21 22 23 24
	fprintf(result, "Exon.Loc\tPos.LEE\tPos.PTC.fs1\tPos.PTC.fs2\tPos.PTC.skip\t");			// 25 26 27 28 29
	fprintf(result, "Pos.LEE.skip\tDist.PTC.LEE.skip\tDist.ESE.5\tDist.ESE.3\tDist.ESE.min\t");	// 30 31 32 33 34
	fprintf(result, "OffT.score\tOffT.list\tOffT.detail.list\n");					// 35 36 37

	hash_start = (start >= 0) ? start : 0;
	hash_end = (end >= 0) ? end : HASH_SIZE;
	hash_end = (hash_end <= HASH_SIZE) ? hash_end : HASH_SIZE;

	for(hash_i = hash_start; hash_i < hash_end; hash_i++) {
		head_cursor = gen->gene_head[hash_i];

		printf(" %d / %d hash bucket is processing.\n", hash_i, HASH_SIZE);

		FOR_LENTRY(head_cursor, gen->gene_head[hash_i]) {
			target_gene = (gene*)head_cursor->entry;
			target_gene->n_functional_trsc = 0;

			if(verbose)
				printf("  Gene: %s%011d (%s)\n", gen->prefix_gene, target_gene->ensembl_id, target_gene->symbol); 

			write_sgrna_gene(gen, ese, bwt, bns, opt, target_gene, result, log, cut_offs, cut_5, cut_3, cut_adjust, verbose);
		}
	}

	printf("Processing for %d - %d hash buckets were finished.\n", hash_start, hash_end - 1);

	fclose(result);
	if(log != NULL)
		fclose(log);

	bwt_destroy(bwt); // free BWT
	bns_destroy(bns);
}

/**** Function: run_interactive_offtarget *******************************************************************************
 *															*
 * Description: This function is for testing off-target searching function.						*
 * 		It accepts a sequence from STDIN, then search its off-target up to 3 mismatches, then 
 * 		report off-targets to STDOUT
 *															*
 ***********************************************************************************************************************/
int run_interactive_offtarget(int argc, char* argv[]) {
	char buf[BUF_LEN], *fasta = NULL, *stop;
	int c, i, j, maxdiff = -1;
	buf[0] = '\0';
	genome* gen;
	crispr_candidate* candidate;
	
	while ((c = getopt(argc, argv, "f:n:")) >= 0) {
		switch (c) {
		case 'f': 
			fasta = optarg;	
			break;
		case 'n': 
			maxdiff = strtol(optarg, &stop, 10);	
			break;
		case '?':
			printf("Unrecognized option: %c\n", c);
		}
	}

	if (fasta == NULL || maxdiff < 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:  CRISPinatoR interactive [options]\n\n");
		fprintf(stderr, "Options: -f     fasta file (bwa indexed)\n");
		fprintf(stderr, "Options: -n     maximum number of mismatch\n");
		fprintf(stderr, "\n");
		return 1;
	}

	printf("Reference: %s\n", fasta);
	printf("Max. no. of diff = %d\n", maxdiff);

	gen = load_genome(fasta);
	if(gen == NULL) 
		return 1;

	// variables for alignment
	bwa_seq_t *seqs, *p;
	int seqid;
	bwt_t *bwt;
	bntseq_t *bns;
	gap_opt_t *opt;
	opt = gap_init_opt();
	opt->max_diff = maxdiff;
	opt->max_seed_diff = maxdiff;
	opt->max_entries = 200000000; // 100 times larger than default
	opt->max_gapo = 0; 
	opt->max_gape = 0;
	opt->fnr = 0.0;

	seqs = (bwa_seq_t*)malloc_cnt(sizeof(bwa_seq_t));
	p = seqs;
	p->name = (char*)malloc_cnt(10 * sizeof(char));
	strncpy(p->name, "@test_seq", 9);
	p->name[9] = '\0';

	candidate = (crispr_candidate*)malloc_cnt(sizeof(crispr_candidate));
	candidate->target_gene = NULL;
	candidate->n_in_gene = (int*)malloc_cnt(sizeof(int) * (opt->max_diff + 1));
	candidate->n_out_gene = (int*)malloc_cnt(sizeof(int) * (opt->max_diff + 1));
	for(j = 0; j <= opt->max_diff; j++) {
		candidate->n_in_gene[j] = 0;
		candidate->n_out_gene[j] = 0;
	}
	candidate->pos_in_trsc = NULL;
	candidate->n_ese = 0;
	candidate->n_aff_trsc = 0;
	candidate->n_aff_functional_trsc = 0;
	candidate->s_pos = 0.0;
	candidate->s_off = 0.0;
	candidate->s_trsc = 0.0;
	candidate->s_ese = 0.0;
	candidate->s_sum = 0.0;
	candidate->aln = p;

	// load BWT
	printf("Load bwt index\n");
	bwt = load_bwt(gen->fasta_path);
	bns = bns_restore(gen->fasta_path);
	srand48(bns->seed);
	printf("Test off-target search function\n");	

	while(1) {
		printf("*************************************************\n");
		printf("Type sequence to be aligned (q to exit): ");
		fgets(buf, BUF_LEN - 1, stdin);

		if((buf[0] == 'q') && ((buf[1] == '\n') || (buf[1] == '0')))
			return 0;

		if(strlen(buf) <= 5) {
			printf("Please type more than 10bp\n");
			continue;
		}

		while(buf[strlen(buf) - 1] == '\n')
			buf[strlen(buf) - 1] = '\0';

		make_str_hc(buf);

		init_seq(p, strlen(buf));
		for(i = 0; i < strlen(buf); i++) {
			p->seq[i] = nst_nt4_table[(int)buf[i]];
			p->qual[i] = 'h';
		}

		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, 0);

		for(j = 0; j <= opt->max_diff; j++)
			candidate->n_in_gene[j] = 0;

		printf("Alignment..");
		do_alignment(gen, bwt, bns, candidate, seqs, 1, opt);
		printf(" finished.\n");

		for(i = 0; i < p->n_multi; i++) {
			bwt_multi1_t *q = p->multi + i;	
			if(q->mm == 0) {
				for(j = i + 1; j < p->n_multi; ++j)
					p->multi[j - 1] = p->multi[j]; 
				p->n_multi--;
				break;
			}
		}

		candidate->s_off = calculate_fs_offtarget_score(gen, candidate, bns, NULL);

		for(j = 0; j < p->n_multi; j++) {
			bwt_multi1_t *q = p->multi + j;	
			candidate->n_in_gene[q->mm]++;	
		}

		printf("     Off-target result (remove one 0-mm match, probably on-target\n");
		printf("       # of mismatch   # of off-targets\n");
		printf("       ------------------------------------\n");
		for(j = 0; j <= opt->max_diff; j++) {
		printf("       %d mismatch     %3d\n", j, candidate->n_in_gene[j]);
		}
		printf("       ------------------------------------\n\n");
		printf("       Off-target score = %f\n", candidate->s_off);

		if (p->n_multi) {
			//int k = p->n_multi < 100 ? p->n_multi : 100;
			int k = p->n_multi;
			printf("       > List of off-targets (at most 100 off-targets are printed)\n");
			printf("         Chr\tpos (strand)\t# of mismatch\n");
			printf("         ------------------------------------\n");
			for (j = 0; j < k; ++j) {
				bwt_multi1_t *q = p->multi + j;
				bns_cnt_ambi(bns, q->pos, p->len, &seqid);
				printf("         %s\t%d (%c)\t", bns->anns[seqid].name,
				        (int)(q->pos - bns->anns[seqid].offset + 1), q->strand? '-' : '+');
				printf("%dmm\n", q->mm);
			}
			printf("\n");
		}

		if(p->aln)	free_cnt(p->aln);
		if(p->md)	free_cnt(p->md); 
		if(p->multi)	free_cnt(p->multi);
		if(p->cigar)	free_cnt(p->cigar);
		free_cnt(p->seq);
		free_cnt(p->rseq);
		free_cnt(p->qual);
	}

	free_cnt(candidate->n_in_gene);
	free_cnt(candidate->n_out_gene);
	free_cnt(candidate);

	free_cnt(p->name);
	free_cnt(p);
}


/**** Function: crispinator_write_library *******************************************************************************
 *															*
 * Description: This function writes generates gRNA for input genes.							*
 * 		In default, it generates all possible gRNAs for genes in the given Ensembl annotation.		 	*
 * 		Unfortunatly multithread is not supported yet.								*
 *              The genes (from Ensembl annotation) are hashed using 2048 buckets in default), 		                *
 *              so multiple instances may be run in parellel using -f (from) and -t (to) options. 			*
 *              (-f in inclusive, -t is exclusive.)									*
 *              For example, -f 0 -t 64 will process bucket 0 to bucket 63.						*
 *															*
 ***********************************************************************************************************************/
int crispinator_write_library(int argc, char* argv[]) {
	int c, start = -1, end = -1, cut_adjust = -3, verbose = 0;
	unsigned int cut5 = 20, cut3 = 10;
	char *fasta = NULL, *ensembl = NULL, *ese_path = NULL, *output = NULL, *stop, *log_path = NULL;
	double cut_offscore = 80.0;
	kmers *ese;
	genome* gen;
	
	while ((c = getopt(argc, argv, "r:d:e:s:o:s:5:3:c:l:f:t:iv")) >= 0) {
		switch (c) {
		case 'r': 
			fasta = optarg;	
			break;
		case 'd': 
			ensembl = optarg; 
			break;
		case 'e':
			ese_path = optarg;
			break;
		case 's':
			cut_offscore = strtod(optarg, &stop);
			break;
		case 'o':
			output = optarg;
			break;
		case '5':
			cut5 = strtol(optarg, &stop, 10);
			break;
		case '3':
			cut3 = strtol(optarg, &stop, 10);
			break;
		case 'c':
			cut_adjust = strtol(optarg, &stop, 10);
			break;
		case 'l':
			log_path = optarg;
			break;
		case 'f':
			start = strtol(optarg, &stop, 10);
			break;
		case 't':
			end = strtol(optarg, &stop, 10);
			break;
		case 'v':
			verbose = 1;
			break;
		case '?':
			printf("Unrecognized option\n");
		}
	}

	if ((fasta == NULL) || (ensembl == NULL) || (ese_path == NULL) || ((output == NULL))) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:  CRISPinatoR library [options]\n\n");
		fprintf(stderr, "Options: -r     fasta file (indexed)\n");
		fprintf(stderr, "         -d     Ensembl database path\n");
		fprintf(stderr, "         -e     ESE list\n");
		fprintf(stderr, "         -o     Outut path\n");
		fprintf(stderr, "         -s     Off-target score cutoff (default: 80.0)\n");
		fprintf(stderr, "         -c     Adjust nt for cut site from PAM. Minus value means upstream, while\n");
		fprintf(stderr, "                plus value means downstream. (default: -3. Max, min: +20, -20)\n");
		fprintf(stderr, "         -5     Cutoff for distance to upstream ESE from the cut site (default: 20)\n");
		fprintf(stderr, "         -3     Cutoff for distance to downstram ESE from the cut site (default: 10)\n");
		fprintf(stderr, "         -l     Path to log file (if it is not specified, no log will be written.)\n");
		fprintf(stderr, "         -f     Index where calculation begin from (default: first of hash)\n");
		fprintf(stderr, "         -t     Index where calculation stop (default: end of hash)\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if((cut_adjust < -20) && (cut_adjust > 20)) {
		fprintf(stderr, "Too large cut_adjust: %d\n", cut_adjust);
		fprintf(stderr, "Please use -20 ~ +20\n");
		return 1;
	}

	gen = load_genome(fasta);
	if(gen == NULL) 
		return 1;
	
	if(load_ensembl_data(gen, ensembl))
		return 1;

	ese = load_kmers(ese_path);

	if(ese == NULL)
		return 1;
	
	run_write_library_batch (gen, ese, output, cut_offscore, cut5, cut3, cut_adjust, start, end, log_path, verbose);

	free_kmers(ese);
	free_genome(gen);
		
	return 0;
}


