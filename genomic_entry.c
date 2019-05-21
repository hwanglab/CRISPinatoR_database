/*
 * genomic_entry.c 	in CRISPinatoR
 *
 *  Created on: 2016. 6. 24.
 *      Author: Yunku Yeu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "genomic_entry.h"
#include "utils.h"
#include "crispinator.h"


// open FASTA file, and load .fai information
genome* load_genome(char* path) {
	int i, j;
	char buf[BUF_LEN];
	FILE* fai = NULL;
	genome* new_genome;

#ifdef DEBUG
	printf("Before load_genome() ==> ");
	print_cnt();
#endif
	printf("Loading fasta file..\n");
	new_genome = (genome*)malloc_cnt(sizeof(struct genome));

	if(new_genome == NULL) {
		fprintf(stderr, "Not enough memory!");
		exit(1);
	}

	new_genome->fasta = fopen(path, "r");
	if(new_genome->fasta == NULL) {
		free_cnt(new_genome);
		fprintf(stderr, "Cannot open fasta file.\n");
		return NULL;
	}
	new_genome->fasta_path = (char*)malloc_cnt(sizeof(char) * 512);
	strncpy(new_genome->fasta_path, path, strlen(path));
	
	strncpy(buf, path, strlen(path));
	buf[strlen(path)] = '\0';
	strcat(buf, ".fai");
	fai = fopen(buf, "r");
	if(fai == NULL) {
		free_cnt(new_genome);
		fprintf(stderr, "Cannot open fai file. Please index the fasta file using samtools faidx.\n");
		return NULL;
	}

	new_genome->num_chrs = 0;
	while(fgets(buf, BUF_LEN, fai) != NULL) {
		new_genome->num_chrs++;
	}
	new_genome->chrs = (chr**)malloc_cnt(sizeof(chr*) * new_genome->num_chrs);
	if(new_genome->chrs == NULL) {
		fprintf(stderr, "Not enough memory!");
		free_cnt(new_genome);
		exit(1);
	}
	new_genome->exon_head = new_genome->trsc_head = new_genome->gene_head = NULL;

	fseek(fai, 0, SEEK_SET);
	i = 0;
	while(fgets(buf, BUF_LEN, fai) != NULL) {
		new_genome->chrs[i]       = (chr*)malloc_cnt(sizeof(struct chromosome));
		new_genome->chrs[i]->name = (char*) malloc_cnt(sizeof(char) * MAX_NAME_LEN);
		new_genome->chrs[i]->num_genes = 0;

		if((new_genome->chrs[i] == NULL) || (new_genome->chrs[i]->name == NULL)) {
			printf("Not enough memory!");
			for(j = 0; j < i; j++) {
				free_cnt(new_genome->chrs[j]->name);
				free_cnt(new_genome->chrs[j]);
			}
			if(new_genome->chrs[i]->name != NULL)
				free_cnt(new_genome->chrs[i]->name);
			if(new_genome->chrs[i] != NULL)
				free_cnt(new_genome->chrs[i]);

			free_cnt(new_genome->chrs);
			free_cnt(new_genome);
			exit(1);
		}

		new_genome->chrs[i]->idx = i;
		new_genome->chrs[i]->genes = NULL;
		sscanf(buf, "%s\t%d\t%d\t%hd\t%hd\n", 	new_genome->chrs[i]->name,
				&(new_genome->chrs[i]->length),    &(new_genome->chrs[i]->offset),
				&(new_genome->chrs[i]->linebases), &(new_genome->chrs[i]->linelength));
		i++;
	}

	new_genome->prefix_gene = new_genome->prefix_trsc = new_genome->prefix_exon = new_genome->prefix_protein = NULL;
	new_genome->plen_gene = new_genome->plen_trsc = new_genome->plen_exon = new_genome->plen_protein = 0;

#ifdef DEBUG
	printf("After load_genome() ==> ");
	print_cnt();
#endif

	return new_genome;
}

// this function convert genomic coordinate to transctipt-based coordinate. (Include "UTRs")
// It returns 1-base coordinate.
unsigned int convert_genomic_coord_to_trsc_coord(transcript* trsc, unsigned int genomic_coord) {
	lentry* cursor;
	exon* cur_exon;
	unsigned int trsc_coord = 1;

	// The original coordinate should be inside of the transcript
	if((genomic_coord < trsc->tr_start) || (trsc->tr_end < genomic_coord))
		return 0;

	cursor = trsc->first_exon;
	while(cursor != NULL) {
		cur_exon = (exon*)cursor->entry;
		
		if((cur_exon->start <= genomic_coord) && (genomic_coord <= cur_exon->end)) {
			if(trsc->strand == FORWARD)
				return trsc_coord + (genomic_coord - cur_exon->start);	// return 1-base coord.
			else	
				return trsc_coord + (cur_exon->end - genomic_coord);	// return 1-base coord.
		}
		
		trsc_coord += (cur_exon->end - cur_exon->start + 1);
		cursor = cursor->next;
	}

	// This function shouldn't come here
	return 0;
}

// This function concatnates sequences of all exons in a transcript
// Parameters ***************************************************************
//   genome* gen: 		genome struct
//   transcript* trsc:		transcript struct
//   int add_buffer:		Sometimes short buffer sequence is required due to phasing 
// **************************************************************************
// Return: DNA sequence of a transcript.
// **************************************************************************
char* concat_exon_seqs(genome* gen, transcript* trsc, int add_buffer) {
	int i, cnt = 0, len;
	char** exon_seqs, *merged_seq;
	lentry* cursor = trsc->first_exon;
	exon* cur_exon;

	// Count number of exons (To be checked: it is stored in trsc->num_exons??)
	while((cursor != NULL) && (cursor->entry != NULL)) {
		cnt++;
		cursor = cursor->next;
	}

	exon_seqs = (char**)malloc_cnt(sizeof(char*) * cnt);

	i = 0;	
	cursor = trsc->first_exon;
	len = 0;
	while((cursor != NULL) && (cursor->entry != NULL)) {
		cur_exon = (exon*)cursor->entry;

		// If add_buffer is set with ADD_3BP_BUFFER, adds +3 bps as buffer at the 5' end of the first exon,
		// (depending on the transcript's strand)
		if(i == 0) {
			if(trsc->strand == FORWARD)
				exon_seqs[i] = get_subseq_by_idx(gen, 
						cur_exon->chr->idx, 
						(add_buffer == ADD_3BP_BUFFER) ? cur_exon->start - 3 : cur_exon->start, 
						cur_exon->end);
			else 
				exon_seqs[i] = get_subseq_by_idx(gen, 
						cur_exon->chr->idx, 
						cur_exon->start, 
						(add_buffer == ADD_3BP_BUFFER) ? cur_exon->end + 3 : cur_exon->end);
		}
		else 
			exon_seqs[i] = get_subseq_by_idx(gen, cur_exon->chr->idx, cur_exon->start, cur_exon->end);

		// swite upper / lower case
		if(i % 2 == 0)
			make_str_hc(exon_seqs[i]);
		else
			make_str_lc(exon_seqs[i]);
		cursor = cursor->next;
		len += strlen(exon_seqs[i]);
		i++;
	}

	// Concatenate all sequences
	merged_seq = (char*)malloc_cnt(sizeof(char) * (len + 1));
	merged_seq[0] = '\0';
	
	if(trsc->strand == FORWARD) {	
		for(i = 0; i < cnt; i++) {
			strcat(merged_seq, exon_seqs[i]);
		}
	}
	else {
		// Put them in reverse order.
		// Please note that, when the transcript is REVERSE strand,
		// the last exon has the smallest genomic coordinate.
		for(i = cnt - 1; i >= 0; i--) {
			strcat(merged_seq, exon_seqs[i]);
		}
		// Convert the concatnated string into Reverse Complement.
		merged_seq = create_rc_free_exist_len(merged_seq, len + 1);
	}

	for(i = 0; i < cnt; i++)
		free_cnt(exon_seqs[i]);
	free_cnt(exon_seqs);

	return merged_seq;
}


// This function read FASTA file and return specified subsequence
// Parameters ***************************************************************
//   genome* gen: 		genome struct
//   int chr_idx:		Chromosome index. get_subseq_by_name() can be used as an alternative if chr name is available.
//   unsigned int start: 	Subsequence range start
//   unsigned int end:		Subsequence range end
// **************************************************************************
// Return: Subsequence from specified region
// **************************************************************************
char* get_subseq_by_idx(genome* gen, int chr_idx, unsigned int start, unsigned int end) {
	char buf[BUF_LEN], *subseq = NULL;
	unsigned int offset, size_new, length_to_read = end - start + 1;
	int pos = 0, read_len;

	if((start > end) || (start < 1) || (end < 1)) {
		fprintf(stderr, "[ERROR] Invalid genomic range (%u-%u)\n", start, end);
		return NULL;
	}
	if(end > gen->chrs[chr_idx]->length) {
		fprintf(stderr, "[ERROR] Genomic range is out of chromosome %s\n", gen->chrs[chr_idx]->name);
		return NULL;
	}
	
	// change 1-base coordinates to 0-base;
	start--;
	end--;
	
	// alloc subseq.
	if(chr_idx < gen->num_chrs) {
		subseq = (char*)malloc_cnt(sizeof(char) * (length_to_read + gen->chrs[chr_idx]->linelength));
		buf[0] = subseq[0] = '\0';

		size_new = gen->chrs[chr_idx]->linelength - gen->chrs[chr_idx]->linebases;
		offset = gen->chrs[chr_idx]->offset + start + (start / gen->chrs[chr_idx]->linebases * size_new);
		//start = start % gen->chrs[chr_idx]->linebases;

		// subsequence is read from FASTA file directly.
		// Move file pointer
		fseek(gen->fasta, offset, SEEK_SET);

		while(fgets(buf, BUF_LEN, gen->fasta) != NULL) {
			read_len = strlen(buf);
			if(buf[read_len - 1] == '\n')
				buf[--read_len] = '\0';
		
			strncpy(subseq + pos, buf, read_len);
			pos += read_len;

			if (length_to_read < pos)
				break;
		}

		subseq[length_to_read] = '\0';
	}
	else {
		fprintf(stderr, "Chromosome doesn't exist\n");
	}

	return subseq;
}

// Wrapper function of get_subseq_by_idx
// Parameters ***************************************************************
//   genome* gen: 		genome struct
//   char* chr_name:		Chromosome name (as written in the .fai file)
//   unsigned int start: 	Subsequence range start
//   unsigned int end:		Subsequence range end
// **************************************************************************
// Return: Subsequence from specified region
// **************************************************************************
char* get_subseq_by_name(genome* gen, char* chr_name, unsigned int start, unsigned int end) {
	unsigned int i;

	for(i = 0; i < gen->num_chrs; i++) {
		if(strcmp(gen->chrs[i]->name, chr_name) == 0) {
			break;
		}
	}
	
	return get_subseq_by_idx(gen, i, start, end);
}


// Utility to get a chromosome struct
// Parameters ***************************************************************
//   genome* gen: 		genome struct
//   int idx:			Chromosome index
// **************************************************************************
// Return: chromosome struct
// **************************************************************************
chr* find_chr_idx(genome* gen, int idx) {
	if((idx >= 0) && (idx < gen->num_chrs))
		return gen->chrs[idx];
	else
		return NULL;
}

// Utility to get a chromosome struct (alternative of find_chr_idx())
// Parameters ***************************************************************
//   genome* gen: 		genome struct
//   char* chr_name:		Chromosome name (as written in the .fai file)
// **************************************************************************
// Return: chromosome struct
// **************************************************************************
chr* find_chr_name(genome* gen, char* name) {
	chr* chr = NULL;
	int i = 0;

	while(i < gen->num_chrs) {
		chr = gen->chrs[i];
		if(strncmp(chr->name, name, strlen(name)) == 0)
			return chr;
		i++;
	}

	return NULL;
}


// Lookup a gene using gene symbol, when Ensembl ID is unavilable.
gene* find_gene_by_name(genome* gen, char* name) {
	lentry* cursor;
	gene* target = NULL;
	int i;

	// make gene name to uppercases
	make_str_hc(name);
	//printf("finding gene in %s\n", name); 
	
	for(i = 0; i < HASH_SIZE; i++) {
		cursor = gen->gene_head[i];
		while((cursor != NULL) && (cursor->entry != NULL)) {
			target = (gene*)cursor->entry;
			//printf("ensembl ID:%d, symbol=%s\n", target->ensembl_id, target->symbol);
			//getchar();
			if(strcmp(target->symbol, name) == 0) 
				return target;	
			cursor = cursor->next;
		}
	}

	return NULL;
}

// Find gene or trsc or exon using their Ensembl ID
// head can be any lentry
void* find_entry(lentry* head, unsigned int id, int type) {
	lentry* cur = head;
	if(type == GENE) {
		while((cur != NULL) && (cur->entry != NULL)) {
			if(((gene*)cur->entry)->ensembl_id == id)
				return cur->entry;
			cur = cur->next;
		}
	}
	else if(type == TRSC) {
		while((cur != NULL) && (cur->entry != NULL)) {
			if(((transcript*)cur->entry)->ensembl_id == id)
				return cur->entry;
			cur = cur->next;
		}
	}
	else if(type == EXON) {
		while((cur != NULL) && (cur->entry != NULL)) {
			if(((exon*)cur->entry)->ensembl_id == id)
				return cur->entry;
			cur = cur->next;
		}
	}

	return NULL;
}

// Wrappers for gene, transcript, and exon respectively
gene* find_gene_by_id(genome* gen, unsigned int ensembl_id) {
	return (gene*)find_entry(gen->gene_head[ensembl_id % HASH_SIZE], ensembl_id, GENE);
}

transcript* find_trsc_by_id(genome* gen, unsigned int ensembl_id) {
	return (transcript*)find_entry(gen->trsc_head[ensembl_id % HASH_SIZE], ensembl_id, TRSC);
}

transcript* find_trsc_by_id_in_gene(gene* gene, unsigned int ensembl_id) {
	return (transcript*)find_entry(gene->first_trsc, ensembl_id, TRSC);
}

exon* find_exon_by_id(genome* gen, unsigned int ensembl_id) {
	return (exon*)find_entry(gen->exon_head[ensembl_id % HASH_SIZE], ensembl_id, EXON);
}

// create an lentry instance and add it next to the head
void create_and_add_entry_last(lentry* head, void* entry) {
	lentry* new_entry = (lentry*)malloc_cnt(sizeof(lentry));
	new_entry->entry = entry;
	if(head == NULL)
		head = new_entry;
	else {
		new_entry->next = head->next;
		head->next = new_entry;
	}
}


// free lentry, but not its entry.
void free_lentry_only(lentry* head) {
	lentry* cur;

	while(head != NULL) {
		cur = head;
		head = head->next;
		free_cnt(cur);
	}
}



// find entry in the list head, if it doesn't exist in the list, add it
void find_entry_add(lentry* head, void* entry) {
	int found = 0;
	lentry* cur = head, *new_entry;

	if(head->entry == NULL) {
		head->entry = entry;
		return;
	}

	while((cur->next != NULL) && (!found)) {
		if(cur->entry == entry)
			found = 1;
		else
			cur = cur->next;
	}

	// if entry exists in the list, do nothing.
	if(!found) {
		// last entry was not checked in previous while loop
		if(cur->entry != entry) {
			new_entry = (lentry*)malloc_cnt(sizeof(lentry));
			new_entry->entry = entry;
			new_entry->next = NULL;
			cur->next = new_entry;
		}
	}
}


// Load Ensembl table. 
// Column names are hard coded.
int load_ensembl_data(genome* gen, char* path) {
	FILE* ensembl_file;
	char* buf, *token, *stop, *tofree, *ptr;
	char temp[MAX_ATTR_LEN];
	char** buf_token;
	int*   idx_of;
	gene* cur_gene = NULL, *prev_gene = NULL;
	transcript* cur_trsc = NULL, *prev_trsc = NULL;
	exon* cur_exon = NULL;
	lentry *cur, *cur_coding, *new_entry, *outer;

	int header_size, num_of_features = 21, cnt = 0, uniq_exon = 0, uniq_trsc = 0, uniq_gene = 0;//, *hash_cnt;
	unsigned int line_gene_id, line_trsc_id, line_exon_id, line_protein_id;
	int gene_id = 0, trsc_id = 1, exon_id = 2, chr_name = 3, strand = 4,
		gene_symbol = 5, exon_start = 6, exon_end = 7, exon_rank = 8,
		gene_start = 9, gene_end = 10, trsc_start = 11, trsc_end = 12,
		cdna_start = 13, cdna_end = 14, start_phase = 15, end_phase = 16, 
		cds_start = 17, cds_end = 18,  protein_id = 19, tr_type = 20,
		i, res = 0, is_coding, cnt_trsc_without_CDS;

#ifdef DEBUG
	printf("Before load_ensembl_data() ==> ");
	print_cnt();
#endif
	
	printf("Loading ensembl data..\n");
	ensembl_file = fopen(path, "r");
	if(ensembl_file == NULL) {
		fprintf(stderr, "Cannot open ensembl file (%s)\n", path);
		return 1;
	}

	buf = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
	tofree = buf;

	header_size = 0;
	fgets(buf, BUF_LEN, ensembl_file);
	while((token = strsep(&buf, "\t")) != NULL)
		header_size++;

	fseek(ensembl_file, 0, SEEK_SET);
	printf("Number of columns = %d\n", header_size);	
	
	buf = tofree;	
	idx_of = (int*)malloc_cnt(sizeof(int) * num_of_features);
	if(idx_of == NULL)
		return 2;

	buf_token = (char**)malloc_cnt(sizeof(char*) * header_size);
	if((buf == NULL) || (buf_token == NULL)) {
		return 2;
	}

	for(i = 0; i < num_of_features; i++) {
		idx_of[i] = -1;
	}
	for(i = 0; i < header_size; i++) {
		buf_token[i] = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
		if(buf_token[i] == NULL) {
			return 2;
		}
	}

	gen->gene_head = (lentry**)malloc_cnt(sizeof(lentry) * HASH_SIZE);
	gen->trsc_head = (lentry**)malloc_cnt(sizeof(lentry) * HASH_SIZE);
	gen->exon_head = (lentry**)malloc_cnt(sizeof(lentry) * HASH_SIZE);

	for(i = 0; i < HASH_SIZE; i++) {
		gen->gene_head[i] = gen->trsc_head[i] = gen->exon_head[i] = NULL;
	}

	if((gen->gene_head == NULL) || (gen->trsc_head == NULL) || (gen->exon_head == NULL)) {
		res = 2;
		goto exit;
	}

	fgets(buf, BUF_LEN, ensembl_file);
	i = -1;
	if(buf[strlen(buf) - 1] == '\n')
		buf[strlen(buf) - 1] = '\0';

	while((i < header_size) && ((token = strsep(&buf, "\t")) != NULL)) {
		i++;
		if(strcmp(token, "Ensembl Gene ID") == 0)		idx_of[gene_id] = i;
		else if(strcmp(token, "Ensembl Transcript ID") == 0)	idx_of[trsc_id] = i;
		else if(strcmp(token, "Ensembl Protein ID") == 0)	idx_of[protein_id] = i;
		else if(strcmp(token, "Chromosome Name") == 0)		idx_of[chr_name] = i;
		else if (strcmp(token, "Strand") == 0)			idx_of[strand] = i;
		else if (strcmp(token, "Exon Chr End (bp)") == 0)	idx_of[exon_end] = i;
		else if (strcmp(token, "Exon Chr Start (bp)") == 0)	idx_of[exon_start] = i;
		else if (strcmp(token, "Exon Rank in Transcript") == 0)	idx_of[exon_rank] = i;
		else if (strcmp(token, "cDNA coding start") == 0)	idx_of[cdna_start] = i;
		else if (strcmp(token, "cDNA coding end") == 0)		idx_of[cdna_end] = i;
		else if (strcmp(token, "end phase") == 0)		idx_of[end_phase] = i;
		else if (strcmp(token, "start phase") == 0)		idx_of[start_phase] = i;
		else if (strcmp(token, "Ensembl Exon ID") == 0)		idx_of[exon_id] = i;
		else if (strcmp(token, "CDS Start") == 0)		idx_of[cds_start] = i;
		else if (strcmp(token, "CDS End") == 0)			idx_of[cds_end] = i;
		else if (strcmp(token, "Associated Gene Name") == 0)	idx_of[gene_symbol] = i;
		else if (strcmp(token, "Transcript Start (bp)") == 0)	idx_of[trsc_start] = i;
		else if (strcmp(token, "Transcript End (bp)") == 0)	idx_of[trsc_end] = i;
		else if (strcmp(token, "Gene Start (bp)") == 0)		idx_of[gene_start] = i;
		else if (strcmp(token, "Gene End (bp)") == 0)		idx_of[gene_end] = i;
		else if (strcmp(token, "Transcript type") == 0)		idx_of[tr_type] = i;
	} 

	buf = tofree;

	for(i = 0; i < num_of_features; i++) {
		if((i == end_phase) && (idx_of[end_phase] == -1)) {
			fprintf(stderr, "Does not contain the end phase info.\n");
			// do nothing. the end phase can be empty according to Ensembl version.
		}
		else if(idx_of[i] == -1) {
			fprintf(stderr, "Some header(%d th) are not found!\n", i);
			fprintf(stderr, "It fixed to 0\n");
			idx_of[i] = 0;
			//goto exit;
		}
	}

	// load ensembl data
	while(fgets(buf, BUF_LEN, ensembl_file)) {
		// remove last newline character
		if(buf[strlen(buf) - 1] == '\n')
			buf[strlen(buf) - 1] = '\0';

		i = 0;
		while((i < header_size) && ((token = strsep(&buf, "\t")) != NULL)) {
			strncpy(buf_token[i], token, strlen(token));
			ptr = buf_token[i];
			ptr[strlen(token)] = '\0';
			i++;
		}

		if(i < header_size) {
			int j;
			fprintf(stderr, "line %s has %d tokens\n", tofree, i);
			for(j = 0; j < i; j++) {
				fprintf(stderr, "token %d: %s\n", j, buf_token[j]);
			}
			getchar();
		}
	
		// Read Chromosome name
		buf = tofree;
		if(find_chr_name(gen, buf_token[idx_of[chr_name]]) == NULL) {
			continue;
		}

		// Read exon ID
		ptr = buf_token[idx_of[exon_id]];
		if(strlen(ptr) > 0) {
			if(gen->prefix_exon == NULL) {	// Prefix can be differed by species
				for(i = 0; i < strlen(ptr); i++) {
					if(('0' <= ptr[i]) && (ptr[i] <= '9'))
						break;
				}
				gen->prefix_exon = (char*)malloc(sizeof(char) * 32);
				strncpy(gen->prefix_exon, ptr, i);
				gen->plen_exon = i;
			}
			strncpy(temp, ptr + gen->plen_exon, strlen(ptr) - gen->plen_exon);	
			temp[strlen(ptr) - gen->plen_exon] = '\0';
			//printf("temp=%s (%d)\n", temp, strlen(temp));	
			line_exon_id = strtoul(temp, &stop, 10);
		}
		else {
			line_exon_id = -1;
		}
		
		// Read Transcript ID
		ptr = buf_token[idx_of[trsc_id]];
		if(strlen(ptr) > 0) {
			if(gen->prefix_trsc == NULL) {	
				for(i = 0; i < strlen(ptr); i++) {
					if(('0' <= ptr[i]) && (ptr[i] <= '9'))
						break;
				}
				gen->prefix_trsc = (char*)malloc(sizeof(char) * 32);
				strncpy(gen->prefix_trsc, ptr, i);
				gen->plen_trsc = i;
			}

			strncpy(temp, ptr + gen->plen_trsc, strlen(ptr) - gen->plen_trsc);
			temp[strlen(ptr) - gen->plen_trsc] = '\0';
			line_trsc_id = strtoul(temp, &stop, 10);
		}
		else {
			line_trsc_id = -1;
		}

		// Read Gene ID
		ptr = buf_token[idx_of[gene_id]];
		if(strlen(ptr) > 0) {
			if(gen->prefix_gene == NULL) {
				for(i = 0; i < strlen(ptr); i++) {
					if(('0' <= ptr[i]) && (ptr[i] <= '9'))
						break;
				}
				gen->prefix_gene = (char*)malloc(sizeof(char) * 32);
				strncpy(gen->prefix_gene, ptr, i);
				gen->plen_gene = i;
			}

			strncpy(temp, ptr + gen->plen_gene, strlen(ptr) - gen->plen_gene);
			temp[strlen(ptr) - gen->plen_gene] = '\0';
			line_gene_id = strtoul(temp, &stop, 10);
		}
		else {
			line_gene_id = -1;
		}

		// Read protein ID
		ptr = buf_token[idx_of[protein_id]];
		if(strlen(ptr) > 0) {
			if(gen->prefix_protein == NULL) {
				for(i = 0; i < strlen(ptr); i++) {
					if(('0' <= ptr[i]) && (ptr[i] <= '9'))
						break;
				}
				gen->prefix_protein = (char*)malloc(sizeof(char) * 32);
				strncpy(gen->prefix_protein, ptr, i);
				gen->plen_protein = i;
			}

			strncpy(temp, ptr + gen->plen_protein, strlen(ptr) - gen->plen_protein);
			temp[strlen(ptr) - gen->plen_protein] = '\0';
			line_protein_id = strtoul(temp, &stop, 10);
		}
		else {
			line_protein_id = 0;
		}

		if((line_gene_id < 0) || (line_trsc_id < 0) || (line_exon_id < 0)) {
			printf("ID NULL: %s\n", buf);
			continue;
		}
		//printf("gene ID = %s -> %u\n", buf_token[idx_of[gene_id]], line_gene_id);


		// check if current exons have been seen already
		cur_exon = find_exon_by_id(gen, line_exon_id);

		if(cur_exon == NULL) {	// assign new exon
			new_entry = (lentry*)malloc_cnt(sizeof(lentry));
			cur_exon = (exon*)malloc_cnt(sizeof(struct exon));

			if((cur_exon == NULL) || (new_entry == NULL)) {
				res = 2;
				goto exit;
			}

			// lentry for track transcripts containing this exon
			cur_exon->first_trsc = (lentry*)malloc_cnt(sizeof(lentry));
			cur_exon->first_trsc->next = cur_exon->first_trsc->entry = NULL;
			
			// Store information
			cur_exon->ensembl_id = line_exon_id;
			cur_exon->chr = find_chr_name(gen, buf_token[idx_of[chr_name]]);
			cur_exon->strand = strtol(buf_token[idx_of[strand]], &stop, 10);
			cur_exon->start = strtoul(buf_token[idx_of[exon_start]], &stop, 10);
			cur_exon->end = strtoul(buf_token[idx_of[exon_end]], &stop, 10);
			cur_exon->start_phase = (short)strtol(buf_token[idx_of[start_phase]], &stop, 10);

			// End phase will be updated later
			if(idx_of[end_phase] >= 0)
				cur_exon->end_phase = (short)strtol(buf_token[idx_of[end_phase]], &stop, 10);
			else
				cur_exon->end_phase = -2;

			// add to new entry in the front of list
			new_entry->entry = cur_exon;
			new_entry->next = gen->exon_head[line_exon_id % HASH_SIZE];	// Hashed for faster lookup
			gen->exon_head[line_exon_id % HASH_SIZE] = new_entry;
			uniq_exon++;
		}
		// Exon was created or found now
		
		
		if(cur_exon->start_phase != (short)strtol(buf_token[idx_of[start_phase]], &stop, 10)) {
			printf("%s%11u has various phase info.\n", gen->prefix_exon, cur_exon->ensembl_id);
		}
		if(idx_of[end_phase] >= 0) {
			if(cur_exon->end_phase != (short)strtol(buf_token[idx_of[end_phase]], &stop, 10)) {
				printf("%s%11u has various phase info.\n", gen->prefix_exon, cur_exon->ensembl_id);
			}
		}
		
		// Look up transcript based on its ID
		if((prev_trsc != NULL) && (prev_trsc->ensembl_id == line_trsc_id))
			cur_trsc = prev_trsc;
		else
			cur_trsc = find_trsc_by_id(gen, line_trsc_id);

		// Create a new transcript struct
		if(cur_trsc == NULL) {
			new_entry = (lentry*)malloc_cnt(sizeof(lentry));
			cur_trsc = (transcript*)malloc_cnt(sizeof(struct transcript));
			cur_trsc->first_exon = (lentry*)malloc_cnt(sizeof(lentry));
			cur_trsc->is_coding_exon = (lentry*)malloc_cnt(sizeof(lentry));
			cur_trsc->trsc_type = (char*)malloc_cnt(sizeof(char) * MAX_ATTR_LEN);

			if((cur_trsc == NULL) || (new_entry == NULL) || (cur_trsc->first_exon == NULL)){
				res = 2;
				goto exit;
			}
		
			cur_trsc->chr = cur_exon->chr;
			cur_trsc->strand = strtol(buf_token[idx_of[strand]], &stop, 10);
			cur_trsc->tr_start = strtoul(buf_token[idx_of[trsc_start]], &stop, 10);
			cur_trsc->tr_end   = strtoul(buf_token[idx_of[trsc_end]], &stop, 10);
			cur_trsc->first_cdna_exon = -2;

			// cDNA_start and end
			if(strlen(buf_token[idx_of[cdna_start]]) > 0) {
				cur_trsc->cdna_start = strtoul(buf_token[idx_of[cdna_start]], &stop, 10);
				cur_trsc->first_cdna_exon = strtol(buf_token[idx_of[exon_rank]], &stop, 10) - 1;
;
			}
			else 
				cur_trsc->cdna_start = 0;

			if(strlen(buf_token[idx_of[cdna_end]]) > 0)
				cur_trsc->cdna_end   = strtoul(buf_token[idx_of[cdna_end]], &stop, 10);
			else
				cur_trsc->cdna_end = 0;

			cur_trsc->trsc_length = 0;
			cur_trsc->aa_length = 0;
			cur_trsc->functional_trsc = 0;
			cur_trsc->num_coding_exons = 0;

			cur_trsc->ensembl_id = line_trsc_id;
			cur_trsc->ensembl_protein_id = line_protein_id;
			cur_trsc->first_exon->entry = cur_trsc->first_exon->next = NULL;
			cur_trsc->is_coding_exon->entry = cur_trsc->is_coding_exon->next = NULL;
			strncpy(cur_trsc->trsc_type, buf_token[idx_of[tr_type]], strlen(buf_token[idx_of[tr_type]]));
			cur_trsc->domains = NULL;

			// add to new entry in the front of list
			new_entry->entry = cur_trsc;
			new_entry->next = gen->trsc_head[line_trsc_id % HASH_SIZE];
			gen->trsc_head[line_trsc_id % HASH_SIZE] = new_entry;
			uniq_trsc++;
		}

		// Store cDNA_start
		if(strlen(buf_token[idx_of[cdna_start]]) > 0) {
			if(cur_trsc->cdna_start == 0) {
				cur_trsc->cdna_start = strtoul(buf_token[idx_of[cdna_start]], &stop, 10);
				cur_trsc->first_cdna_exon = strtol(buf_token[idx_of[exon_rank]], &stop, 10) - 1;
			}
			else if((strlen(buf_token[idx_of[cdna_start]]) > 0) && (strtoul(buf_token[idx_of[cdna_start]], &stop, 10) < cur_trsc->cdna_start)) {
				cur_trsc->cdna_start = strtoul(buf_token[idx_of[cdna_start]], &stop, 10);
				cur_trsc->first_cdna_exon = strtol(buf_token[idx_of[exon_rank]], &stop, 10) - 1;
			}
			cur_trsc->num_coding_exons++;
			
			is_coding = 1;
		}
		else {
			is_coding = 0;
		}

		if(!strcmp(cur_trsc->trsc_type, "nonsense_mediated_decay")) {
			is_coding = 0;
		}
		else {
			if(cur_trsc->ensembl_protein_id > 0) {
				cur_exon->is_functional = 1;
			}
		}

		// store cDNA_end
		if((strlen(buf_token[idx_of[cdna_end]]) > 0) && (strtoul(buf_token[idx_of[cdna_end]], &stop, 10) > cur_trsc->cdna_end)) {
			cur_trsc->cdna_end = strtoul(buf_token[idx_of[cdna_end]], &stop, 10);
		}

		// add cur_exon to the cur_trsc
		i = 1;
		cur = cur_trsc->first_exon;
		cur_coding = cur_trsc->is_coding_exon;
		while(i < strtol(buf_token[idx_of[exon_rank]], &stop, 10)) {
			if(cur->next == NULL) {
				new_entry = (lentry*)malloc_cnt(sizeof(lentry));
				new_entry->entry = new_entry->next = NULL;
				cur->next = new_entry;
			}
			if(cur_coding->next == NULL) {
				new_entry = (lentry*)malloc_cnt(sizeof(lentry));
				new_entry->entry = new_entry->next = NULL;
				cur_coding->next = new_entry;
			}
			cur = cur->next;
			cur_coding = cur_coding->next;
			i++;
		}
		cur->entry = cur_exon;
		if(is_coding)	cur_coding->entry = cur_exon;
		else		cur_coding->entry = NULL;
	
		// update transcripts list of exon.
		find_entry_add(cur_exon->first_trsc, cur_trsc);
			
		cur_trsc->trsc_length += (cur_exon->end - cur_exon->start + 1);

		// Lookup gene based on its ID
		if((prev_gene != NULL) && (prev_gene->ensembl_id == line_gene_id))
			cur_gene = prev_gene;
		else
			cur_gene = find_gene_by_id(gen, line_gene_id);

		if(cur_gene == NULL) {
			new_entry = (lentry*)malloc_cnt(sizeof(lentry));
			cur_gene = (gene*)malloc_cnt(sizeof(struct gene));
			cur_gene->symbol = (char*)malloc_cnt(sizeof(char) * MAX_NAME_LEN);
			cur_gene->first_exon = (lentry*)malloc_cnt(sizeof(lentry));
			cur_gene->first_trsc = (lentry*)malloc_cnt(sizeof(lentry));
			if((cur_gene == NULL) || (cur_gene->symbol == NULL) || (new_entry == NULL) ||
				(cur_gene->first_exon == NULL) || (cur_gene->first_trsc == NULL)) {
				res = 2;
				goto exit;
			}

			cur_gene->ensembl_id = line_gene_id;
			strncpy(cur_gene->symbol, buf_token[idx_of[gene_symbol]], strlen(buf_token[idx_of[gene_symbol]])+1);
			make_str_hc(cur_gene->symbol);
			
			cur_gene->chr = cur_exon->chr;
			cur_gene->strand = strtol(buf_token[idx_of[strand]], &stop, 10);
			cur_gene->start = strtoul(buf_token[idx_of[gene_start]], &stop, 10);
			cur_gene->end   = strtoul(buf_token[idx_of[gene_end]], &stop, 10);
			cur_gene->first_trsc->next = cur_gene->first_exon->next = cur_gene->first_trsc->entry = cur_gene->first_exon->entry = NULL;

			// add to new gene in the front of list
			new_entry->entry = cur_gene;
			new_entry->next = gen->gene_head[line_gene_id % HASH_SIZE];
			gen->gene_head[line_gene_id % HASH_SIZE] = new_entry;
			uniq_gene++;
		}
		
		// Add transcript and exon to the gene
		find_entry_add(cur_gene->first_trsc, cur_trsc);
		find_entry_add(cur_gene->first_exon, cur_exon);
		cur_exon->gene = cur_gene;
		cur_trsc->gene = cur_gene;

		prev_trsc = cur_trsc;
		prev_gene = cur_gene;
		cnt++;
		if((cnt % 100000) == 0) {
			printf("%d lines were processed.. (%d genes, %d trscs, %d exons)\n", cnt, uniq_gene, uniq_trsc, uniq_exon);
			/*for(j = 0; j < HASH_SIZE; j++) {
				if((j > 0) && (j % 32 == 0))
					printf("\n");
				printf("%5d", hash_cnt[j]);
			}
			printf("\n");*/
		}
	}
	printf("Finish loading. (%d lines, %d genes, %d trscs, %d exons)\n", cnt, uniq_gene, uniq_trsc, uniq_exon);
exit:
	cnt_trsc_without_CDS = 0;

	// Postprocess entries (it may be deleted)
	for(i = 0; i < HASH_SIZE; i++) {
		int utr5;
		exon* next_exon;
		char* trsc_seq;

		outer = gen->exon_head[i];
		while((outer != NULL) && (outer->entry != NULL)) {
			cur_exon = (exon*)outer->entry;
			if((cur_exon->start < 30) || (cur_exon->end > cur_exon->chr->length - 30)) {
				fprintf(stderr, "Exon %011u is on boundary of chromosome\n", cur_exon->ensembl_id);
			}
			outer = outer->next;
		}

		outer = gen->trsc_head[i];
		while((outer != NULL) && (outer->entry != NULL)) {
			cur_trsc = (transcript*)outer->entry;

			cur_trsc->num_exons = 0;
			utr5 = cur_trsc->cdna_start;
			cur = cur_trsc->first_exon;
			while((cur != NULL) && (cur->entry != NULL)) {
				cur_exon = (exon*)cur->entry;
				// Update end phase information here
				if(cur_exon->end_phase == -2) {	// GRCh37 data
					if((utr5 == 1) && (cur_exon->start_phase == -1)) // coding region is start at the first base of a exon (no 5'UTR in this exon)
						cur_exon->start_phase = 0;

					if(cur->next != NULL) {
						next_exon = (exon*)cur->next->entry;
						if(((cur_trsc->num_exons + 2) == cur_trsc->first_cdna_exon) && (next_exon->start_phase == 0))	//next of cur_exon is the first coding exon and its start phase is 0. first_cdna_exon is 1-based number.
							cur_exon->end_phase = -1;
						else
							cur_exon->end_phase = next_exon->start_phase;
					}
					else	
						cur_exon->end_phase = -1;
				}

				utr5 -= cur_exon->end - cur_exon->start + 1;
				cur = cur->next;
				cur_trsc->num_exons++;
			}

			//if(cur_trsc->ensembl_protein_id > 0)
			//	cur_trsc->aa_length = calculate_aa_len(gen, cur_trsc);
			

			// Check transcripts if they have correct start and stop codons
			trsc_seq = concat_exon_seqs(gen, cur_trsc, NO_BUFFER);
			if( 
		            (cur_trsc->cdna_start > 0) &&
			    (
			      ((trsc_seq[cur_trsc->cdna_start-1] == 'a') || (trsc_seq[cur_trsc->cdna_start-1] == 'A')) &&
			      ((trsc_seq[cur_trsc->cdna_start]   == 't') || (trsc_seq[cur_trsc->cdna_start]   == 'T')) && 
			      ((trsc_seq[cur_trsc->cdna_start+1] == 'g') || (trsc_seq[cur_trsc->cdna_start+1] == 'G'))
			    )
			  ) {
				// Start codon is okay. Then check stop codon				
				if( ( (trsc_seq[cur_trsc->cdna_end-3] == 't') || (trsc_seq[cur_trsc->cdna_end-3] == 'T') ) &&
	 			    ( (trsc_seq[cur_trsc->cdna_end-2] == 'a') || (trsc_seq[cur_trsc->cdna_end-2] == 'A') ) &&
				    ( (trsc_seq[cur_trsc->cdna_end-1] == 'g') || (trsc_seq[cur_trsc->cdna_end-1] == 'G') ) ) {
				    	// TAG. Okay
					cur_trsc->has_CDS = 1;
				}
				else if( ( (trsc_seq[cur_trsc->cdna_end-3] == 't') || (trsc_seq[cur_trsc->cdna_end-3] == 'T') ) &&
				  	 ( (trsc_seq[cur_trsc->cdna_end-2] == 'g') || (trsc_seq[cur_trsc->cdna_end-2] == 'G') ) &&
					 ( (trsc_seq[cur_trsc->cdna_end-1] == 'a') || (trsc_seq[cur_trsc->cdna_end-1] == 'A') ) ) {
					// TGA. Okay
					cur_trsc->has_CDS = 1;
				}
				else if( ( (trsc_seq[cur_trsc->cdna_end-3] == 't') || (trsc_seq[cur_trsc->cdna_end-3] == 'T') ) &&
					 ( (trsc_seq[cur_trsc->cdna_end-2] == 'a') || (trsc_seq[cur_trsc->cdna_end-2] == 'A') ) &&
					 ( (trsc_seq[cur_trsc->cdna_end-1] == 'a') || (trsc_seq[cur_trsc->cdna_end-1] == 'A') ) ) {
					// TAA. Okay
					cur_trsc->has_CDS = 1;
				}
				else {
					// Proper start codon, but no stop codon.
					cnt_trsc_without_CDS++;
					cur_trsc->has_CDS = 0;
				}
			}
			else {
				// no start codon
				cnt_trsc_without_CDS++;
				cur_trsc->has_CDS = 0;
			}

			free_cnt(trsc_seq);

			outer = outer->next;
		}

		outer = gen->gene_head[i];
		while((outer != NULL) && (outer->entry != NULL)) {
			cur_gene = (gene*)outer->entry;
			create_and_add_entry_last(cur_gene->chr->genes, cur_gene);
			cur_gene->chr->num_genes++;
			cur_gene->num_exons = cur_gene->num_transcripts = cur_gene->n_functional_trsc = 0;
			cur = cur_gene->first_exon;
			while((cur != NULL) && (cur->entry != NULL)) {
				cur = cur->next;
				cur_gene->num_exons++;
			}
			
			cur = cur_gene->first_trsc;
			while((cur != NULL) && (cur->entry != NULL)) {
				cur_trsc = (transcript*)cur->entry;
				// Number of coding transctions --> counted in crispinator.c
				//if(cur_trsc->ensembl_protein_id > 0)
				//	cur_gene->n_coding_trsc++;
				cur = cur->next;
				cur_gene->num_transcripts++;
			}
			outer = outer->next;
		}
	}

	printf("%d transcripts were marked due to lack of proper start codon.\n", cnt_trsc_without_CDS);

	for(i = 0; i < header_size; i++) {
		if(buf_token[i] != NULL) 
			free_cnt(buf_token[i]);
	}
	if(buf_token != NULL)
		free_cnt(buf_token);
	if(idx_of != NULL)
		free_cnt(idx_of);
	if(tofree != NULL)
		free_cnt(tofree);

#ifdef DEBUG
	printf("After load_ensembl_data() ==> ");
	print_cnt();
#endif

	return res;
}


// Load ensemble domain information. Also column names are hard-coded.
// This function should be called after load_ensembl_data() was called.
int load_ensembl_domain(genome* gen, char* path) {
	FILE *ensembl_file;
	char *buf, *temp_buf, *token, *stop, *tofree, *ptr;
	char temp[MAX_ATTR_LEN];
	char** buf_token;
	int*   idx_of;
	transcript* cur_trsc = NULL, *prev_trsc = NULL;
	domain* cur_dom = NULL;

	int header_size, num_of_features = 8;
	unsigned int line_trsc_id, line_interPro_id;
	int gene_id = 0, trsc_id = 1, interPro_id = 2, description = 3,
		start = 4, end = 5, hmm_start = 6, hmm_end = 7, i, cnt = 0, uniq_trsc = 0, found;

	lentry* cursor;

	printf("Loading ensembl protein domain data..\n");
	ensembl_file = fopen(path, "r");
	if(ensembl_file == NULL) {
		fprintf(stderr, "Cannot open ensembl file (%s)\n", path);
		return 1;
	}

	buf = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
	temp_buf = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
	tofree = buf;

	header_size = 0;
	fgets(buf, BUF_LEN, ensembl_file);
	while((token = strsep(&buf, "\t")) != NULL)
		header_size++;

	fseek(ensembl_file, 0, SEEK_SET);
	printf("Number of columns = %d\n", header_size);	
	
	buf = tofree;	
	idx_of = (int*)malloc_cnt(sizeof(int) * num_of_features);
	if(idx_of == NULL)
		return 2;

	buf_token = (char**)malloc_cnt(sizeof(char*) * header_size);
	if((buf == NULL) || (buf_token == NULL)) {
		return 2;
	}

	for(i = 0; i < num_of_features; i++) {
		idx_of[i] = -1;
	}
	for(i = 0; i < header_size; i++) {
		buf_token[i] = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
		if(buf_token[i] == NULL) {return 2;
			return 2;
		}
	}

	fgets(buf, BUF_LEN, ensembl_file);
	i = -1;
	if(buf[strlen(buf) - 1] == '\n')
		buf[strlen(buf) - 1] = '\0';

	while((i < header_size) && ((token = strsep(&buf, "\t")) != NULL)) {
		i++;
		if(strcmp(token, "Gene ID") == 0)				idx_of[gene_id] = i;
		else if(strcmp(token, "Transcript ID") == 0)			idx_of[trsc_id] = i;
		else if (strcmp(token, "Interpro ID") == 0)				idx_of[interPro_id] = i;
		else if (strcmp(token, "Interpro Description") == 0)			idx_of[description] = i;
		else if (strcmp(token, "Interpro start") == 0)				idx_of[start] = i;
		else if (strcmp(token, "Interpro end") == 0)				idx_of[end] = i;
		else if (strcmp(token, "HMMPanther start") == 0)			idx_of[hmm_start] = i;
		else if (strcmp(token, "HMMPanther end") == 0)				idx_of[hmm_end] = i;
	} 

	buf = tofree;

	for(i = 0; i < num_of_features; i++) {
		if(idx_of[i] == -1) {
			fprintf(stderr, "Some header(%d th) are not found!\n", i);
			fprintf(stderr, "It fixed to 0\n");
			idx_of[i] = 0;
		}
	}

	// load ensembl data
	while(fgets(buf, BUF_LEN, ensembl_file)) {
		// remove last newline character
		if(buf[strlen(buf) - 1] == '\n')
			buf[strlen(buf) - 1] = '\0';

		strncpy(temp_buf, buf, strlen(buf));
		temp_buf[strlen(buf) - 1] = '\0';

		i = 0;
		while((i < header_size) && ((token = strsep(&buf, "\t")) != NULL)) {
			strncpy(buf_token[i], token, strlen(token));
			ptr = buf_token[i];
			ptr[strlen(token)] = '\0';
			//printf("Token %d = %s\n", i, buf_token[i]);
			i++;
		}

		if(i < header_size) {
			buf = tofree;
			continue;
		}

		if(!strcmp(buf_token[idx_of[hmm_start]], buf_token[idx_of[start]]) &&
		   !strcmp(buf_token[idx_of[hmm_end]], buf_token[idx_of[end]])) {
			buf = tofree;
			continue;
		}

		buf = tofree;
		ptr = buf_token[idx_of[interPro_id]];
		if(strlen(ptr) > 0) {
			strncpy(temp, ptr + 3, strlen(ptr) - 3);	
			temp[strlen(ptr)-3] = '\0';
			line_interPro_id = strtoul(temp, &stop, 10);
		}
		else {
			line_interPro_id = 0;
		}
		
		ptr = buf_token[idx_of[trsc_id]];
		if(strlen(ptr) > 0) {
			strncpy(temp, ptr + 4, strlen(ptr) - 4);	
			temp[strlen(ptr)-4] = '\0';
			line_trsc_id = strtoul(temp, &stop, 10);
		}
		else {
			line_trsc_id = 0;
		}

		//printf("%s --> %u\n", temp_buf, line_trsc_id);

		if((line_trsc_id == 0) || (line_interPro_id == 0)) {
			//printf("ID NULL: %s\n", temp_buf);
			continue;
		}

		if((prev_trsc != NULL) && (prev_trsc->ensembl_id == line_trsc_id))
			cur_trsc = prev_trsc;
		else {
			cur_trsc = find_trsc_by_id(gen, line_trsc_id);
			uniq_trsc++;
		}

		if(cur_trsc == NULL) {
			//printf("Cannot find transcript: [%s]\n", temp_buf);
			continue;
		}

		found = 0;
		int new_start, new_end;
		new_start = strtol(buf_token[idx_of[start]], &stop, 10);
		new_end = strtol(buf_token[idx_of[end]], &stop, 10);
		//printf("For domain %d (%d-%d)\n", line_interPro_id, new_start, new_end);
		FOR_LENTRY(cursor, cur_trsc->domains) {
			domain* cur_dom = (domain*)cursor->entry;
			//printf("  exist domain: %d (%d-%d)\n", cur_dom->interProID, cur_dom->start, cur_dom->end);
			if( ( cur_dom->interProID == line_interPro_id ) && 
			    ( if_any_overlap(new_start, new_end, cur_dom->start, cur_dom->end ) ) ){
				
				if (cur_dom->end - cur_dom->start < new_end - new_start) {
					cur_dom->start = new_start;
					cur_dom->end = new_end;
				}

				found = 1;
				break;
			}
		}

		if(found == 1) {
			continue;
		}
		
		cur_dom = (domain*)malloc_cnt(sizeof(domain));
		cur_dom->description = (char*)malloc_cnt(sizeof(char) * MAX_ATTR_LEN);
		cur_dom->interProID = line_interPro_id;
		cur_dom->start = strtol(buf_token[idx_of[start]], &stop, 10);
		cur_dom->end   = strtol(buf_token[idx_of[end]], &stop, 10);
		strncpy(cur_dom->description, buf_token[idx_of[description]], strlen(buf_token[idx_of[description]]));
		cur_dom->description[strlen(buf_token[idx_of[description]])] = '\0';
			
		if(cur_trsc->domains == NULL) {
			cur_trsc->domains = (lentry*)malloc_cnt(sizeof(lentry));
			cur_trsc->domains->next = NULL;
			cur_trsc->domains->entry = cur_dom;
		}
		else {
			find_entry_add(cur_trsc->domains, cur_dom);
		}

		{
			//printf("Line: %s\n", temp_buf);
			//printf("TRSC: ENST%011d (AA len: %u)\n", cur_trsc->ensembl_id, cur_trsc->aa_length);
			lentry* cursor = cur_trsc->domains;
			while(cursor != NULL) {
				cur_dom = (domain*)cursor->entry;
				//printf("\tIPR%06d\t%d-%d\t%s\n", cur_dom->interProID, cur_dom->start, cur_dom->end, cur_dom->description);
				cursor = cursor->next;
			}
			//getchar();
		}

		cnt++;
		if((cnt % 100000) == 0) {
			printf("%d lines were processed.. (%d trsc)\n", cnt, uniq_trsc);
		}

		prev_trsc = cur_trsc;
	}

	free_cnt(buf);
	free_cnt(temp_buf);
	for(i = 0; i < header_size; i++) {
		free_cnt(buf_token[i]);
	}
	free_cnt(buf_token);
	free_cnt(idx_of);

	return 0;
}

// Free all struct
void free_genome(genome* gen) {
	lentry *cur, *inner;
	exon* cur_exon;
	gene* cur_gene;
	transcript* cur_trsc;
	domain* cur_dom;

#ifdef DEBUG
	printf("Before free_genome() ==> ");
	print_cnt();
#endif

	if(gen != NULL) {
		int i;
		for(i = 0; i < gen->num_chrs; i++) {
			if(gen->chrs[i]->genes != NULL)
				free_cnt(gen->chrs[i]->genes);
			if(gen->chrs[i]->name != NULL)
				free_cnt(gen->chrs[i]->name);
			if(gen->chrs[i] != NULL)
				free_cnt(gen->chrs[i]);
		}
		for(i = 0; i < HASH_SIZE; i++) {
			while(gen->exon_head[i] != NULL) {
				cur = gen->exon_head[i];
				gen->exon_head[i] = gen->exon_head[i]->next;
				cur_exon = (exon*)cur->entry;
				free_cnt(cur);
				if(cur_exon != NULL) {
					free_cnt(cur_exon);
				}
			}
		}
		if(gen->exon_head != NULL) {
			free_cnt(gen->exon_head);
		}
	
		for(i = 0; i < HASH_SIZE; i++) {
			while(gen->trsc_head[i] != NULL) {
				cur = gen->trsc_head[i];
				gen->trsc_head[i] = gen->trsc_head[i]->next;
				cur_trsc = (transcript*)cur->entry;
				free_cnt(cur);

				if (cur_trsc != NULL) {
					free_lentry_only(cur_trsc->first_exon);
					free_lentry_only(cur_trsc->is_coding_exon);
					free_cnt(cur_trsc->trsc_type);
					if(cur_trsc->domains != NULL) {
						inner = cur_trsc->domains;
						while(inner != NULL) {
							cur_dom = (domain*)inner->entry;
							free_cnt(cur_dom->description);
							free_cnt(cur_dom);
							inner = inner->next;
						}
					}
					free_cnt(cur_trsc);
				}
			}
		}
		if(gen->trsc_head != NULL) {
			free_cnt(gen->trsc_head);
		}

		for(i = 0; i < HASH_SIZE; i++) {
			while(gen->gene_head[i] != NULL) {
				cur = gen->gene_head[i];
				gen->gene_head[i] = gen->gene_head[i]->next;
				cur_gene = (gene*)cur->entry;
				free_cnt(cur);

				if(cur_gene != NULL) {
					if(cur_gene->symbol != NULL)
						free_cnt(cur_gene->symbol);
					free_lentry_only(cur_gene->first_exon);
					free_lentry_only(cur_gene->first_trsc);
					free_cnt(cur_gene);
				}
			}
		}
		if(gen->gene_head != NULL) {
			free_cnt(gen->gene_head);
		}

		if(gen->chrs != NULL) {
			for(i = 0; i < gen->num_chrs; i++)
				free_lentry_only(gen->chrs[i]->genes);
			free_cnt(gen->chrs);
		}
		fclose(gen->fasta);
		free_cnt(gen->fasta_path);
		free_cnt(gen);
	}

#ifdef DEBUG
	printf("After free_genome() ==> ");
	print_cnt();
#endif
}
