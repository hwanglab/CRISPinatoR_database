#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwtgap.h"
#include "bwtaln.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define aln_score(m,o,e,p) ((m)*(p)->s_mm + (o)*(p)->s_gapo + (e)*(p)->s_gape)

gap_stack_t *gap_init_stack2(int max_score)
{
	gap_stack_t *stack;
	stack = (gap_stack_t*)calloc(1, sizeof(gap_stack_t));
	stack->n_stacks = max_score;
	stack->stacks = (gap_stack1_t*)calloc(stack->n_stacks, sizeof(gap_stack1_t));
	return stack;
}

gap_stack_t *gap_init_stack(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt)
{
	return gap_init_stack2(aln_score(max_mm+1, max_gapo+1, max_gape+1, opt));
}

void gap_destroy_stack(gap_stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_stacks; ++i) free(stack->stacks[i].stack);
	free(stack->stacks);
	free(stack);
}

static void gap_reset_stack(gap_stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_stacks; ++i)
		stack->stacks[i].n_entries = 0;
	stack->best = stack->n_stacks;
	stack->n_entries = 0;
}

static inline void gap_push(gap_stack_t *stack, int i, bwtint_t k, bwtint_t l, int n_mm, int n_gapo, int n_gape, int n_ins, int n_del,
							int state, int is_diff, const gap_opt_t *opt)
{
	int score;
	gap_entry_t *p;
	gap_stack1_t *q;
	score = aln_score(n_mm, n_gapo, n_gape, opt);
	q = stack->stacks + score;
	if (q->n_entries == q->m_entries) {
		q->m_entries = q->m_entries? q->m_entries<<1 : 4;
		q->stack = (gap_entry_t*)realloc(q->stack, sizeof(gap_entry_t) * q->m_entries);
	}
	p = q->stack + q->n_entries;
	p->info = (u_int32_t)score<<21 | i; p->k = k; p->l = l;
	p->n_mm = n_mm; p->n_gapo = n_gapo; p->n_gape = n_gape;
	p->n_ins = n_ins; p->n_del = n_del;
	p->state = state; 
	p->last_diff_pos = is_diff? i : 0;
	++(q->n_entries);
	++(stack->n_entries);
	if (stack->best > score) stack->best = score;
}

static inline void gap_pop(gap_stack_t *stack, gap_entry_t *e)
{
	gap_stack1_t *q;
	q = stack->stacks + stack->best;
	*e = q->stack[q->n_entries - 1];
	--(q->n_entries);
	--(stack->n_entries);
	if (q->n_entries == 0 && stack->n_entries) { // reset best
		int i;
		for (i = stack->best + 1; i < stack->n_stacks; ++i)
			if (stack->stacks[i].n_entries != 0) break;
		stack->best = i;
	} else if (stack->n_entries == 0) stack->best = stack->n_stacks;
}

static inline void gap_shadow(int x, int len, bwtint_t max, int last_diff_pos, bwt_width_t *w)
{
	int i, j;
	for (i = j = 0; i < last_diff_pos; ++i) {
		if (w[i].w > x) w[i].w -= x;
		else if (w[i].w == x) {
			w[i].bid = 1;
			w[i].w = max - (++j);
		} // else should not happen
	}
}

static inline int int_log2(uint32_t v)
{
	int c = 0;
	if (v & 0xffff0000u) { v >>= 16; c |= 16; }
	if (v & 0xff00) { v >>= 8; c |= 8; }
	if (v & 0xf0) { v >>= 4; c |= 4; }
	if (v & 0xc) { v >>= 2; c |= 2; }
	if (v & 0x2) c |= 1;
	return c;
}

// This function was modified to find off-targets for a sgRNA.
// It assumes that the input sequences includes PAM sequence and were converted into NGG form already. (No CCN sgRNA sequence)
bwt_aln1_t *bwt_match_gap(bwt_t *const bwt, int len, const ubyte_t *seq, bwt_width_t *width,
						  bwt_width_t *seed_width, const gap_opt_t *opt, int *_n_aln, gap_stack_t *stack)
{ // $seq is the reverse complement of the input read
	//int best_score = aln_score(opt->max_diff+1, opt->max_gapo+1, opt->max_gape+1, opt);
	int best_score = aln_score(opt->max_diff+1, 0, 0, opt);	// we didn't consider gap
	int best_diff = opt->max_diff + 1, max_diff = opt->max_diff;
	int best_cnt = 0;
	int max_entries = 0, j, _j, n_aln, m_aln;
	bwt_aln1_t *aln;

	m_aln = 4; n_aln = 0;
	aln = (bwt_aln1_t*)calloc(m_aln, sizeof(bwt_aln1_t));

	// check whether there are too many N
	for (j = _j = 0; j < len; ++j)
		if (seq[j] > 3) ++_j;
	if (_j > max_diff) {
		*_n_aln = n_aln;
		return aln;
	}

	//for (j = 0; j != len; ++j) printf("#0 %d: [%d,%u]\t[%d,%u]\n", j, w[0][j].bid, w[0][j].w, w[1][j].bid, w[1][j].w);
	gap_reset_stack(stack); // reset stack
	gap_push(stack, len, 0, bwt->seq_len, 0, 0, 0, 0, 0, 0, 0, opt);

	while (stack->n_entries) {
		gap_entry_t e;
		int i, m, m_seed = 0, hit_found, allow_diff, allow_M;
		bwtint_t k, l, cnt_k[4], cnt_l[4];//, occ;

		if (max_entries < stack->n_entries) max_entries = stack->n_entries;
		if (stack->n_entries > opt->max_entries) break;
		gap_pop(stack, &e); // get the best entry
		k = e.k; l = e.l; // SA interval
		i = e.info&0xffff; // length
		//if (!(opt->mode & BWA_MODE_NONSTOP) && e.info>>21 > best_score + opt->s_mm) break; // no need to proceed

		m = max_diff - e.n_mm;
		if(i >= 2)	// increase m by 1 for N of NGG.
			m++;
		else if((i < 2) && (m > 0))	// no difference for GG of NGG
			m = 0;

		if (m < 0) continue;
		if (seed_width) { // apply seeding
			m_seed = opt->max_seed_diff - (e.n_mm + e.n_gapo);
			if (opt->mode & BWA_MODE_GAPE) m_seed -= e.n_gape;
		}
		
		// modified by yyk.
		// seq[0-2] is PAM sequence. seq[0-1] is 'GG', and mismatch on seq[2] will be ignored.
		if (i > 0 && m < width[i-1].bid) continue;

		// check whether a hit is found
		hit_found = 0;
		if ((i == 0) && (e.n_mm <= max_diff)) hit_found = 1;

		if (hit_found) { // action for found hits
			int score = aln_score(e.n_mm, e.n_gapo, e.n_gape, opt);
			//printf("#2 hits found: %d:(%u,%u)\n", e.n_mm+e.n_gapo, k, l);
			if (n_aln == 0) {
				best_score = score;
				best_diff = e.n_mm + e.n_gapo;
				if (opt->mode & BWA_MODE_GAPE) best_diff += e.n_gape;
			}
			if (score == best_score) best_cnt += l - k + 1;
			
			bwt_aln1_t *p;
			gap_shadow(l - k + 1, len, bwt->seq_len, e.last_diff_pos, width);
			if (n_aln == m_aln) {
				m_aln <<= 1;
				aln = (bwt_aln1_t*)realloc(aln, m_aln * sizeof(bwt_aln1_t));
				memset(aln + m_aln/2, 0, m_aln/2*sizeof(bwt_aln1_t));
			}
			p = aln + n_aln;
			p->n_mm = e.n_mm; p->n_gapo = e.n_gapo; p->n_gape = e.n_gape;
			p->n_ins = e.n_ins; p->n_del = e.n_del;
			p->k = k; p->l = l;
			p->score = score;
			//fprintf(stderr, "*** n_mm=%d,n_gapo=%d,n_gape=%d,n_ins=%d,n_del=%d\n", e.n_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del);
			++n_aln;

			continue;
		}

		--i;
		bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l); // retrieve Occ values
		//occ = l - k + 1;
		// test whether diff is allowed
		allow_diff = allow_M = 1;
		if (i < 2)	// no mismatch in the PAM sequences (NGG)
			allow_diff = allow_M = 0;
		else if (i > 0) {
			int ii = i - (len - opt->seed_len);
			if (width[i-1].bid > m - 1) allow_diff = 0;
			else if (width[i-1].bid == m - 1 && width[i].bid == m - 1 && width[i-1].w == width[i].w) allow_M = 0;
			if (seed_width && ii > 0) {
				printf("Enter here?\n");
				if (seed_width[ii-1].bid > m_seed-1) allow_diff = 0;
				else if (seed_width[ii-1].bid == m_seed-1 && seed_width[ii].bid == m_seed-1
						 && seed_width[ii-1].w == seed_width[ii].w) allow_M = 0;
			}
		}
		// indel is deleted by yyk
		// mismatches
		if ((i == 2) || (allow_diff && allow_M)) { // mismatch is allowed
			for (j = 1; j <= 4; ++j) {
				int c = (seq[i] + j) & 3;
				int is_mm = (j != 4 || seq[i] > 3);
				k = bwt->L2[c] + cnt_k[c] + 1;
				l = bwt->L2[c] + cnt_l[c];
				if (k <= l) {
					if((i == 2) && (seq[i] <= 3)) // any nt at 'N'GG is not counted as a misatch
						gap_push(stack, i, k, l, e.n_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del, STATE_M, 0, opt);
					else
						gap_push(stack, i, k, l, e.n_mm + is_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del, STATE_M, is_mm, opt);
				}
			}
		} else if (seq[i] < 4) { // try exact match only
			int c = seq[i] & 3;
			k = bwt->L2[c] + cnt_k[c] + 1;
			l = bwt->L2[c] + cnt_l[c];
			if (k <= l) gap_push(stack, i, k, l, e.n_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del, STATE_M, 0, opt);
		}
	}

	*_n_aln = n_aln;
	//fprintf(stderr, "max_entries = %d\n", max_entries);
	return aln;
}
