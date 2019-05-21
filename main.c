/*
 * main.c     CRISPinatoR
 *
 *  Date: 2016. 6. 27.
 *  Author: Yunku Yeu
 */

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "kstring.h"
#include "crispinator.h"
#include "tools.h"
#include "bwtindex.h"


static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   CRISPinatoR <command> [options]\n\n");
	fprintf(stderr, "Command: index                   index sequences in the FASTA format\n");
	fprintf(stderr, "         offtarget               Run an interactive off-target test module\n");
	fprintf(stderr, "         library                 Write sgRNA list for whole genome\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char* argv[]) {
	int ret = 0;

	if (argc < 2) return usage();

	if      (strcmp(argv[1], "index") == 0)	{ret = bwa_index(argc-1, argv+1);}
	else if (strcmp(argv[1], "offtarget") == 0)	{ret = run_interactive_offtarget(argc-1, argv+1);}
	else if (strcmp(argv[1], "library") == 0) 	{ret = crispinator_write_library(argc-1, argv+1);}
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}

	err_fflush(stdout);
	err_fclose(stdout);

	return ret;
}
