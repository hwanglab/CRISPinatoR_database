#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "motif.h"

queue* init_queue(void) {
	queue* q = (queue*)malloc(sizeof(queue));
	q->first = q->last = NULL;
	q->length = 0;
	return q;
}

void free_queue(queue* q) {
	node* n;
	while(q->first != NULL) {
		n = q->first;
		q->first = q->first->next;
		free(n);
	}
	free(q);
}

void enqueue(queue* q, pwm* motif, int new_val) {
	node *n = (node*)malloc(sizeof(node));
	n->value = new_val;
	n->motif = motif;
	n->next  = NULL;
	if(q->length == 0)
		q->first = q->last = n;
	else {
		q->last->next = n;
		q->last = n;
	}
	q->length++;
}

node* dequeue_nofree(queue* q) {
	node *n;

	if(q->length == 0)
		return NULL;
	n = q->first;
	q->first = q->first->next;
	q->length--;
	
	return n;
}

int dequeue(queue* q) {
	int val;
	node* n;
	if(q->length == 0)
		return -1000000;
	n = q->first;
	q->first = q->first->next;
	q->length--;
	
	val = n->value;
	free(n);
	return val;
}

// load ESS or ESS list
kmers* load_kmers(char* path) {
	int i;
	FILE* input = NULL;
	char *buf, *token, *tobuf, *stop;
	kmers* set;
	
	buf = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
	tobuf = buf;

	input = fopen(path, "r");
	if(input == NULL) {
		fprintf(stderr, "Cannot open %s\n", path);
		return NULL;
	}

	set = (kmers*)malloc_cnt(sizeof(kmers));
	set->length = 0;
		
	while(fgets(buf, BUF_LEN, input) != NULL) {
		if(strlen(buf) > 0) 
			set->length++;
	}
	
	set->mers    = (char**)malloc_cnt(sizeof(char*) * set->length);
	set->notes   = (char**)malloc_cnt(sizeof(char*) * set->length);
	set->k_size  = (int*)malloc_cnt(sizeof(int) * set->length);
	set->weights = (double*)malloc_cnt(sizeof(double) * set->length);

	fseek(input, 0, SEEK_SET);

	i = 0;
	while(fgets(buf, BUF_LEN, input) != NULL) {
		if(strlen(buf) > 0) {
			while((buf[strlen(buf) - 1] == '\n') || (buf[strlen(buf) - 1] == '\r'))
				buf[strlen(buf) - 1] = '\0';

			token = strsep(&buf, "\t");

			set->k_size[i] = strlen(token);
			while((token[set->k_size[i] - 1] == '\n') || (token[set->k_size[i] - 1] == '\r'))
				set->k_size[i]--;
			set->mers[i] = (char*)malloc_cnt(sizeof(char) * (set->k_size[i] + 1));
			strncpy(set->mers[i], token, set->k_size[i]);
			set->mers[i][set->k_size[i]] = '\0';

			token = strsep(&buf, "\t");
			if(strlen(token) > 0) 
				set->weights[i] = strtod(token, &stop);
			else
				set->weights[i] = 0;

			token = strsep(&buf, "\t");
			set->notes[i] = (char*)malloc_cnt(sizeof(char) * 24);
			strncpy(set->notes[i], token, strlen(token));
			set->notes[i][strlen(token)] = '\0';
			
			i++;		
			buf = tobuf;
		}
	}

	fclose(input);
	free_cnt(buf);

	return set;
}

void free_kmers(kmers* set) {
	int i;
	for(i = 0; i < set->length; i++) {
		free_cnt(set->mers[i]);
		free_cnt(set->notes[i]);
	}
	free_cnt(set->weights);
	free_cnt(set->mers);
	free_cnt(set->k_size);
	free_cnt(set);
}

void compare_print_mers_list(FILE* stream, kmers* list, char* mer1, char* mer2, int tap) {
	int i, j, size1 = 0, size2 = 0, *mer1_int, *mer2_int, exist;
	char *token, *stop;

	for(i = 0; i < strlen(mer1); i++) 
		if(mer1[i] == '/')
			size1++;
	if(strlen(mer1) > 0)
		size1++;

	for(i = 0; i < strlen(mer2); i++) 
		if(mer2[i] == '/')
			size2++;
	if(strlen(mer2) > 0)
		size2++;

	mer1_int = (int*)malloc_cnt(sizeof(int) * size1);
	mer2_int = (int*)malloc_cnt(sizeof(int) * size2);

	for(i = 0; i < size1; i++) {
		token = strsep(&mer1, "/");
		mer1_int[i] = strtol(token, &stop, 10);
	}

	for(i = 0; i < size2; i++) {
		token = strsep(&mer2, "/");
		mer2_int[i] = strtol(token, &stop, 10);
	}

	for(j = 0; j < tap; j++) 
		fputc(' ', stream);
	fprintf(stream, "            Original                           Modified\n");
	for(j = 0; j < tap; j++) 
		fputc(' ', stream);
	fprintf(stream, " Motif    Weight       Sources          Motif    Weight       Sources\n");
	for(j = 0; j < tap; j++) 
		fputc(' ', stream);
	fprintf(stream, "-------------------------------------------------------------------------\n");

	for(i = 0; i < size1; i++) {
		for(j = 0; j < tap; j++) 
			fputc(' ', stream);

		fprintf(stream, "%-10s", list->mers[mer1_int[i]]);
		if(list->weights[mer1_int[i]] != 0)
			fprintf(stream, "%1.4f", list->weights[mer1_int[i]]);
		else
			fprintf(stream, "  NA  ");
		fprintf(stream, "%16s", list->notes[mer1_int[i]]);

		exist = 0;
		for(j = 0; (j < size2) && (!exist); j++) {
			if(mer1_int[i] == mer2_int[j]) {
				exist = 1;
				mer2_int[j] = -1;
			}
		}

		if(exist) {
			fprintf(stream, "\t-> preserved\n");
		}
		else {
			fprintf(stream, "\t-> destroyed\n");
		}
	}
	for(i = 0; i < size2; i++) {
		if(mer2_int[i] >= 0) {
			for(j = 0; j < tap; j++) 
				fputc(' ', stream);
			fprintf(stream, "                               \t-> %-10s", list->mers[mer2_int[i]]);
			if(list->weights[mer2_int[i]] != 0)
				fprintf(stream, "%1.4f", list->weights[mer2_int[i]]);
			else
				fprintf(stream, "  NA  ");
			fprintf(stream, "%16s\n", list->notes[mer2_int[i]]);
		}
	}
	free_cnt(mer1_int);
	free_cnt(mer2_int);
}

pwm_list* create_pwm_from_file(const char* path) {
	int i, j, k;
	//double tmp_min, tmp_max, min, max, ratio = 0.9;
	FILE* input = NULL;
	char *buf, *token, *tobuf, *stop;

	input = fopen(path, "r");
	if(input == NULL) {
		fprintf(stderr, "Cannot open %s\n", path);
		return NULL;
	}

	buf = (char*)malloc_cnt(sizeof(char) * BUF_LEN);
	tobuf = buf;

	pwm_list* list = (pwm_list*)malloc_cnt(sizeof(pwm_list) * 1);
	list->length = 0;
		
	while(fgets(buf, BUF_LEN, input) != NULL) {
		if(strlen(buf) > 0) 
			list->length++;
	}
	
	list->length /= 5;	// 1 for name and length, 4 for each NT.

	//printf("Number of PWMs = %d\n", list->length);

	list->pwms = (pwm*)malloc_cnt(sizeof(pwm) * list->length);
	
	fseek(input, 0, SEEK_SET);

	i = 0;
	while(fgets(buf, BUF_LEN, input) != NULL) {
		if(strlen(buf) > 0) {
			while((buf[strlen(buf) - 1] == '\n') || (buf[strlen(buf) - 1] == '\r'))
				buf[strlen(buf) - 1] = '\0';

			// name of PWM
			list->pwms[i].name = (char*)malloc_cnt(sizeof(char) * 12);
			token = strsep(&buf, "\t");
			strncpy(list->pwms[i].name, token, strlen(token));
			list->pwms[i].name[strlen(token)] = '\0';
			
			// length of PWM
			token = strsep(&buf, "\t");
			list->pwms[i].length = strtol(token, &stop, 10);

			// length of PWM
			token = strsep(&buf, "\t");
			list->pwms[i].threshold = strtod(token, &stop);
			
			//printf("PWM: %s (len=%d)", list->pwms[i].name, list->pwms[i].length);

			list->pwms[i].weights = (double**)malloc_cnt(sizeof(double*) * list->pwms[i].length);
			for(j = 0; j < list->pwms[i].length; j++) 
				list->pwms[i].weights[j] = (double*)malloc_cnt(sizeof(double) * 4);

			for(k = 0; k < 4; k++) {
				buf = tobuf;
				if(fgets(buf, BUF_LEN, input) == NULL) {
					fprintf(stderr, "Unexpected PWM format.\n");
					return NULL;
				}
				while((buf[strlen(buf) - 1] == '\n') || (buf[strlen(buf) - 1] == '\r'))
					buf[strlen(buf) - 1] = '\0';

				for(j = 0; j < list->pwms[i].length; j++) {
					token = strsep(&buf, "\t");
					if(token == NULL) {
						fprintf(stderr, "Unexpected PWM format.\n");
						return NULL;
					}
					list->pwms[i].weights[j][k] = strtod(token, &stop);
				}
			}

			/*min = max = 0.0;
			for(j = 0; j < list->pwms[i].length; j++) {
				tmp_min = 999999.0;
				tmp_max = -999999.0;
				for(k = 0; k < 4; k++) {
					if(list->pwms[i].weights[j][k] < tmp_min)
						tmp_min = list->pwms[i].weights[j][k];
					if(list->pwms[i].weights[j][k] > tmp_max)
						tmp_max = list->pwms[i].weights[j][k];
				}
				min += tmp_min;
				max += tmp_max;
			}

			list->pwms[i].threshold = (max - min) * ratio;*/
			//printf(" Thres = %f\n", list->pwms[i].threshold);
	
			i++;		
			buf = tobuf;
		}
	}

	fclose(input);
	free_cnt(buf);

	return list;
}


pwm* create_predefined_pwm(int PWM_NAME) {
	int i;
	pwm *p = malloc_cnt(sizeof(pwm));
	double **w;
	
	if(PWM_NAME == SF2_ASF) {
		p->length = 7;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SF2_ASF\0", 8);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-1.14; w[0][1]= 1.37; w[0][2]=-0.21; w[0][3]=-1.58;
		w[1][0]= 0.62; w[1][1]=-1.10; w[1][2]= 0.17; w[1][3]=-0.50;
		w[2][0]=-1.58; w[2][1]= 0.73; w[2][2]= 0.48; w[2][3]=-1.58;
		w[3][0]= 1.32; w[3][1]= 0.33; w[3][2]=-1.58; w[3][3]=-1.13;
		w[4][0]=-1.58; w[4][1]= 0.94; w[4][2]= 0.33; w[4][3]=-1.58;
		w[5][0]=-1.58; w[5][1]=-1.58; w[5][2]= 0.99; w[5][3]=-1.13;
		w[6][0]= 0.62; w[6][1]=-1.58; w[6][2]=-0.11; w[6][3]= 0.27;
		p->threshold = 1.956;
	}
	else if(PWM_NAME == SF2_ASF_LM) {
		p->length = 7;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SF2_ASF_LM\0", 11);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-1.58; w[0][1]= 1.55; w[0][2]=-1.35; w[0][3]=-1.55;
		w[1][0]= 0.15; w[1][1]=-0.53; w[1][2]= 0.44; w[1][3]=-0.28;
		w[2][0]=-0.97; w[2][1]= 0.79; w[2][2]= 0.41; w[2][3]=-1.28;
		w[3][0]= 0.74; w[3][1]= 0.33; w[3][2]=-0.98; w[3][3]=-0.92;
		w[4][0]=-1.19; w[4][1]= 0.72; w[4][2]= 0.51; w[4][3]=-1.09;
		w[5][0]=-0.75; w[5][1]=-0.62; w[5][2]= 1.03; w[5][3]=-0.52;
		w[6][0]= 0.43; w[6][1]=-0.99; w[6][2]= 0.00; w[6][3]= 0.20;
		p->threshold = 1.867;
	}
	else if(PWM_NAME == SC35) {
		p->length = 8;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SC35\0", 5);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-0.88; w[0][1]=-1.16; w[0][2]= 0.87; w[0][3]=-1.18;
		w[1][0]= 0.09; w[1][1]=-1.58; w[1][2]= 0.45; w[1][3]=-0.20;
		w[2][0]=-0.06; w[2][1]= 0.95; w[2][2]=-1.36; w[2][3]= 0.38;
		w[3][0]=-1.58; w[3][1]= 1.11; w[3][2]=-1.58; w[3][3]= 0.88;
		w[4][0]= 0.09; w[4][1]= 0.56; w[4][2]=-0.33; w[4][3]=-0.20;
		w[5][0]=-0.41; w[5][1]= 0.86; w[5][2]=-0.05; w[5][3]=-0.86;
		w[6][0]=-0.06; w[6][1]= 0.32; w[6][2]=-1.36; w[6][3]= 0.96;
		w[7][0]= 0.23; w[7][1]=-1.58; w[7][2]= 0.68; w[7][3]=-1.58;
		p->threshold = 2.383;
	}
	else if(PWM_NAME == SRP40) {
		p->length = 7;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SRP40\0", 6);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-0.13; w[0][1]= 0.56; w[0][2]=-1.58; w[0][3]= 0.92;
		w[1][0]=-1.58; w[1][1]= 0.68; w[1][2]=-0.14; w[1][3]= 0.37;
		w[2][0]= 1.28; w[2][1]=-1.12; w[2][2]=-1.33; w[2][3]= 0.23;
		w[3][0]=-0.33; w[3][1]= 1.24; w[3][2]=-0.48; w[3][3]=-1.14;
		w[4][0]= 0.97; w[4][1]=-0.77; w[4][2]=-1.58; w[4][3]= 0.72;
		w[5][0]=-0.13; w[5][1]= 0.13; w[5][2]= 0.44; w[5][3]=-1.58;
		w[6][0]=-1.58; w[6][1]=-0.05; w[6][2]= 0.80; w[6][3]=-1.58;
		p->threshold = 2.67;
	}
	else if(PWM_NAME == SRP55) {
		p->length = 6;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SRP55\0", 6);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-0.66; w[0][1]= 0.39; w[0][2]=-1.58; w[0][3]= 1.22;
		w[1][0]= 0.11; w[1][1]=-1.58; w[1][2]= 0.72; w[1][3]=-1.58;
		w[2][0]=-0.66; w[2][1]= 1.48; w[2][2]=-1.58; w[2][3]=-0.07;
		w[3][0]= 0.11; w[3][1]=-1.58; w[3][2]= 0.72; w[3][3]=-1.58;
		w[4][0]=-1.58; w[4][1]=-1.58; w[4][2]= 0.21; w[4][3]= 1.02;
		w[5][0]= 0.61; w[5][1]= 0.98; w[5][2]=-0.79; w[5][3]=-1.58;
		p->threshold = 2.676;
	}
	else if(PWM_NAME == SIRONI_1) {
		p->length = 8;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SIRONI_1\0", 9);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-0.2268; w[0][1]= 1.1095; w[0][2]=-0.7412; w[0][3]=-0.778 ;
		w[1][0]=-0.1553; w[1][1]=-0.5850; w[1][2]=-0.5059; w[1][3]= 0.6623;
		w[2][0]= 1.4676; w[2][1]=-1.5850; w[2][2]=-1.5850; w[2][3]=-1.5850;
		w[3][0]=-1.5850; w[3][1]=-1.5850; w[3][2]= 1.7227; w[3][3]=-1.5850;
		w[4][0]= 1.4676; w[4][1]=-1.5850; w[4][2]=-1.5850; w[4][3]=-1.5850;
		w[5][0]=-1.5850; w[5][1]=-1.5850; w[5][2]= 1.7227; w[5][3]=-1.5850;
		w[6][0]=-1.1089; w[6][1]=-0.3666; w[6][2]= 1.3850; w[6][3]=-1.2632;
		w[7][0]=-0.7516; w[7][1]=-0.3666; w[7][2]=-0.7412; w[7][3]= 0.9062;
		p->threshold = 2.379;
	}
	else if(PWM_NAME == SIRONI_2) {
		p->length = 7;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SIRONI_2\0", 9);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-0.5967; w[0][1]=-1.1920; w[0][2]= 0.1377; w[0][3]= 0.7395;
		w[1][0]=-1.5850; w[1][1]=-1.5850; w[1][2]= 1.7227; w[1][3]=-1.5850;
		w[2][0]= 0.2724; w[2][1]=-1.5850; w[2][2]=-1.5850; w[2][3]= 0.8693;
		w[3][0]=-1.5850; w[3][1]=-1.5850; w[3][2]= 1.3787; w[3][3]=-0.1786;
		w[4][0]=-1.5850; w[4][1]=-1.5850; w[4][2]= 1.7227; w[4][3]=-1.5850;
		w[5][0]=-1.5850; w[5][1]=-1.5850; w[5][2]= 1.7227; w[5][3]=-1.5850;
		w[6][0]=-0.7876; w[6][1]=-0.6295; w[6][2]= 1.3221; w[6][3]=-1.0276;
		p->threshold = 1.629;
	}
	else if(PWM_NAME == SIRONI_3) {
		p->length = 8;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "SIRONI_3\0", 9);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]=-0.8348; w[0][1]=-0.2771; w[0][2]=-0.5862; w[0][3]= 0.8602;
		w[1][0]=-1.0901; w[1][1]= 1.5403; w[1][2]=-1.2635; w[1][3]=-0.9174;
		w[2][0]=-1.5850; w[2][1]=-1.5850; w[2][2]=-1.5850; w[2][3]= 1.4143;
		w[3][0]=-1.5850; w[3][1]= 1.7779; w[3][2]=-1.5850; w[3][3]=-1.5850;
		w[4][0]=-1.5850; w[4][1]= 1.7779; w[4][2]=-1.5850; w[4][3]=-1.5850;
		w[5][0]=-1.5850; w[5][1]= 1.7779; w[5][2]=-1.5850; w[5][3]=-1.5850;
		w[6][0]= 1.1942; w[6][1]=-0.9801; w[6][2]=-1.4697; w[6][3]=-0.6474;
		w[7][0]= 0.5997; w[7][1]=-0.0441; w[7][2]=-0.1271; w[7][3]=-0.7494;
		p->threshold = 2.302;
	}
	else if(PWM_NAME == HSF_1) {
		p->length = 6;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "HSF_1\0", 6);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]= -11; w[0][1]=  -9; w[0][2]=  27; w[0][3]= 148;
		w[1][0]= 349; w[1][1]= -23; w[1][2]= -22; w[1][3]= -23;
		w[2][0]=   7; w[2][1]= -12; w[2][2]= 190; w[2][3]= -13;
		w[3][0]=  26; w[3][1]= -12; w[3][2]= 148; w[3][3]= -11;
		w[4][0]=  31; w[4][1]=  -3; w[4][2]= 112; w[4][3]= -10;
		w[5][0]=  33; w[5][1]=  11; w[5][2]=   2; w[5][3]=  23;
		p->threshold = 561.2;
	}
	else if(PWM_NAME == HSF_2) {
		p->length = 5;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "HSF_2\0", 6);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]= 261; w[0][1]= 15; w[0][2]= 15; w[0][3]= -16;
		w[1][0]= 261; w[1][1]= 15; w[1][2]= 15; w[1][3]= -16;
		w[2][0]=  -5; w[2][1]=  5; w[2][2]= 53; w[2][3]=  38;
		w[3][0]= 261; w[3][1]= 16; w[3][2]= 15; w[3][3]= -16;
		w[4][0]= 125; w[4][1]=  6; w[4][2]= 12; w[4][3]=  -9;
		p->threshold = 551.8;
	}
	else if(PWM_NAME == HSF_3) {
		p->length = 6;
		p->weights = (double**)malloc_cnt(sizeof(double*) * p->length);
		p->name = (char*)malloc_cnt(sizeof(char) * 12);
		strncpy(p->name, "HSF_3\0", 6);
		for(i = 0; i < p->length; i++)
			p->weights[i] = (double*)malloc_cnt(sizeof(double) * 4);
		w = p->weights;
		w[0][0]= -18; w[0][1]=   8; w[0][2]= 266; w[0][3]= -14;
		w[1][0]= 208; w[1][1]=  16; w[1][2]=  38; w[1][3]= -16;
		w[2][0]=  62; w[2][1]=  51; w[2][2]=   0; w[2][3]=   5;
		w[3][0]=  -2; w[3][1]=  18; w[3][2]= 268; w[3][3]=  -3;
		w[4][0]= 236; w[4][1]=  14; w[4][2]=  25; w[4][3]= -14;
		w[5][0]=  21; w[5][1]= 182; w[5][2]= -14; w[5][3]=   8;
		p->threshold = 706.4;
	}
	else {
		fprintf(stderr, "CANNOT FIND PWM NAME!\n");
		free_cnt(p);
		return NULL;
	}

	return p;
}

void free_pwm(pwm *p) {
	int i;
	for(i = 0; i < p->length; i++)
		free_cnt(p->weights[i]);
	free_cnt(p->weights);
	free_cnt(p->name);
	free_cnt(p);
}

void free_pwm_list(pwm_list* pl) {
	int i, j;
	for(i = 0; i < pl->length; i++) {
		for(j = 0; j < pl->pwms[i].length; j++)
			free_cnt(pl->pwms[i].weights[j]);
		free_cnt(pl->pwms[i].weights);
		free_cnt(pl->pwms[i].name);
	}
	free_cnt(pl->pwms);
	free_cnt(pl);
}

int check_motif(char* str, pwm* p) {
	int i, nn;
	double weight = 0.0; 

	if(strlen(str) < p->length)
		return 0;

	for(i = 0; i < p->length; i++) {
		if(str[i] == 'A') 	nn = 0;
		else if(str[i] == 'C') 	nn = 1;
		else if(str[i] == 'G') 	nn = 2;
		else if(str[i] == 'T') 	nn = 3;
		else			return 0;
		weight += p->weights[i][nn];
	}
	if(weight >= p->threshold)
		return 1;
	else
		return 0;
}

queue* find_motif(char* str, pwm* p) {
	int i, j, nn;
	double weight = 0.0;
	queue* q = init_queue();

	for(i = 0; i < strlen(str) - p->length + 1; i++) {
		weight = 0.0;
		for(j = 0; j < p->length; j++) {
			if(str[i + j] == 'A') 		nn = 0;
			else if(str[i + j] == 'C') 	nn = 1;
			else if(str[i + j] == 'G') 	nn = 2;
			else if(str[i + j] == 'T') 	nn = 3;
			else 				return NULL;
			weight += p->weights[j][nn];
		}
		if(weight >= p->threshold)
			enqueue(q, p, i);
	}
	
	return q;
}

queue* find_motif_list(char* str, pwm_list* pl) {
	int i, j, k, nn;
	double weight = 0.0;
	queue* q = init_queue();

	for(j = 0; j < pl->length; j++) {
		for(i = 0; i < strlen(str) - pl->pwms[j].length + 1; i++) {
			weight = 0.0;
			for(k = 0; k < pl->pwms[j].length; k++) {
				if(str[i + k] == 'A') 		nn = 0;
				else if(str[i + k] == 'C') 	nn = 1;
				else if(str[i + k] == 'G') 	nn = 2;
				else if(str[i + k] == 'T') 	nn = 3;
				else 				return NULL;
				weight += pl->pwms[j].weights[k][nn];
			}
			if(weight >= pl->pwms[j].threshold)
				enqueue(q, &(pl->pwms[j]), i);
		}
	}
	
	return q;
}



