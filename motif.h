#ifndef _MOTIF_
#define _MOTIF_
// struct for ESE & ESS list

#define SF2_ASF 	1
#define SF2_ASF_LM	2
#define SC35		3
#define SRP40		4
#define SRP55		5
#define HSF_1		6
#define HSF_2		7
#define HSF_3		8
#define SIRONI_1	9
#define SIRONI_2	10
#define SIRONI_3	11


typedef struct kmers kmers;
struct kmers {
	int length;
	char** mers;
	int* k_size;
	double* weights;
	char** notes;
};

typedef struct pwm pwm;
struct pwm {
	int length;
	double** weights;
	double threshold;
	char* name;
};

typedef struct pwm_list pwm_list;
struct pwm_list {
	int length;
	pwm* pwms;
};
	

typedef struct node node;
struct node {
	int value;
	pwm* motif;
	node *next;
};

typedef struct queue queue;
struct queue {
	node *first, *last;
	int length;
};


kmers* load_kmers(char*);
void free_kmers(kmers*);
void compare_print_mers_list(FILE*, kmers*, char*, char*, int);
pwm* create_pwm(int);
void free_pwm(pwm*);
void free_pwm_list(pwm_list*);
int check_motif(char*, pwm*);
queue* find_motif(char*, pwm*);
queue* init_queue(void);
void free_queue(queue*);
void enqueue(queue*, pwm*, int);
int dequeue(queue*);
node* dequeue_nofree(queue*);
pwm_list* create_pwm_from_file(const char* path);
queue* find_motif_list(char* str, pwm_list* pl);

#endif
