/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef LH3_UTILS_H
#define LH3_UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#ifdef __GNUC__
// Tell GCC to validate printf format string and args
#define ATTRIBUTE(list) __attribute__ (list)
#else
#define ATTRIBUTE(list)
#endif

#define err_fatal_simple(msg) _err_fatal_simple(__func__, msg)
#define err_fatal_simple_core(msg) _err_fatal_simple_core(__func__, msg)

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)

#define xassert(cond, msg) if ((cond) == 0) _err_fatal_simple_core(__func__, msg)

typedef struct {
	uint64_t x, y;
} pair64_t;

typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; pair64_t *a; } pair64_v;

#ifdef __cplusplus
extern "C" {
#endif

	void err_fatal(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void err_fatal_core(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void _err_fatal_simple(const char *func, const char *msg) ATTRIBUTE((noreturn));
	void _err_fatal_simple_core(const char *func, const char *msg) ATTRIBUTE((noreturn));
	FILE *err_xopen_core(const char *func, const char *fn, const char *mode);
	FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp);
	gzFile err_xzopen_core(const char *func, const char *fn, const char *mode);
    size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
	size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream);

	int err_gzread(gzFile file, void *ptr, unsigned int len);
	int err_fseek(FILE *stream, long offset, int whence);
#define err_rewind(FP) err_fseek((FP), 0, SEEK_SET)
	long err_ftell(FILE *stream);
	int err_fprintf(FILE *stream, const char *format, ...)
        ATTRIBUTE((format(printf, 2, 3)));
	int err_printf(const char *format, ...)
        ATTRIBUTE((format(printf, 1, 2)));
	int err_fputc(int c, FILE *stream);
#define err_putchar(C) err_fputc((C), stdout)
	int err_fputs(const char *s, FILE *stream);
	int err_puts(const char *s);
	int err_fflush(FILE *stream);
	int err_fclose(FILE *stream);
	int err_gzclose(gzFile file);

	double cputime();
	double realtime();

	void ks_introsort_64 (size_t n, uint64_t *a);
	void ks_introsort_128(size_t n, pair64_t *a);

#ifdef __cplusplus
}
#endif

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

// added by yyk
#define BUF_LEN 		40960
#define NAME_LEN		1024
void make_str_hc(char* str);
void make_str_lc(char* str);
int seqncmp_ignore_case(const char*, const char*, int);
char* create_rc_free_exist(char*);
char* create_rc_free_exist_len(char*, int);
char* create_rc_without_free(char*);
void convert_reverse_complement(char* to, const char* from, int size);
int if_any_overlap(int s1, int e1, int s2, int e2);
char* strncpy_rc(char* dest, const char* src, size_t n);

#define FOR_LENTRY(cursor, head)	for((cursor) = (head); (cursor) != NULL; (cursor) = (cursor)->next)
#define IF_OVERLAP(start1, end1, start2, end2)	if((((start1) <= (start2)) && ((start2) <= (end1))) || (((start1) <= (end2)) && ((end2) <= (end1))) || (((start2) <= (start1)) && ((start1) <= (end2))) || (((start2) <= (end1)) && ((end1) <= (end2))))

//#define DEBUG

#ifdef DEBUG

void free_cnt(void* x);
void* malloc_cnt(int size);
void print_cnt();

#else
#define free_cnt(x) 			free(x)
#define malloc_cnt(size) 		malloc(size)

#endif


#endif
