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
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <sys/resource.h>
#include <sys/time.h>
#include "utils.h"

#include "ksort.h"
#define pair64_lt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y < (b).y))
KSORT_INIT(128, pair64_t, pair64_lt)
KSORT_INIT(64,  uint64_t, ks_lt_generic)

#include "kseq.h"
KSEQ_INIT2(, gzFile, err_gzread)

int alloc_cnt = 0;

/********************
 * System utilities *
 ********************/

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) err_fatal(func, "Out of memory");
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s\n", func, msg);
	exit(EXIT_FAILURE);
}

void _err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb) 
		_err_fatal_simple("fwrite", strerror(errno));
	return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

int err_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf", strerror(saveErrno));
	return done;
}

int err_fputc(int c, FILE *stream)
{
	int ret = putc(c, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputc", strerror(errno));
	}

	return ret;
}

int err_fputs(const char *s, FILE *stream)
{
	int ret = fputs(s, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputs", strerror(errno));
	}

	return ret;
}

int err_puts(const char *s)
{
	int ret = puts(s);
	if (EOF == ret)
	{
		_err_fatal_simple("puts", strerror(errno));
	}

	return ret;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) _err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
	/* Calling fflush() ensures that all the data has made it to the
	   kernel buffers, but this may not be sufficient for remote filesystems
	   (e.g. NFS, lustre) as an error may still occur while the kernel
	   is copying the buffered data to the file server.  To be sure of
	   catching these errors, we need to call fsync() on the file
	   descriptor, but only if it is a regular file.  */
	{
		struct stat sbuf;
		if (0 != fstat(fileno(stream), &sbuf))
			_err_fatal_simple("fstat", strerror(errno));
		
		if (S_ISREG(sbuf.st_mode))
		{
			if (0 != fsync(fileno(stream)))
				_err_fatal_simple("fsync", strerror(errno));
		}
	}
#endif
    return ret;
}

int err_fclose(FILE *stream) 
{
	int ret = fclose(stream);
	if (ret != 0) _err_fatal_simple("fclose", strerror(errno));
	return ret;
}

int err_gzclose(gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		_err_fatal_simple("gzclose", Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}

	return ret;
}

/*********
 * Timer *
 *********/

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}



/*** Utility functions written by Yunku yeu***/

// Make the input str to upper case string
void make_str_hc(char* str) {
	int i;
	for(i = 0; i < strlen(str); i++)
		if((str[i] >= 'a') && (str[i] <= 'z'))
			str[i] -= 'a'-'A';
}

// Make the input str to lower case string
void make_str_lc(char* str) {
	int i;
	for(i = 0; i < strlen(str); i++)
		if((str[i] >= 'A') && (str[i] <= 'Z'))
			str[i] += 'a'-'A';
}

// Convert DNA sequence (from) into its reverse complement and store it into (to)
// Parameters ***************************************************************
//   char* to: 		Reverse complemented string. It should be alloced with (len + 1) size
//   char* from:	Original string
//   int len:		length of original string
// **************************************************************************
// Return: void
// **************************************************************************
void convert_reverse_complement(char* to, const char* from, int len) {
	int i;
	for(i = 0; i < len; i++) {
		if     (from[i] == 'A')	to[len - 1 - i] = 'T';
		else if(from[i] == 'a')	to[len - 1 - i] = 't';
		else if(from[i] == 'T')	to[len - 1 - i] = 'A';
		else if(from[i] == 't')	to[len - 1 - i] = 'a';
		else if(from[i] == 'G')	to[len - 1 - i] = 'C';
		else if(from[i] == 'g')	to[len - 1 - i] = 'c';
		else if(from[i] == 'C')	to[len - 1 - i] = 'G';
		else if(from[i] == 'c')	to[len - 1 - i] = 'g';
		else 			to[len - 1 - i] = from[i];
	}	

	to[len] = '\0';
	return;
}

// Alias of convert_reverse_complement. to match with strncpy()
char* strncpy_rc(char* dest, const char* src, size_t n) {
	convert_reverse_complement(dest, src, n);
	return dest;
}


// Wrapper function of convert_reverse_complement(), string length is given.
// It make an RC of the original sequence, then free the original, return the RC.
// Parameters ***************************************************************
//   char* original: 	Original string to be converted into RC
//   int size:		Length of the original string
// **************************************************************************
// Return: void
// **************************************************************************
char* create_rc_free_exist_len(char* original, int size) {
	char* rc = (char*)malloc_cnt(sizeof(char) * (size + 1));

	convert_reverse_complement(rc, original, size);

	free_cnt(original);
	return rc;
}


// create_rc_free_exist_len() without length
// [WARNING] This function is valid only when original sequence has no empty space.
// Otherwise, KERNEL ERROR might be occured because memory size of input string can be changed.
char* create_rc_free_exist(char* original) {
	return create_rc_free_exist_len(original, strlen(original));
}

// Wrapper of convert_reverse_complement(). the original sequence is not freed.
char* create_rc_without_free(char* original) {
	char* rc = (char*)calloc(sizeof(char), strlen(original) + 1);

	convert_reverse_complement(rc, original, strlen(original));

	return rc;
}


// strcmp function for DNA sequence, ignoring upper / lower character
int seqncmp_ignore_case(const char* seq1, const char* seq2, int cmplen) {
	int i;
	char c1, c2;
	for(i = 0; i < cmplen; i++) {
		c1 = ((seq1[i] >= 'a') && (seq1[i] <= 'z')) ? seq1[i] + 'A' - 'a' : seq1[i];
		c2 = ((seq2[i] >= 'a') && (seq2[i] <= 'z')) ? seq2[i] + 'A' - 'a' : seq2[i];

		if(c1 < c2)
			return -1;
		else if(c1 > c2)
			return 1;
	}

	return 0;
}

// return nonzero if there is any overlap between two range, [s1, e1] and [s2, e2]
// s1 <= e1 and s2 <= e2 are required.
int if_any_overlap(int s1, int e1, int s2, int e2) {
	if(s1 > e1) {
		int temp = e1;
		e1 = s1;
		s1 = temp;
	}
	if(s2 > e2) {
		int temp = e2;
		e2 = s2;
		s2 = temp;
	}

	if((e1 < s2) || (e2 < s1))
		return 0;
	else
		return 1;
}

#ifdef DEBUG

void print_cnt() {
	printf("current alloc_cnt = %d\n", alloc_cnt);
}
void free_cnt(void* x) {
	free(x);
	alloc_cnt--;
}

void* malloc_cnt(int size) {
	void* out = malloc(size);
	if(out != NULL)
		alloc_cnt++;
	return out;
}

#endif

