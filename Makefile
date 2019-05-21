SHELL         = /bin/sh
PROG          = CRISPinatoR
CC            = gcc
CFLAGS=		-g -Wall -Wno-unused-function -O2
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
LOBJS=		utils.o kthread.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o bwamem_extra.o malloc_wrap.o \
			QSufSort.o bwt_gen.o rope.o rle.o is.o bwtindex.o 
#AOBJS=		bwashm.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o \
#			bwape.o kopen.o pemerge.o maxk.o \
#			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
#			bwtsw2_chain.o fastmap.o bwtsw2_pair.o
AOBJS=		genomic_entry.o	crispinator.o bwtaln.o bwtgap.o bntseq.o bwase.o motif.o tools.o
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

CRISPinatoR:libbwa.a $(AOBJS) main.o 
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -lbwa $(LIBS)

bwamem-lite:libbwa.a example.o
		$(CC) $(CFLAGS) $(DFLAGS) example.o -o $@ -L. -lbwa $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )


main.o: utils.h crispinator.h main.c
crispinator.o: utils.h genomic_entry.h crispinator.h motif.h
genomic_entry.o: utils.h genomic_entry.c bwt.h kseq.h bwtaln.h
motif.o: utils.h motif.c
utils.o: utils.c utils.h ksort.h malloc_wrap.h kseq.h
tools.o: utils.h genomic_entry.h tools.h motif.h

bntseq.o: bntseq.h utils.h kseq.h malloc_wrap.h khash.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kvec.h
bwa.o: kseq.h
bwamem.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.o: ksort.h utils.h kbtree.h
bwamem_extra.o: bwa.h bntseq.h bwt.h bwamem.h kstring.h malloc_wrap.h
bwamem_pair.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h kvec.h
bwamem_pair.o: utils.h ksw.h
bwase.o: bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h malloc_wrap.h
bwase.o: bwa.h ksw.h
bwt.o: utils.h bwt.h kvec.h malloc_wrap.h
bwt_gen.o: QSufSort.h malloc_wrap.h
bwtaln.o: bwtaln.h bwt.h bwtgap.h utils.h bwa.h bntseq.h malloc_wrap.h
bwtindex.o: bntseq.h bwa.h bwt.h utils.h rle.h rope.h malloc_wrap.h
bwtgap.o: bwtgap.h bwt.h bwtaln.h malloc_wrap.h
is.o: malloc_wrap.h
kstring.o: kstring.h malloc_wrap.h
ksw.o: ksw.h malloc_wrap.h
malloc_wrap.o: malloc_wrap.h
QSufSort.o: QSufSort.h
rle.o: rle.h
rope.o: rle.h rope.h
