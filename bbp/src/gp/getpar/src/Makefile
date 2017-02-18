# Copyright 1990 Science Applications International Corporation
#
# Makefile for libget library.
#		cloned from Makefile.lib
# SccsId:  @(#)Makefile	56.1	10/25/93	Copyright 1990 Science Applications International Corporation
#
 
CSRCS		= defpar.c \
		  getpar.c \
		  getarg.c \
		  getlocation.c

COBJS		= defpar.o \
		  getpar.o \
		  getarg.o \
		  getlocation.o

FSRCS		= getpar.c \
		  getarg.c

FOBJS		= fgetpar.o \
		  fgetarg.o

PUBLIC_HDRS	= libget.h
PRIVATE_HDRS	=

CHDRS           = $(PUBLIC_HDRS) $(PRIVATE_HDRS)
FHDRS		=
HDRS		= $(CHDRS) $(FHDRS)

MAKEFILE        = Makefile

ALLSRCS         = $(CSRCS) $(MAKEFILE) $(HDRS) $(CHANGES) $(DOC)

LIBS		= 
STDLIBS		= 
LIBRARY		= libget.a
#FLIBRARY	= libfpar.a
MAN1		=
MAN3		= getpar.3 getarg.3
CHANGES		= Changes.libget
DOC		= $(MAN1) $(MAN3) $(CHANGES)
SUBDIRS		= 
INCLUDES	= 
CFLAGS		= $O $(INCLUDES) -DENVIRONMENT
FFLAGS		= $O
LDFLAGS		= $O
LFLAGS		= -bchu $(INCLUDES)

DESTDIR		= ../../lib
INCLDIR		= ../../include
MAN1DIR		= ../../doc/man/man1
MAN3DIR		= ../../doc/man/man3
CHANGESDIR	= ../../doc/Changes

#CC		= cc
CC		= gcc
FC		= f77
O		= -xO3 -xdepend -xlibmil -fsimple
O		= -O3
#O		= -g
LD		= $(CC)
LINT		= lint
YACC		= yacc
LEX		= lex

LMODE		= 664
PMODE		= 775
GROUP		= scm

GET		= sccs get
GFLAGS		= -s -r56.1
SCCSGETFLAGS	= -s -r56.1
SCLEAN		= sccs clean

RM		= rm
INSTALL		= install
MAKE		= make
PRINT		= lpr -p
AR		= ar
RANLIB		= ranlib
MKDIR		= mkdir -p

FLAGS		= $(MFLAGS) "CC=$(CC)" "FC=$(FC)" "O=$O" \
		  "GET=$(GET)" "GFLAGS=$(GFLAGS)" "SCLEAN=$(SCLEAN)" \
		  "LMODE=$(LMODE)" "PMODE=$(PMODE)" "GROUP=$(GROUP)" \
		  "PRINT=$(PRINT)" "MAKE=$(MAKE)" "LDFLAGS=$(LDFLAGS)" \
		  "PCC=$(PCC)"

all:		$(LIBRARY)  $(MAN3)

$(COBJS) $(FOBJS): libget.h

fgetpar.o:	getpar.c
		cp getpar.c fgetpar.c
		$(CC) $(CFLAGS) -DFORTRAN -c fgetpar.c
		$(RM) -f fgetpar.c

fgetarg.o:	getarg.c
		cp getarg.c fgetarg.c
		$(CC) $(CFLAGS) -DFORTRAN -c fgetarg.c
		$(RM) -f fgetarg.c

$(LIBRARY):	$(CHDRS) $(COBJS) $(FOBJS)
		-$(RM) -f $(LIBRARY)
		$(AR) cq $(LIBRARY) $(COBJS) $(FOBJS)
		if [ -x /usr/bin/ranlib -o -x /bin/ranlib ]; then \
			$(RANLIB) $(LIBRARY); \
		fi
		cp libget.a ../lib/

#$(FLIBRARY):	$(FHDRS) $(FOBJS)
#		-$(RM) -f $(FLIBRARY)
#		$(AR) cq $(FLIBRARY) $(FOBJS)
#		$(RANLIB) $(FLIBRARY)

includes:	$(HDRS)
		$(MKDIR) $(INCLDIR)
		for i in $(PUBLIC_HDRS); do \
			$(INSTALL) -m $(LMODE) -g $(GROUP) $$i $(INCLDIR);\
		done

install:	all $(DOC)
		$(MKDIR) $(DESTDIR)
		$(INSTALL) -m $(LMODE) -g $(GROUP) $(LIBRARY) \
		$(DESTDIR)
		if [ -x /usr/bin/ranlib -o -x /bin/ranlib ]; then \
			$(RANLIB) $(DESTDIR)/$(LIBRARY); \
		fi
#			$(RANLIB) $(DESTDIR)/$(FLIBRARY); \
#			Uncomment below as necessary.
#		$(MKDIR) $(MAN1DIR)
#		for i in $(MAN1); do \
#			$(INSTALL) -m $(LMODE) -g $(GROUP) $$i $(MAN1DIR);\
#		done
		$(MKDIR) $(MAN3DIR)
		for i in $(MAN3); do \
			$(INSTALL) -m $(LMODE) -g $(GROUP) $$i $(MAN3DIR);\
		done
#		$(MKDIR) $(CHANGESDIR)
#		for i in $(CHANGES); do \
#			$(INSTALL) -m $(LMODE) -g $(GROUP) $$i $(CHANGESDIR);\
#		done

clean:
		-$(RM) -f $(COBJS) $(FOBJS) libget.a core

xclean:
		-$(RM) -f $(COBJS) $(FOBJS) libget.a core
		-$(RM) -f $(PROGRAM) $(LIBRARY) 
		-$(SCLEAN)
		
lint:		$(CSRCS) $(HDRS) $(FSRCS)
		$(LINT) $(LFLAGS) $(CSRCS) $(FSRCS)

print:		$(CSRCS) $(HDRS) Makefile
		$(PRINT) $(HDRS) Makefile

sccs:
		$(GET) $(GFLAGS) $(CSRCS) $(HDRS) $(DOC) Makefile

