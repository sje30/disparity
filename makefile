#############################################################################
###
### makefile
### 
### Stephen Eglen
### COGS, University of Sussex.
###
### Created 09 Nov 95
###
### $Revision: 1.4 $
### $Date: 1995/11/21 23:33:29 $
#############################################################################
CC = gcc

# Place to look for includes of the form "file.h" (not <file.h>)
CHEADERS = -I$(HOME)/Clib

#Place to look for libraries.
CLIBDIRS = -L$(HOME)/Clib

# Default rule for turning .c files into .o files uses the $CFLAGS 
# So, to get make to include my library, just add the -I and -L paths
# onto the CFLAGS
CFLAGS = -ansi -gstabs -g $(CHEADERS) $(CLIBDIRS)

testnet : testnet.o readnet.o readnet.h dispscan.o dispinputs.o \
	 dispnet.o dispwts.o dispmasks.o  convolve.o testconvolve.o \
	disperrors.o cg_williams_module.o bp_check_deriv.o
	$(LINK.c) -o $(HOME)/bin/testnet \
	testnet.o readnet.o dispscan.o dispwts.o \
	dispinputs.o dispnet.o dispmasks.o convolve.o testconvolve.o \
	disperrors.o cg_williams_module.o bp_check_deriv.o \
	-lm -lmygen $(CHEADERS) $(CLIBDIRS)

testsig: testsig.o
	$(LINK.c) -o testsig testsig.o


testconvolve: testconvolve.o convolve.o
	$(LINK.c) -o testconvolve testconvolve.o convolve.o  \
	-lm -lmygen

nagconvolve: nagconvolve.c
	cc $(CHEADERS) $(CLIBDIRS) -o nagconvolve \
	nagconvolve.c \
	-Wl,-z,nodefs /usr/local/lib/libnag.so \
	-lM77 -lF77 -lsunmath  -lm


# let make automatically build .o file from .l file
dispscan.c : dispscan.l dispvars.h 
	lex -t dispscan.l > dispscan.c

dispscan.o : dispscan.c

clean:
	rm -f *.o *~


# Use my version of tags.
SRCTAGS = readnet.c readnet.h testnet.c dispnet.h dispnet.c dispwts.c \
dispinputs.c dispinputs.h dispmasks.c dispglobals.h convolve.c \
disperrors.c cg_williams_module.c bp_check_deriv.c

TAGS:   $(SRCTAGS)
	etags $(SRCTAGS)

tags.testconvolve: testconvolve.c convolve.c
	etags -f tags.testconvolve testconvolve.c convolve.c

############################ Version Log #############################
#
# $Log: makefile,v $
# Revision 1.4  1995/11/21  23:33:29  stephene
# About to include CG Code
#
# Revision 1.3  1995/11/17  00:06:34  stephene
# # Daily update
#
# Revision 1.2  1995/11/10  22:02:09  stephene
# *** empty log message ***
#
# Revision 1.1  1995/11/09  20:49:41  stephene
# Initial revision
#
#
