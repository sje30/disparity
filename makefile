#############################################################################
###
### makefile
### 
### Stephen Eglen
### COGS, University of Sussex.
###
### Created 09 Nov 95
###
### $Revision: 1.8 $
### $Date: 1995/12/14 21:42:36 $
#############################################################################
CC = gcc

# Place to look for includes of the form "file.h" (not <file.h>)
#CHEADERS = -I$(HOME)/Clib

#Place to look for libraries.
#CLIBDIRS = -L$(HOME)/Clib
#OPTS = -msupersparc -O3
OPTS = -O3
# Default rule for turning .c files into .o files uses the $CFLAGS 
# So, to get make to include my library, just add the -I and -L paths
# onto the CFLAGS
CFLAGS = -ansi -gstabs -g $(OPTS) $(CHEADERS) $(CLIBDIRS)


# By default the executable testnet will be placed into the file
# ~/bin/testnet

testnet : testnet.o readnet.o readnet.h dispscan.o dispinputs.o \
	 dispnet.o dispwts.o dispmasks.o  convolve.o testconvolve.o \
	disperrors.o cg_williams_module.o bp_check_deriv.o
	$(LINK.c) $(OPTS) -o $(HOME)/bin/testnet \
	testnet.o readnet.o dispscan.o dispwts.o \
	dispinputs.o dispnet.o dispmasks.o convolve.o testconvolve.o \
	disperrors.o cg_williams_module.o bp_check_deriv.o \
	-lm  $(CHEADERS) $(CLIBDIRS)


# testnet2 is same as testnet, it just allows me to create new
# executable while a program is still running.

testnet2 : testnet.o readnet.o readnet.h dispscan.o dispinputs.o \
	 dispnet.o dispwts.o dispmasks.o  convolve.o testconvolve.o \
	disperrors.o cg_williams_module.o bp_check_deriv.o
	$(LINK.c) $(OPTS) -o $(HOME)/bin/testnet2 \
	testnet.o readnet.o dispscan.o dispwts.o \
	dispinputs.o dispnet.o dispmasks.o convolve.o testconvolve.o \
	disperrors.o cg_williams_module.o bp_check_deriv.o \
	-lm  $(CHEADERS) $(CLIBDIRS)



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
SRCTAGS = readnet.c readnet.h testnet.c dispnet.h dispnet.c \
dispvars.h dispglobals.h dispwts.c rnd.h cg_williams_module.h \
dispinputs.c dispinputs.h dispmasks.c dispglobals.h convolve.c \
disperrors.c cg_williams_module.c bp_check_deriv.c

TAGS:   $(SRCTAGS)
	etags $(SRCTAGS)

tags.testconvolve: testconvolve.c convolve.c
	etags -f tags.testconvolve testconvolve.c convolve.c

############################ Version Log #############################
#
# $Log: makefile,v $
# Revision 1.8  1995/12/14 21:42:36  stephene
# updated tags rule
#
# Revision 1.7  1995/12/11  06:27:10  stephene
# testnet2 put in for checks.
# Removed -lmygen and CLIBDIRS, CHEADERDIRS for portability
#
# Revision 1.6  1995/12/07  15:43:05  stephene
# post nips sort out
#
# Revision 1.5  1995/11/23  16:20:26  stephene
# CG now installed
#
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
