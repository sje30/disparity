#############################################################################
###
### makefile
### 
### Stephen Eglen
### COGS, University of Sussex.
###
### Created 09 Nov 95
###
### $Revision$
### $Date$
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

testnet : testnet.o readnet.o readnet.h
	$(LINK.c) -o testnet testnet.o readnet.o -lm  $(CHEADERS) $(CLIBDIRS)


clean:
	rm -f *.o *~


# Use my version of tags.
TAGS:   readnet.c readnet.h testnet.c
	etags readnet.c readnet.h testnet.c



############################ Version Log #############################
#
# $Log$
#
