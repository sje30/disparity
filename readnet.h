/****************************************************************************
***
*** Time-stamp: <>
***
*** readnet.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/09 20:48:50 $
****************************************************************************/


#ifndef _READNET_H
#define _READNET_H

/* Function definitions for readnet.c
 * This file can be included by all other source files. */



void freeActivationsArray();
void readNet(char *fname);
int cellArea(int tlx, int tly, int brx, int bry);
void createWeights(int len);
void freeWeights();
Real *nextFreeWeight();



#endif


/*************************** Version Log ****************************/
/* $Log: readnet.h,v $
 * Revision 1.1  1995/11/09  20:48:50  stephene
 * Initial revision
 *
 *
 */
