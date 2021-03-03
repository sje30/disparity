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
*** $Revision: 1.3 $
*** $Date: 1995/11/17 00:07:24 $
****************************************************************************/


#ifndef _READNET_H
#define _READNET_H

/* Function definitions for readnet.c
 * This file can be included by all other source files. */



void	freeActivationsArray();
void	readNet(char *fname);
int	cellArea(int tlx, int tly, int brx, int bry);
void 	freeCellInfo();
void	freePreCellInfo();
void	createPreCellInfo();
void	printPreCellInfo();

#endif


/*************************** Version Log ****************************/
/* $Log: readnet.h,v $
 * Revision 1.3  1995/11/17  00:07:24  stephene
 * # daily update
 *
 * Revision 1.2  1995/11/10  22:04:41  stephene
 * about to make a snapshot
 *
 * Revision 1.1  1995/11/09  20:48:50  stephene
 * Initial revision
 *
 *
 */
