/****************************************************************************
***
*** Time-stamp: <12 Nov 95 22:16:44 stephene>
***
*** dispwts.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 12 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/12 22:16:31 $
****************************************************************************/


#ifndef _DISPWTS_H
#define _DISPWTS_H

void createWeights(int len);
void freeWeights();
Real *nextFreeWeight();
void noMoreWeights();		
void initWtsRnd();
void writeWts(char *fname);
#endif


/*************************** Version Log ****************************/
/* $Log: dispwts.h,v $
 * Revision 1.1  1995/11/12  22:16:31  stephene
 * Initial revision
 *
 *
 */
