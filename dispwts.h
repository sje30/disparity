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
*** $Revision: 1.4 $
*** $Date: 1995/12/11 06:26:33 $
****************************************************************************/


#ifndef _DISPWTS_H
#define _DISPWTS_H

void createWeights(int len);
void freeWeights();
Real *nextFreeWeight(int preCellNum, int postCellNum);
void noMoreWeights();		
void initWtsRnd();
void writeWts(char *fname);
void readWts(char *fname);
void printWtsInfo();
#endif


/*************************** Version Log ****************************/
/* $Log: dispwts.h,v $
 * Revision 1.4  1995/12/11  06:26:33  stephene
 * new function readWts created
 *
 * Revision 1.3  1995/11/21  23:33:20  stephene
 * About to include CG Code
 *
 * Revision 1.2  1995/11/17  00:06:21  stephene
 * Daily update
 *
 * Revision 1.1  1995/11/12  22:16:31  stephene
 * Initial revision
 *
 *
 */
