/****************************************************************************
***
*** Time-stamp: <>
***
*** dispwts.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 12 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef _DISPWTS_H
#define _DISPWTS_H

void createWeights(int len);
void freeWeights();
Real *nextFreeWeight();
void noMoreWeights();		
void initWtsRnd();

#endif


/*************************** Version Log ****************************/
/* $Log$
 *
 */
