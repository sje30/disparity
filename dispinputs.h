/****************************************************************************
***
*** Time-stamp: <12 Nov 95 23:07:31 stephene>
***
*** dispinputs.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 12 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/12 23:07:36 $
****************************************************************************/


#ifndef _DISPINPUTS_H
#define _DISPINPUTS_H

#include "dispnet.h"

void getInputVector();
void setBiases();
void clearActivationArray();
void createInputVectorsAndShifts();
void testInputVectors();
void extractInputVector(Array source, Array dest,
			int colstart, int rowstart);
void readInputFile(char *inputFile, Array arr);

void getNextInputVector(int vecnum);
#endif


/*************************** Version Log ****************************/
/* $Log: dispinputs.h,v $
 * Revision 1.1  1995/11/12  23:07:36  stephene
 * Initial revision
 *
 *
 */
