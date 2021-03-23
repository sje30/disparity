/****************************************************************************
***
***
*** dispmasks.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 20 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef _DISPMASKS_H
#define _DISPMASKS_H


/* Function definitions */


void freeMasks();
void testMasks();
void createMasks();

void createDerivMask(Mask mask, Mask *derivMask);
void copyMask(Mask mask, Mask *rMask);
void writeMask(Mask mask, char *fname);

double half2lambda(double half);
void testMasks2();
void createMasks2();
#endif


/*************************** Version Log ****************************/
/* $Log$
 *
 */
