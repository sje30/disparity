/****************************************************************************
***
***
*** testnet.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 22 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/22 00:54:50 $
****************************************************************************/


#ifndef _TESTNET_H
#define _TESTNET_H


double	evalFn(Real *wts);
void	evalPartials(Real *wts, Real *derivs);
char	finishedFn(Real *wts, int iteration);


void calcMeritAndPartials();
void vecNegate(Real *vec, int len);
void copyVec(Real *wdest, Real *wsrc, int numWeights);
#endif


/*************************** Version Log ****************************/
/* $Log: testnet.h,v $
 * Revision 1.1  1995/11/22  00:54:50  stephene
 * Initial revision
 *
 *
 */
