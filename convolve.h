/****************************************************************************
***
***
*** convolve.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 17 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef _CONVOLVE_H
#define _CONVOLVE_H

#include "dispnet.h"

/** Main calling function **/
void double_convolve_wrap(Array input,
			  Mask mask,
			  Array output);


void double_convolve1d_wrap(Real *input, int inputWid,
			    Real *mask, int maskwid, int maskextent,
			    Real *output);

void double_convolve2d_wrap(Real *input, int inputWid, int inputHt,
			    Real *mask, int maskwid, int maskht, int maskextent,
			    Real *output);


#endif


/*************************** Version Log ****************************/
/* $Log$
 *
 */
