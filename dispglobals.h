/****************************************************************************
***
***
*** dispglobals.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 20 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/12/14 21:19:50 $
****************************************************************************/

/*
 * Global Variables used by the disparity program that are not
 * parameters of the system.
 */


#ifndef _DISPGLOBALS_H
#define _DISPGLOBALS_H

#include "dispnet.h"

extern Mask uMask, uMaskD;	/* Short Range Mask and its version for derivatives. */
extern Mask vMask, vMaskD;	/* Long Range Mask and derivative version. */

extern Array shifts; /* Stores the shifts array. */
extern Array inputs; /* Store the input vectors to the network */

extern allActns_t allActns;	/* Activation info stored for each input
			 * vector
			 */

/* 
 * z is the array of output values, arranged in either a 1d or 2d set up.
 * ztilde is the corresponding array of short range values, and 
 * zbar is the long range average values.
 *
 * ztildeminz, zbarminz are useful auxiliary variables.
 * ztildeminz[i] = ztilde[i] - z[i]. (useful for computing wt change vec)
 */	

extern Array z, zbar, ztilde;
extern Array ztildeminz, zbarminz;


/*
 * See maths section Computing dU/dx_a
 *
 *dudx[a] = du/dx_a =  partial diff of U with respect to
 *		       activation of output unit x 
 *dvdx[a] = dv/dx_a =  partial diff of V ...
 *
 * da[i]  = 1/V dv/dx_a - 1/U du/dx_a = \delta_a
 *
 * Hence da is the array of delta's for the top layer of cells.
 *
 */

extern Array		dudx, dvdx, da;


/* dw[i] = weight change for weight i, over all input vectors.
 * onedw[i] = weight change for weight i, given just one input vector.
 */
extern Array dw; /* Stores the partials for each weight, summed over all inputs */
extern Array onedw; /* Partials for just one input vector. */



extern FILE	*opfp; /* Diagnostic file pointer for outputting various bits of info */

#endif


/*************************** Version Log ****************************/
/* $Log: dispglobals.h,v $
 * Revision 1.1  1995/12/14  21:19:50  stephene
 * Initial revision
 *
 *
 */
