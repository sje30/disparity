/****************************************************************************
***
***
*** convolve.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 17 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header$";
#endif



/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "convolve.h"
#include "dispnet.h"

/* - Defines - */

/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */

void double_convolve2d_wrap(Real *input, int inputWid, int inputHt,
			  Real *mask, int maskwid, int maskht, int maskextent,
			  Real *output)

{
  /* Input array is of size wid * ht 
   * Output array is of size wid * ht - same as the input array.
   * mask is of size (2 * maskextent + 1) * ( 2* maskextent + 1)
   */

  int	x,y;	/* Where the mask is currently centred */
  int	mx, my;	/* Mask x,y indexes */ 
  Real	sum;
  Real	*tmask;		/* Pointer to the mask - used in convolution */
  Real	*result;	/* Pointer to the output array. */
  int	startrow, startcol, offset;

  result = output;
  
  for(y=0; y < inputHt; y++) {
    for( x=0; x< inputWid; x++) {

      /* Find out the result of convolving the image array with the
       * mask, centred around the point x,y in the input image.  The
       * result is placed in the outptu array at the position x,y.
       */

      /* work out the starting row to start convolving */
      startrow = y -  maskextent;
      if (startrow < 0)
	startrow += inputHt;
      
      sum =0.0;
      tmask = mask;
      for(my=0; my < maskht; my++) {

	offset = inputWid*startrow;
	
	startcol = x - maskextent;
	if (startcol < 0)
	  startcol += inputWid;
	

	for(mx=0; mx < maskwid; mx++) {
	  printf( "mask %lf  * input %lf = %lf\n", *tmask,
		 input[offset+startcol],
		 (*tmask * input[startcol+offset]));
	  sum += (*tmask * input[offset+startcol]);
	  tmask++;
	  startcol++;
	  if (startcol == inputWid) {
	    startcol = 0; /* should just wrap round nicely. */
	  }
	}
	printf("\n");
	/* move onto next row of input to mask */
	startrow++;
	if ( startrow == inputHt) {
	  startrow = 0;
	}

      }

      /* Store the result of the convolution */
      *result = sum;
      result++;
    } /* next x */
  } /* next y */
}




void double_convolve1d_wrap(Real *input, int inputWid,
			  Real *mask, int maskwid, int maskextent,
			  Real *output)

{
  /* Input array is of size wid
   * Output array is of size wid  - same as the input array.
   * mask is 1d dimensional of size maskwid = (2*maskextent) + 1
   */

  int	x;	/* Where the mask is currently centred */
  int	mx;	/* Mask x index */ 
  Real	sum;
  Real	*tmask;	/* Pointer to the mask - used in convolution */
  Real  *result;	/* Pointer to the output array */
  int	startcol;

  result = output;
  for( x=0; x< inputWid; x++) {

    /* Find out the result of convolving the image array with the
     * mask, centred around the point x in the input image.  The
     * result is placed in the outptu array at the position x.
     */

    startcol = x - maskextent;
    if (startcol < 0)
      startcol += inputWid;
	
    tmask = mask;
    sum=0.0;
    for(mx=0; mx < maskwid; mx++) {
      printf( "mask %lf  * input %lf = %lf\n", *tmask, input[startcol],
	     (*tmask * input[startcol]));
      sum += (*tmask * input[startcol]);
      tmask++;
      startcol++;
      if (startcol == inputWid) {
	startcol = 0; /* should just wrap round nicely. */
      }
    }
    /* Store the result of the convolution */
    printf("next sum is %lf\n", sum);
    *result = sum;
    result++;
  } /* next x */

}


/*************************** Version Log ****************************/
/*
 * $Log$
 */

