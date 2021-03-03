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
*** $Revision: 1.3 $
*** $Date: 1995/12/09 21:11:18 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/convolve.c,v 1.3 1995/12/09 21:11:18 stephene Exp $";
#endif



/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "convolve.h"
#include "dispnet.h"

/* - Defines - */
#define MASK1D 1


/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */

void double_convolve_wrap(Array input,
			  Mask mask,
			  Array output)
{

  /* Call the appropriate convolution code, according to whether the
   * array is one dimensional or two dimensional.  As long as the
   * dimensions are correct, the input is convolved with the mask to
   * produce the output.  All arrays must be allocated before calling
   * this routine.
   *
   * If the mask is bigger than the input, then this routine will
   * produce an error.  */
  
  /*** Local Variables ***/
  int	inputWid, inputHt;
  Real	*inputData;
  int	maskWid, maskHt;
  Real	*maskData;
  int	maskExtent;
  
  int	outputWid, outputHt;
  Real	*outputData;

  /* Get the sizes of the inputs and the masks. */
  inputWid = input.wid; inputHt = input.ht; inputData = input.data;
  maskWid = mask.wid; maskHt = mask.ht; maskData = mask.data;
  maskExtent = mask.maskExtent;
  outputData = output.data;
  
  if ( maskHt == MASK1D) {
    if (inputHt == 1) {
      /* One dimensional convolution */
      if ( maskWid > inputWid ) {
	printf("%s: maskWid (%d) must not be greater than inputWid (%d)\n",
	       __FUNCTION__, maskWid, inputWid);
	exit(-1);
      }
      double_convolve1d_wrap( inputData, inputWid,
			      maskData, maskWid, maskExtent,
			      outputData);
    } else {
      printf("Error: inputData is one dimensional, but mask is not\n");
      exit(-1);
    }
  } else {
    /* 2d convolution */
    if ( maskWid > inputWid ) {
      printf("%s: maskWid (%d) must not be greater than inputWid (%d)\n",
	     __FUNCTION__, maskWid, inputWid);
      exit(-1);
    }
    if ( maskHt > inputHt ) {
      printf("%s: maskHt (%d) must not be greater than inputHt (%d)\n",
	     __FUNCTION__, maskHt, inputHt);
      exit(-1);
    }
		
    double_convolve2d_wrap( inputData, inputWid, inputHt,
			    maskData, maskWid, maskHt,  maskExtent,
			    outputData);
  }
}

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
	  /*
	  printf( "mask %lf  * input %lf = %lf\n", *tmask,
		 input[offset+startcol],
		 (*tmask * input[startcol+offset]));
		 */
	  sum += (*tmask * input[offset+startcol]);
	  tmask++;
	  startcol++;
	  if (startcol == inputWid) {
	    startcol = 0; /* should just wrap round nicely. */
	  }
	}
/* 	printf("\n"); */
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
      /*
       printf( "mask %lf  * input %lf = %lf\n", *tmask, input[startcol],
	     (*tmask * input[startcol])); */
      sum += (*tmask * input[startcol]);
      tmask++;
      startcol++;
      if (startcol == inputWid) {
	startcol = 0; /* should just wrap round nicely. */
      }
    }
    /* Store the result of the convolution */
/*     printf("next sum is %lf\n", sum); */
    *result = sum;
    result++;
  } /* next x */

}


/*************************** Version Log ****************************/
/*
 * $Log: convolve.c,v $
 * Revision 1.3  1995/12/09  21:11:18  stephene
 * commented out the printf statements in convolve 2d.
 *
 * Revision 1.2  1995/11/21  23:31:37  stephene
 * About to include CG Code
 *
 * Revision 1.1  1995/11/17  21:03:55  stephene
 * Initial revision
 *
 */

