/****************************************************************************
***
***
*** dispmasks.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 20 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/20 14:43:28 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispmasks.c,v 1.1 1995/11/20 14:43:28 stephene Exp stephene $";
#endif


/* Functions relating to the masks used for convolutions. */


/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "dispnet.h"
#include "dispglobals.h"
#include "dispmasks.h"
#include "dispvars.h"
#include "rnd.h"

/* - Defines - */
#define MASK1D 1	/* Should be 1 to indicate that ht of array is 1. */



/* - Function Declarations - */
static void normaliseSumD(double *data, int datasize, double sum);


/* - Global Variables - */ 



/* - Start of Code  - */


void create1dExpMask(double lambda, Mask *rmask)
{
  /* Create a 1D exponential mask.
   * The size of the mask will be 2extent+1,
   * where extent = 4 half lives of the exponential.
   */

  /*** Local Variables ***/
  int	offset;
  Real	*data;
  Mask	mask;
  int	extent, datasize, i;
  Real	coeff;
  int	maxsize;

  
  extent =(int)(4 * log(2)/ lambda);
  printf("%s: Lambda %lf Extent %d\n", __FUNCTION__, lambda, extent);

  maxsize = ( outputWid -1) / 2;

  
  /* In the 1d case, the maximum size of the mask should be outputWid,
   * (eg 1000), and so extent must be no greater than (1000-1)/2.  If
   * the extent is bigger than this value, then extent should be
   * truncated to this size.
   */
  
  if ( extent > maxsize) {
    printf("Extent (%d) must be less than maxsize (%d)\n", extent, maxsize);
    printf("Truncating extent to maxsize\n");

    extent=maxsize;
  }
	   
  /* Create memory for the data */
  datasize = 2 * extent+1;

  data = (Real*)calloc(datasize, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }


  /* now fill in elements of the mask */

  /* Offset is the location of the centre element in the array. */
  offset = extent;
  data[offset+0]=0.0;  /* set centre element to 0. */
  
  /* The 1D mask is symmetrical, so we can just do both sides of
   * the exponential at once. */

  for(i=1; i<= extent; i++) {
    coeff = exp( -lambda * i);
    data[offset+i] = coeff;
    data[offset-i] = coeff;
  }

  /* Normalise the mask so that the sum of the elements is 1.0 */
  normaliseSumD(data, datasize, 1.0);

  /* Save the mask details and return the object. */
  mask.maskExtent = extent;
  mask.wid = datasize;
  mask.ht = MASK1D;
  mask.data = data;
  mask.centre = offset;

  /* Return the mask. */
  *rmask = mask;
}



void  create2dExpMask(double lambda, Mask *rmask)
{
  /* Create a 2D exponential mask.
   * extent of the mask = 4 half lives.
   * width of mask = height of mask = 2extent + 1.
   */
  
  /*** Local Variables ***/
  int		offset, x,y;
  Real		*data, *cdata;
  Mask		mask;
  int		width, height;
  int		extent, datasize, i, maxsize;
  double	dist;
  
  /* There must be a check in creating 1d and 2d masks that the
   * size of the masks do not exceed the size of the input
   * arrays that they will be convolved with.
   */
  
  extent =(int)((4 * log(2.0))/ lambda);
  printf("%s: Lambda %lf Extent %d\n", __FUNCTION__, lambda, extent);


  
  maxsize = ( outputWid -1) / 2;	
  
  
  /* In the 1d case, the maximum size of the mask should be outputWid,
   * (eg 1000), and so extent must be no greater than (1000-1)/2.  If
   * the extent is bigger than this value, then extent should be
   * truncated to this size.
   */
  
  if ( extent > maxsize) {
    printf("Extent (%d) must be less than maxsize (%d)\n", extent, maxsize);
    printf("Truncating extent to maxsize\n");

    extent=maxsize;
  }

  /* Check extent. */

  
  /* Create memory for the data */
  width = 2*extent+1;
  height = width;
  datasize = width * height;


  data = (Real*)calloc(datasize, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }

  /* now fill in elements of the mask */
  /* Offset is the location of the centre element in the array. */
  
  offset = (width*extent) + extent;
  data[offset+0]=0.0;  /* set centre element to 0. */
  
  /* The data array is an array that is bounded from [-extent:extent]
   * in both dimensions.   We will be sweeping across the array on
   * a row by row basis.
   */

  cdata = data; /* Take a copy of the pointer to the data. */
  for(y=-extent; y<=extent; y++) {
    /* Do one row of the array */
    for(x=-extent; x <= extent; x++) {
      dist = sqrt ((double) (x*x) + (y*y));
      *cdata = exp(-dist * lambda);
      cdata++;
    }	
  } /* next row */
  
  /* Normalise the mask so that the sum of the elements is 1.0 */
  normaliseSumD(data, datasize, 1.0);
  
  /* Save the mask details and return the object. */
  mask.maskExtent = extent;
  mask.wid = width;
  mask.ht = height;
  mask.data = data;
  mask.centre = offset;
  *rmask = mask;
}



void freeMasks()
{
  /* Free the mask data structures */
  cfree(uMask.data);
  cfree(vMask.data);
  cfree(uMaskD.data);
  cfree(vMaskD.data);
}




void testMasks()
{
  ulambda =0.1; vlambda = 0.001;
  outputWid = 1000;
  createMasks();
  writeMask(uMask, "umask");
  writeMask(vMask, "vmask");
  writeMask(uMaskD, "umaskD");
  writeMask(vMaskD, "vmaskD");

  freeMasks();
}


void testMasks2()
{
  ulambda =0.1; vlambda = 0.001;
  outputWid = 30; outputHt = 30;
  createMasks2();
  writeMask(uMask, "umask");
  writeMask(vMask, "vmask");
  writeMask(uMaskD, "umaskD");
  writeMask(vMaskD, "vmaskD");

  freeMasks();
}


void createMasks()
{
  /* Create the masks for the convolutions */
  create1dExpMask(ulambda, &uMask);
  create1dExpMask(vlambda, &vMask);

  createDerivMask(uMask, &uMaskD);
  createDerivMask(vMask, &vMaskD);
	  
}


void createMasks2()
{
  /* Create the masks for the convolutions */
  create2dExpMask(ulambda, &uMask);
  create2dExpMask(vlambda, &vMask);

  createDerivMask(uMask, &uMaskD);
  createDerivMask(vMask, &vMaskD);
	  
}

void createDerivMask(Mask mask, Mask *derivMask)
{
  /* Create the mask to be used for creating the derivatives.
   * This is done by subtracting 1.0 from the centre element of 
   * the mask.
   */
  
  /*** Local Variables ***/
  int centre;
  
  copyMask(mask, derivMask);
  centre = derivMask->centre;
  
  /* Subtract 1.0 from centre element */
  derivMask->data[centre] = derivMask->data[centre] - 1.0;
}



void copyMask(Mask mask, Mask *rMask)
{
  /* Copy the MASK into the NEWMASK, creating the memory as needed. */

  /*** Local Variables ***/  
  int wid, ht;
  int extent;
  int datasize;
  int centre, i;
  Real *data;
  Real	*newdata;
  
  wid = mask.wid; ht = mask.ht; extent = mask.maskExtent;
  data = mask.data; centre = mask.centre;
  
  datasize = wid * ht;

  newdata = (Real*)calloc(datasize, sizeof(Real));
  if (! newdata) { 
    printf("%s: could not allocate space for newdata\n", __FUNCTION__);
    exit(-1);
  }

  
  
  rMask->data = newdata; rMask->centre = centre;
  rMask->wid = wid; rMask->ht = ht; rMask->maskExtent = extent;

  /* Copy the data elements over to the new array. */
  
  for(i=datasize; i-->0; ) {
    *newdata = *data;
    newdata++; data++;
  } 
}


void writeMask(Mask mask, char *fname)
{
  /* write the mask MASK out to the file called FNAME.  If the file
   * is given as "-", then it is written to stdout. */
	
  /*** Local Variables ***/
  FILE	*fp;
  int x,y;
  Real *data;
  int	usingstdout=0;
  
  data = mask.data;
  /* get  the relevant file pointer. */
  if (strcmp(fname, "-") == 0) {
    fp = stdout;
    usingstdout=1;
  } else {
    
    fp = fopen( fname, "w");
    if (! fp ) {
      printf("%s: %s could not be opened for writing",
	     __FUNCTION__, fname);
      exit(-1);
    }
  }

  for(y=mask.ht; y-->0; ) {
    for(x=mask.wid; x-->0; ) {
      fprintf(fp, "%lf ", *data++);
    }
    fprintf(fp, "\n");
  }
  /* Close the file. */
  if (!usingstdout) {
    fclose(fp);
  }
}


void normaliseSumD(double *data, int datasize, double sum)
{
  /* Normalise DATA[0..datasize-1] so that the sum of the elements are
   *  equal to SUM.
   *
   * Normalisation is done divisively by dividing each element by
   * VAL. */

  /*** Local Variables ***/
  double	total;
  double	*cdata; /* Copy of the data */
  double	val;
  int		i;
	
  cdata = data; total = 0.0;
  for(i=datasize; i-->0; ) {
    total += *cdata;
    cdata++;
  }
  val = total/sum;
  printf("Total is %lf\tVal is %lf\n to make sum %lf\n",
	 total, val, sum);
  
  /* now divide each element by VAL. */
  cdata = data;
  for(i=datasize; i-->0; ) {
    *cdata /= val;
    cdata++;
  }
}

void testNormalisation()
{
  /* Create a random array and then normalise it. */
  /*** Local Variables ***/  
  double	sum =2.0;
  Array		arr;
  int		i, wid, ht, size;
  Real		*data;

  wid = 4; ht = 3;
  size = wid*ht;

  createRndArray(wid, ht, &arr);
  
  printf("Before normalisation to %lf\n", sum);
  writeArray(arr, "-");
  normaliseSumD(arr.data, size, sum);
  printf("After norm....\n");
  printf("After normalisation to %lf\n", sum);
  writeArray(arr, "-");

  freeArray(arr);
}


/*************************** Version Log ****************************/
/*
 * $Log: dispmasks.c,v $
 * Revision 1.1  1995/11/20  14:43:28  stephene
 * Initial revision
 *
 */

