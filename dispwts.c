/****************************************************************************
***
*** Time-stamp: <12 Nov 95 22:17:01 stephene>
***
*** dispwts.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 12 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/12 23:06:13 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispwts.c,v 1.1 1995/11/12 23:06:13 stephene Exp stephene $";
#endif

/* Functions to create and manipulate weights for disparity program */


  
/* -  Include Files - */
#include <stdio.h>
#include "rnd.h"
#include "dispnet.h"
#include "dispwts.h"
/* - Defines - */



/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */


/********************************************/
/*** Functions for allocating the weights ***/
/********************************************/

void createWeights(int len)
{
  /* Create the weight vector of length len. Initialise other relevant
   * structures */

  Real	*data;

  data = (Real*)calloc(len, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }

  weightInfo.data = data;
  weightInfo.nextFreeWeight = 0;
  weightInfo.maxIndex = len-1;
}


void freeWeights()
{
  /* Clear up the weights structure */
  cfree(weightInfo.data);
  weightInfo.nextFreeWeight = -999;
  weightInfo.maxIndex = -999;
}
  

Real *nextFreeWeight()
{
  /* Allocate the next free weight */
  Real *nextwt;
  nextwt = &(weightInfo.data[weightInfo.nextFreeWeight]);
  weightInfo.nextFreeWeight++;
  return nextwt;
}

void noMoreWeights()
{
  /* Have decided that there are no more weights to be allocated. So
     we can now set the value of numWts. */

  weightInfo.numWts = weightInfo.nextFreeWeight;
}
  

void initWtsRnd1()
{
  /* Initialise weights to random values */
  int i;
  Real *data;
  double r;
  data = weightInfo.data;
  i = weightInfo.numWts;
  for(;i-- > 0; ) {
    *data = 0.3 * drnd();
    /* Choose at random to make the weight negative. */
    r = drnd();
    if (r>0.5) {
      *data *= -1.0;
    }
    data++;
  }
}


void initWtsRnd()
     /* Code for initialising weights from Jim. */
{
  /* according to PW's neural comp paper  2/95 */
  int		num_wts, i;
  Real	*wts;
  Real	r,  A;

  /* What does this give the weights - zero mean, and a certain variance? */

  wts = weightInfo.data;
  num_wts = weightInfo.numWts;
	
  A = 1.0 / sqrt((double)(2.0*num_wts));
  /* numerator = 1 for input wts, and 1.6 for op wts  */ 	
	
  for (i=0;i<num_wts;i++)
    {
      /* r = 0 --> 1 */
      r = drnd();
      
      wts[i] = A * log(r);
      
      r = drnd();
      if (r<=0.5)  r=1; else  r=(-1);
      wts[i] *= r;
    }
}


void writeWts(char *fname)
{
  /* Write the weights out the file called FNAME */
  /*** Local Variables ***/  
  int i;
  Real *data;
  FILE	*fp;

  fp = fopen( fname, "w");
  if (! fp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, fname);
    exit(-1);
  }

  i = weightInfo.numWts;

  data = weightInfo.data;
  for(;i-->0;) {
    fprintf(fp, "%lf\n", *data++);
  }	
  fclose(fp);
}

/*************************** Version Log ****************************/
/*
 * $Log: dispwts.c,v $
 * Revision 1.1  1995/11/12  23:06:13  stephene
 * Initial revision
 *
 */

