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
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header$";
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
  

void initWtsRnd()
{
  /* Initialise weights to random values */
  int i;
  Real *data;
  data = weightInfo.data;
  i = weightInfo.numWts;
  for(;i-- > 0; ) {
    *data = drnd();
    data++;
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
 * $Log$
 */

