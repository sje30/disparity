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
*** $Revision: 1.6 $
*** $Date: 1996/01/16 01:27:42 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /home/stephen/disparity/dispwts.c,v 1.6 1996/01/16 01:27:42 stephene Exp stephen $";
#endif

/********************************************************************/
/* Functions to create and manipulate weights for disparity program */
/********************************************************************/
  
/* -  Include Files - */
#include <stdio.h>
#include <string.h>
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
  int	*preCells, *postCells;
  
  data = (Real*)calloc(len, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }

  preCells = (int*)calloc(len, sizeof(int));
  if (! preCells) { 
    printf("%s: could not allocate space for preCells\n", __FUNCTION__);
    exit(-1);
  }

  postCells = (int*)calloc(len, sizeof(int));
  if (! postCells) { 
    printf("%s: could not allocate space for PostCells\n", __FUNCTION__);
    exit(-1);
  }

  
  weightInfo.data = data;
  weightInfo.nextFreeWeight = 0;
  weightInfo.maxIndex = len-1;
  weightInfo.preCell = preCells;
  weightInfo.postCell = postCells;
}

void showWeights()
{
  /* Print out the weights Info */
  int i;
  int numWeights;

  numWeights = weightInfo.numWts;

  printf("There are %d weights : \n");
  for(i=0; i<numWeights; i++) {
    printf("Weight %d (%lf) connects cell %d to cell %d\n",
	   i, weightInfo.data[i],
	   weightInfo.preCell[i], weightInfo.postCell[i]);
  }
}

  
  

void freeWeights()
{
  /* Clear up the weights structure */
  cfree(weightInfo.data);
  cfree(weightInfo.preCell);
  cfree(weightInfo.postCell);  
  weightInfo.nextFreeWeight = -999;
  weightInfo.maxIndex = -999;
}
  

Real *nextFreeWeight(int preCellNum, int postCellNum)
{
  /* Allocate the next free weight. This weight is to be used to
     connect cell number preCellNum to cell number postCellNum. */
  
  Real *nextwt;
  int nextfree = weightInfo.nextFreeWeight;

  if (nextfree >= weightInfo.maxIndex ) {
    printf("%s: Have exceeded maximum number of weights %d - change dispdefines.h and recompile\n",
	   __FUNCTION__, weightInfo.maxIndex);
    exit(-1);
  }
    
  weightInfo.preCell[nextfree] = preCellNum;
  weightInfo.postCell[nextfree] = postCellNum;
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


void readWts(char *fname)
{
  /* Read in the weights from the file FNAME into the weights data structure. */
  
  /*** Local Variables ***/  
  int i;
  Real *data, val;
  FILE	*fp;

  fp = fopen( fname, "r");
  if (! fp ) {
    printf("%s: %s could not be opened for reading",
	   __FUNCTION__, fname);
    exit(-1);
  }

  i = weightInfo.numWts;

  data = weightInfo.data;
  for(;i-->0;) {
    fscanf(fp, "%lf\n", &val);
    *data++ = val;
  }	
  fclose(fp);
}



void printWtsInfo()
{
  /* Print out the weight information to the file weight.info */
  /* This function adapted from createPartials() */
  
  /*** Local Variables ***/  
  int weight, numWeights;
  char wtInfo[30];
  FILE	*fp;
  int preCell, postCell;

  numWeights = weightInfo.numWts;
  strcpy(wtInfo, "weight.info");  /* file for storing weight info */


  fp = fopen( wtInfo, "w");
  if (! fp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, wtInfo);
    exit(-1);
  }
  
  for(weight=0; weight<numWeights; weight++) {
    preCell = weightInfo.preCell[weight];
    postCell = weightInfo.postCell[weight];

    
    fprintf(fp,"weight %d connects cell %d to cell %d\n",
	     weight, preCell, postCell);

  }
  
  fclose(fp);

}
/*************************** Version Log ****************************/
/*
 * $Log: dispwts.c,v $
 * Revision 1.6  1996/01/16 01:27:42  stephene
 * now have the code in place so that weight sharing can be done or left
 * out.
 *
 * Revision 1.5  1995/12/13  04:04:47  stephene
 * Added new function void printWtsInfo() to create weight.info file.
 *
 * Revision 1.4  1995/12/11  06:26:16  stephene
 * *** empty log message ***
 *
 * Revision 1.3  1995/11/21  23:33:09  stephene
 * About to include CG Code
 *
 * Revision 1.2  1995/11/17  00:06:02  stephene
 * New randomise procedure for the weights provided by Jim.
 *
 * Revision 1.1  1995/11/12  23:06:13  stephene
 * Initial revision
 *
 */

