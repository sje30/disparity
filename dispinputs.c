/****************************************************************************
***
*** Time-stamp: <12 Nov 95 22:31:25 stephene>
***
*** dispinputs.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 12 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/12 23:06:05 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispinputs.c,v 1.1 1995/11/12 23:06:05 stephene Exp stephene $";
#endif

/* Functions to provide the input to the disparity network */

/* -  Include Files - */
#include <stdio.h>
#include "dispnet.h"
#include "rnd.h"

/* - Defines - */
#define INPUTLAYER 0		/* Layer number of the input layer. */

/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */


void getInputVector()
{
  /* Store the next inputVector into the activations array.  */
  /* Bias is set by the function setBiases(), which is called
     separately */

  /* For the moment, we can just create random numbers. */
  /*** Local Variables ***/
  int i;
  Real *actn, *op;
  int offset;
  i = layerInfo[ INPUTLAYER ].ncells;

  offset = actInfo.startLayer[INPUTLAYER]; /* Should be zero. */ 
  actn = &(actInfo.actn[offset]);
  op = &(actInfo.op[offset]);
  for(; i-->0; ) {
    *actn = rnd();
    *op = *actn;
    actn++;
    op++;
  }
}



void setBiases()
{
  /* If bias is required for a layer, this function will set the
   * appropriate bias elements in the activation array and the output
   * array to 1.0 .  */
  int layer;
  int lastInputLayer;
  
  lastInputLayer = netInfo.nLayers-1;
  for(layer=0; layer <lastInputLayer; layer++) {
    if ( layerInfo[layer].bias == Bias) {
      /* Set the activation level */

      /* This will need checking  */
      actInfo.actn[actInfo.biasIndex[layer] ] = 1.0;
      actInfo.op[actInfo.biasIndex[layer] ] = 1.0;
    }
  }
}



void clearActivationArray()
{
  /* Reset the activations array, so that all of the elements are 0.0 */
  int i;
  Real *actn;			/* Pointer to the activations array */
  i = actInfo.size;
  actn = actInfo.actn;
  

  for(;i-- > 0; ) {
    *actn++ = 0.0;
  }
}

void readInData()
{
  /* Read in the data from the files, along with the shift values. */
  ;
}
  
/*************************** Version Log ****************************/
/*
 * $Log: dispinputs.c,v $
 * Revision 1.1  1995/11/12  23:06:05  stephene
 * Initial revision
 *
 */

