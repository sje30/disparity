/****************************************************************************
***
*** Time-stamp: <12 Nov 95 23:09:11 stephene>
***
*** dispnet.c
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

/* Main code file for the Disparity net. */

/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "dispinputs.h"
#include "dispnet.h"

/* - Defines - */

#define WTFORMATn "%lf\n"


/* - Function Declarations - */
void showGnuplot();
void setUpGnuplot();
void showActivations (int layer, char *fname);
void showAllActivations();


/* - Global Variables - */ 



/* - Start of Code  - */




void calcAllActivations()
{
  int layer;
  int firstOPlayer;
  int nLayers;
  firstOPlayer = 1;
  
  nLayers = netInfo.nLayers;
  /* clear activations array from last iteration */
  /* only do this if the next input is not in place */	
  clearActivationArray();
  getInputVector();
  setBiases();
  for(layer=firstOPlayer; layer< nLayers; layer++) {
    calcActivation(layer);
  }
}


void calcActivation(int layer)
{
  /* Calculate the activation levels of cells in layer LAYER. */

  /*** Local Variables ***/
  Real **ptrInputs;
  Real op;
  int actFN;
  int cell, offset;
  int ncells;
  Real *wtStart;
  double sum;
  Real wt, in, res;
  cellInfo_t cellInfo;
  int numInputs;
  int	*inputs;		/* the array of inputs to a cell */
  
  ncells = layerInfo[layer].ncells;
  offset = actInfo.startLayer[layer]; /* starting location of
  					 activation units for this
  					 layer */
  printf("Calculating activation of cells in layer %d\n", layer);
  /* Loop over each cell in the destLayer. */
  for(cell=0; cell<ncells; cell++) {
    cellInfo = layerInfo[layer].cellInfo[cell];
    wtStart = cellInfo.wtsStart	;
    numInputs = cellInfo.numInputs;
    inputs = cellInfo.inputs;
#undef slowsum
#ifdef slowsum    
    sum = 0.0;
    for(;numInputs-- > 0; ) {	
/*       in = **ptrInputs; */
      in = actInfo.actn[*inputs];
      wt = *wtStart;
      res = in * wt;
      printf ("input %lf * wt %lf = %lf\n", in, wt, res);
      sum += res;
      wtStart++;
/*       ptrInputs++; */
      inputs++;
    }
#else
    sum = 0.0;
    ptrInputs = cellInfo.ptrInputs;
    for(;numInputs-- > 0; ) {	
      sum += *wtStart++ * **ptrInputs;
      ptrInputs++;
    }
#endif
    printf("Actn cell %d is %lf\n", cell, sum);

    actFN = layerInfo[layer].actfn;
    /* Now pass through the activation function of the cells */
    switch (actFN) {

    case Linearfn:	
      /* do nothing for a linear unit. */
      /* do nothing */
      ;
      break;
    case    Tanhfn:
      /* Tanh activation function */

      op = tanh(sum);
      sum = op;
      break;
      
    }
    printf("Op cell %d is %lf\n", cell, sum);
  
  /* Store the activation of the cell */

    actInfo.actn[offset+cell] = sum;
  } /* next unit*/
}

void iteration()
{
  setUpGnuplot();
  clearActivationArray();
  setBiases();
  getInputVector();
  calcAllActivations();
  showAllActivations();
  showGnuplot();
}

void showAllActivations()
{
  /* print out all the activation info */
  /* Produces files layerN.dat, where N =0..nlayers-1 */

  char fname[80];
  int layer, nlayers;
  nlayers = netInfo.nLayers;

  for(layer=0; layer < nlayers; layer++) {
    sprintf(fname,"layer%d.act", layer);
    showActivations(layer, fname);
  }
}


void showActivations2(int layer, char *fname)
{
  /* Simply output all activation values */
  int i;
  Real *actn;
  FILE	*fp;
  char op[80];
  i = actInfo.size;

  strcpy(op, "op");
  fp = fopen( op, "w");
  if (! fp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, op);
    exit(-1);
  }
  actn =actInfo.actn;
  for(;i--> 0; ) {
    fprintf(fp, "%lf ", *actn);
    actn++;
  }
  fclose(fp);
}

void showActivations(int layer, char *fname)
{
  /* Write out the activations of layer LAYER to the file called
     fname. */

  /*** Local Variables ***/

  int unit;
  int ncells;
  FILE	*fp;
  int offset;


  fp = fopen( fname, "w");
  if (! fp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, fname);
    exit(-1);
  }
  offset = actInfo.startLayer[layer];
  ncells = layerInfo[layer].ncells;
  for(unit=0; unit<ncells; unit++) {
    fprintf(fp,"%lf\n", actInfo.actn[offset+unit]);
  }

  if (layerInfo[layer].bias ==Bias) {
    fprintf(fp,"%lf\n", actInfo.actn[ actInfo.biasIndex[layer] ]);
  }
  /* Finish up */
  fclose(fp);
}

void clearUpMemory()
{
  /* Clear up the memory at the end of the program */

  /* Should clear up all the space allocated to the cellInfo structures */
  freeWeights();
  cfree(layerInfo);
  freeActivationsArray();
}

void setUpGnuplot()
{
  /* Set up gnuplot for displaying activations */
  system("gplot -create l0 -geometry 400x200+0+0");
  system("gplot -create l1 -geometry 400x200+0+226");
  system("gplot -create l2 -geometry 400x200+0+455");
  system("gplot -create wts -geometry 400x200+0--15");
}
     

void showGnuplot()
{
  system("gplot -to l0 -with linespoints layer0.act");
  system("gplot -to l1 -with linespoints layer1.act");
  system("gplot -to l2 -with linespoints layer2.act");
  system("gplot -to wts -with linespoints wts0.wts");	
}
/*************************** Version Log ****************************/
/*
 * $Log$
 */

