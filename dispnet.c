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
*** $Revision: 1.2 $
*** $Date: 1995/11/13 22:14:48 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispnet.c,v 1.2 1995/11/13 22:14:48 stephene Exp stephene $";
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
void showOutputs(int layer, char *fname);


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
#define slowsum
#ifdef slowsum    
    sum = 0.0;
    for(;numInputs-- > 0; ) {	
/*       in = **ptrInputs; */
      in = actInfo.op[*inputs];
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

    /* Store the activation of the cell */
    actInfo.actn[offset+cell] = sum;
    
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
  
    /* Store the output of the cell */
    actInfo.op[offset+cell] = sum;

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

    sprintf(fname,"layer%d.op", layer);
    showOutputs(layer, fname);

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


void showOutputs(int layer, char *fname)
{
  /* Write out the outputs of layer LAYER to the file called
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
    fprintf(fp,"%lf\n", actInfo.op[offset+unit]);
  }

  if (layerInfo[layer].bias ==Bias) {
    fprintf(fp,"%lf\n", actInfo.op[ actInfo.biasIndex[layer] ]);
  }
  /* Finish up */
  fclose(fp);
}

void clearUpMemory()
{
  /* Clear up the memory at the end of the program */

  /* Should clear up all the space allocated to the cellInfo structures */
  freeWeights();

  freeActivationsArray();
  freePreCellInfo();
  freeCellInfo();		/* Free up cellInfo before layerInfo */
  cfree(layerInfo);
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
  system("gplot -to l1 -with linespoints layer1.act layer1.op");
  system("gplot -to l2 -with linespoints layer2.act layer2.op");
  system("gplot -to wts -with linespoints wts0.wts");	
}

/*** Conversion between half life and lambda ***/

/*** Jim - why have you used log2 rather than just natural logarithm
     for these functions ? ***/

Real half_life2lambda(Real h)
{		
  Real lambda;
	
  lambda = pow(2.72,(-1.0/h));
  return(lambda);
}

Real lambda2half_life(Real lam)
{	
  Real h, L;
  L = log(lam);
  h = -1.0/L;
  return(h);
}

/*************************** Version Log ****************************/
/*
 * $Log: dispnet.c,v $
 * Revision 1.2  1995/11/13  22:14:48  stephene
 * Daily Change
 *
 * Revision 1.1  1995/11/12  23:37:48  stephene
 * Initial revision
 *
 */

