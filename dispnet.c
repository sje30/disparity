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
*** $Revision: 1.4 $
*** $Date: 1995/11/21 02:32:37 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispnet.c,v 1.4 1995/11/21 02:32:37 stephene Exp stephene $";
#endif

/* Main code file for the Disparity net. */

/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "dispinputs.h"
#include "dispvars.h"
#include "dispglobals.h"
#include "dispnet.h"
#include "rnd.h"

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


#undef slowsum

#ifdef slowsum
  printf("Calculating activation of cells in layer %d\n", layer);
#endif
  
  /* Loop over each cell in the destLayer. */
  for(cell=0; cell<ncells; cell++) {
    cellInfo = layerInfo[layer].cellInfo[cell];
    wtStart = cellInfo.wtsStart	;
    numInputs = cellInfo.numInputs;
    inputs = cellInfo.inputs;

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

#ifdef slowsum
    printf("Actn cell %d is %lf\n", cell, sum);
#endif
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
#ifdef slowsum
    printf("Op cell %d is %lf\n", cell, sum);
#endif
  
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
  freeAllActns();
  freePreCellInfo();
  freedw();
  freeZs();
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


/********************* Array Creation Functions *********************/
void createArray( int wid, int ht, Array *rarr)
{
  /* Create an array of size wid*ht. Allocate memory, and return array
   * in *RARR. */
	
  int size;
  Real *data;


  size = wid*ht;
  data = (Real*)calloc(size, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }
  
  rarr-> data = data; rarr->ht = ht; rarr->wid = wid;
}


void createRndArray( int wid, int ht, Array *rarr)
{
  /* Create an array and fill it with random values */
  Real *data;
  Array	arr;
  int	i, size;

  createArray(wid, ht, &arr);
  
  data = arr.data;
  size = wid*ht;
  
  for(i=0; i<size; i++) {
    *data++ = drnd(); 
  }	
  /* Return the array */
  *rarr = arr;
}


void freeArray(Array arr)
{
  /* free the memory allocated to arr. */
  cfree(arr.data);
}


void writeArray(Array arr, char *fname)
{
  /* write the array ARR out to the file called FNAME.  If the file
   * is given as "-", then it is written to stdout. */
	
  /*** Local Variables ***/
  FILE	*fp;
  int x,y;
  Real *data;
  int	usingstdout=0;
  
  data = arr.data;
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

  for(y=arr.ht; y-->0; ) {
    for(x=arr.wid; x-->0; ) {
      fprintf(fp, "%lf ", *data++);
    }
    fprintf(fp, "\n");
  }
  /* Close the file. */
  if (!usingstdout) {
    fclose(fp);
  }
}

/********************** Functions for allActns **********************/

void createAllActns ()
{
  /* Create the arrays to store the activation information for each
   * input vector. */
  /*** Local Variables ***/
  Real	*actdata;
  Real	*opdata;
  Real	*errordata;
  Real	**allOps;
  Real	**allActs;
  Real	**errors;  
  int	numUnits, num, i;

  
  /* each Element is a pointer to an array of activations and a
     pointer to an array of outputs. */

  num = numInputVectors;
  numUnits = actInfo.size;
  
  allActs = (Real **)calloc(num, sizeof(Real *));
  if (! allActs) { 
    printf("%s: could not allocate space for allActs\n", __FUNCTION__);
    exit(-1);
  }

  allOps = (Real **)calloc(num, sizeof(Real *));
  if (! allOps) { 
    printf("%s: could not allocate space for allOps\n", __FUNCTION__);
    exit(-1);
  }


  errors = (Real **)calloc(num, sizeof(Real *));
  if (! errors) { 
    printf("%s: could not allocate space for errors\n", __FUNCTION__);
    exit(-1);
  }



  for(i=0; i< numInputVectors; i++) {
    actdata = (Real*)calloc(numUnits, sizeof(Real));
    if (! actdata) { 
      printf("%s: could not allocate space for actdata\n", __FUNCTION__);
      exit(-1);
    }
    
    opdata = (Real*)calloc(numUnits, sizeof(Real));
    if (! opdata) { 
      printf("%s: could not allocate space for opdata\n", __FUNCTION__);
      exit(-1);
    }

    errordata = (Real*)calloc(numUnits, sizeof(Real));
    if (! errordata) { 
      printf("%s: could not allocate space for errordata\n", __FUNCTION__);
      exit(-1);
    }

    allActs[i] = actdata;
    allOps[i] = opdata;
    errors[i] = errordata;
  }
  allActns.allActs = allActs;
  allActns.allOps= allOps;
  allActns.errors= errors;
  allActns.num = num;
  allActns.numUnits = numUnits;
}


void printAllActns(char *fname)
{


  /* print all the activations out the file fname.
   * If fname = "-", then we use stdout.
   */
  /*** Local Variables ***/
  Real	*ops;
  int numVecs, nCells, vec;
  FILE *fp;
  int	usingstdout=0;
  int	i;

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


  numVecs = allActns.num;
  nCells = allActns.numUnits;
  
  for(vec=0; vec< numVecs; vec++) {

    /* Print the outputs */
    fprintf(fp, "Outputs: \n");
    ops = allActns.allOps[vec];
    for(i=0; i<nCells; i++) {
      fprintf(fp, "%.4lf ", *ops++);
    }
    fprintf(fp, "\n");


    /* Print the activations */
    fprintf(fp, "Activations: \n");
    ops = allActns.allActs[vec];
    for(i=0; i<nCells; i++) {
      fprintf(fp, "%.4lf ", *ops++);
    }
    fprintf(fp, "\n");

    /* Print the errors */
    fprintf(fp, "Errors: \n");
    ops = allActns.errors[vec];
    for(i=0; i<nCells; i++) {
      fprintf(fp, "%.4lf ", *ops++);
    }
    fprintf(fp, "\n");


  } /* next vector */

  
  if (!usingstdout) {
    fclose(fp);
  }

}
     
void freeAllActns()
{
  /* Clear up the allActns data structure. */
  
  int i;
  int num = allActns.num;

  for(i=0; i<num; i++) {
    cfree(allActns.allActs[i]);
    cfree(allActns.allOps[i]);
  }	

  cfree(allActns.allOps);
  cfree(allActns.allActs);
}


void storeActivations(int input)
{
  /* Store the activation information for the input vector numbered INPUT
     into the allActns array. */

  /*** Local Variables ***/
  int	i;
  int	numUnits;
  Real	*op;
  Real	*acts;

  numUnits = allActns.numUnits;

  if (input>allActns.num) {
    printf("%s: Error - input (%d) must be less than allActns.num (%d)\n",
	   __FUNCTION__, input, allActns.num);
    exit(-1);
  }
  op = allActns.allOps[input];
  acts = allActns.allActs[input];

  /* Copy across the activation and output values. */
  for(i=0; i<numUnits; i++) {
    acts[i] = actInfo.actn[i];
    op[i] = actInfo.op[i];
  }
}


void getZ()
{
  /* After all the inputs have been presented, get the Z array */

  /*** Local Variables ***/  
  int vecnum;
  int indexOpCell;		/* Index to the output cell. */
  int oplayer;

  oplayer = netInfo.nLayers - 1;
  
  indexOpCell = actInfo.startLayer[oplayer];
  printf("Op cell is found at location %d\n", indexOpCell);
  
  for(vecnum=0; vecnum < numInputVectors; vecnum++) {
    z.data[vecnum] = allActns.allOps[vecnum][indexOpCell];
  }
}




void createZs()
{
  /* create the z, zbar, and ztilde arrays. */

  createArray( outputWid, outputHt, &z);
  createArray( outputWid, outputHt, &zbar);
  createArray( outputWid, outputHt, &ztilde);
  createArray( outputWid, outputHt, &ztildeminz);
  createArray( outputWid, outputHt, &zbarminz);


  createArray( outputWid, outputHt, &dudx);
  createArray( outputWid, outputHt, &dvdx);
  createArray( outputWid, outputHt, &da);
  
}


void freeZs()
{
  /* free the z, zbar, and ztilde arrays. */

  freeArray(z); freeArray(zbar); freeArray(ztilde);
  freeArray(zbarminz);
  freeArray(ztildeminz);

  freeArray(dudx);   freeArray(dvdx);   freeArray(da);

}

void createdw()
{
  /* create dw and onedw arrays */
  int numWeights;
  numWeights = weightInfo.numWts;
  createArray( numWeights,1, &dw);
  createArray( numWeights,1, &onedw);
}


void freedw()
{
  /* free the dw and onedw arrays */
  freeArray(dw); freeArray(onedw);
}

double arrayDist(Array a1, Array a2)
{
  /* a1 and a2 have the same size.
   * Return value = sum_k (a1[k] - a2[k])^2, k=0..size-1
   */

  /*** Local Variables ***/
  double	sum, diff;
  Real		*a1data, *a2data;
  int		a1size, a2size;
  int		i;
  
  a1data = a1.data; a2data = a2.data;
  a1size = a1.wid * a1.ht;
  a2size = a2.wid * a2.ht;
  if (a1size != a2size) {
    printf("%s: Error - size of arrays are different %d v %d\n",
	   a1size, a2size);
    exit(-1);
  }

  sum = 0.0;
  for(i=a1size; i-->0; ) {
    diff = (*a1data - *a2data);
    sum += diff*diff;
    a1data++; a2data++;
  }
  return sum;
  
}


void testArrayDist()
{
  /* test procedure for arrayDist. */
  Array a1, a2;
  double dist;
  
  int wid =3;
  int ht =2;

  createRndArray(wid, ht, &a1);
  createRndArray(wid, ht, &a2);

  writeArray(a1, "-");
  writeArray(a2, "-");

  dist = arrayDist(a1, a2);
  printf("dist is %lf\n", dist);
}

void subArray(Array a1, Array a2, Array result)
{
  /* result[k] = a1[k] - a2[k] for all elements k */

  /*** Local Variables ***/
  int size = a1.wid * a1.ht;
  Real *a1data, *a2data, *resultdata;

  a1data = a1.data;   a2data = a2.data;  resultdata = result.data;

  for(; size-- > 0; ) {
    *resultdata = (*a1data - *a2data);
    a1data++; a2data++;
    resultdata++;
  }

}


void testSubArray()
{
  /* test procedure for subArray */
  Array a1, a2;
  Array diff;
  
  int wid =3;
  int ht =2;

  createRndArray(wid, ht, &a1);
  createRndArray(wid, ht, &a2);
  createRndArray(wid, ht, &diff);

  subArray( a1, a2, diff); writeArray(a1, "-"); writeArray(a2, "-");
  writeArray(diff, "-");

  freeArray(a1);
  freeArray(a2);
  freeArray(diff);
}


void multArrayInPlace(Array a1, double k)
{
  /* Multiply all array elements by k. */

  /*** Local Variables ***/
  int size = a1.wid * a1.ht;
  Real *data;

  data = a1.data;

  for(; size-- > 0; ) {
    *data *= k;
    data++;
  }

}

void testMult()
{
  /* test the mult procedure for arrays. */
  
  Array a1;
  
  
  int wid =3;
  int ht =2;

  createRndArray(wid, ht, &a1);
  writeArray(a1, "-");

  multArrayInPlace(a1, 2.0);
  printf("After:\n");
  writeArray(a1, "-");

  freeArray(a1);
}



void setArray(Array a1, double k)
{
  /* Set all elements of array to the value k. */

  /*** Local Variables ***/
  int size = a1.wid * a1.ht;
  Real *data;

  data = a1.data;

  for(; size-- > 0; ) {
    *data = k;
    data++;
  }

}


void testSetArray()
{
  /* test the setArray() function */
  /* working ok */
  
  Array a1;
  
  int wid =3;
  int ht =2;

  createRndArray(wid, ht, &a1);
  writeArray(a1, "-");

  setArray(a1, 2.0);
  printf("After:\n");
  writeArray(a1, "-");

  setArray(a1, 0.0);
  printf("After:\n");
  writeArray(a1, "-");

  freeArray(a1);

}


void addArrayInPlace(Array a1, Array a2)
{
  /* new a1[k] + a2[k] -> a1[k], for each element k.
   * old value of a1 is overwritten.
   */
  

  /*** Local Variables ***/
  int size = a1.wid * a1.ht;
  Real *a1data, *a2data;

  a1data = a1.data;   a2data = a2.data;

  for(; size-- > 0; ) {
    *a1data += *a2data;
    a1data++; a2data++;
  }
}


void testAddArrayInPlace()
{
  /* test procedure for addArryInPlace() */
  /* fine */
  Array a1, a2;

  int wid =3;
  int ht =2;

  createRndArray(wid, ht, &a1);
  createRndArray(wid, ht, &a2);


  writeArray(a1, "-"); writeArray(a2, "-");
  addArrayInPlace(a1, a2);
  printf("New a1 :\n");
  writeArray(a1, "-");

  freeArray(a1);
  freeArray(a2);
}


/*************************** Version Log ****************************/
/*
 * $Log: dispnet.c,v $
 * Revision 1.4  1995/11/21  02:32:37  stephene
 * Update - moving towards a merit function
 *
 * Revision 1.3  1995/11/17  00:04:58  stephene
 * Daily update
 *
 * Revision 1.2  1995/11/13  22:14:48  stephene
 * Daily Change
 *
 * Revision 1.1  1995/11/12  23:37:48  stephene
 * Initial revision
 *
 */

