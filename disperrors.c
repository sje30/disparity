/****************************************************************************
***
***
*** disperrors.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 21 Nov 95
***
*** $Revision: 1.6 $
*** $Date: 1996/01/16 01:27:35 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/disperrors.c,v 1.6 1996/01/16 01:27:35 stephene Exp $";
#endif


/* Store the errors for eachinput and propagate them back to earlier levels.
 */


/* -  Include Files - */

#include <stdio.h>
#include <math.h>
#include "dispnet.h"
#include "disperrors.h"
#include "dispglobals.h"

/* - Defines - */

/* - Function Declarations - */
void copyErrorFromDa(int tlx, int tly, int block, int oplayerht, int oplayerwid);


/* - Global Variables - */ 



/* - Start of Code  - */
#ifdef oldstoretoplayer
void storeTopLayerErrors(Array da)
{
  /* da is the error vector storing the errors for the output cell for
   * each input vector. */
  
  int opLayer = netInfo.nLayers - 1;
  int opCell; /* Index to the output cell in the activation info */
  int vec, numVecs;
  
  opCell = actInfo.startLayer[opLayer];

/*   printf("%s: output is layer %d, op cell is index %d\n", */
/* 	 __FUNCTION__, opLayer, opCell); */

  numVecs = allActns.num;

  for(vec=0; vec< numVecs; vec++) {
    allActns.errors[vec][opCell] = da.data[vec];
  }
}
#else

void storeTopLayerErrors(Array da)
{
  /* da is the error vector storing the errors for the output cell for
   * each input vector. */

  /* Here we will have to loop over each virtual block, and copy a
  ** group of output values back into the allActns.error vectors, -
  ** previously, we just copied one element.
  */

     
  int opLayer = netInfo.nLayers - 1;
  int opCell; /* Index to the output cell in the activation info */
  int vec, numVecs;
  int oplayerht, oplayerwid; /* Dimensions of the output layer. */
  int tlx, tly, block, oplayer;
  int nrows, ncols, i,j;
  
  oplayer = netInfo.nLayers - 1;
  oplayerht = layerInfo[oplayer].nrows;
  oplayerwid = layerInfo[oplayer].ncols;

  nrows = da.ht / oplayerht;
  ncols = da.wid / oplayerwid;
  tlx = 0; tly=0;		/* top left corner of where to start
				 * copying from the da array */
  block =0;			/* virtual block number */
  for(j=0; j< nrows; j++) {
    for(i=0; i< ncols; i++) {
      copyErrorFromDa(tlx, tly, block, oplayerht, oplayerwid);
      tlx += oplayerwid;
      block++;
    }
    /* next row */
    tlx = 0;
    tly += oplayerht;
  }
}

void copyErrorFromDa(int tlx, int tly, int block, int oplayerht, int oplayerwid)
{
  /* Copy part of the array of error vectors da into the error part of
  ** the top layer of cells in the net for virtual block number BLOCK.
  ** The size of the array copied from da is oplayerwid*oplayerht, and
  ** with the top left hand corner place at tlx, tly in the array
  */

  Real *dadata, *nextcellerror;
  int  firstopCell, oplayer;
  int  i,j, nextline;

  
  oplayer = netInfo.nLayers - 1;
  firstopCell = actInfo.startLayer[oplayer];
  nextline = da.wid - oplayerwid;
  
  dadata 		= &(da.data[(tly*da.wid)+tlx]);
  nextcellerror 	= &(allActns.errors[block][firstopCell]);

  for(j=0; j< oplayerht; j++) {
    for(i=0; i< oplayerwid; i++) {
      *nextcellerror++ = *dadata++;
    }
    dadata += nextline; /* move onto the start of the next line */
  }
  
}


#endif


void propagateErrors()
{
  /* Propagate the error measures to the previous layers. */

  /*** Local Variables ***/

  int lastInputLayer = netInfo.nLayers - 2;

  /* Number of the last input layer, ie. the one before the output
   * layer.  It is minus two, because -1 is the number of the output
   * layer.*/
  
  int firstOutputLayer = 1;
  int vec, layer;
  int numVecs;

  numVecs = allActns.num;

  for(vec=0; vec < numVecs; vec++) {
    
    /** Calculate the errors for the input numbered VEC
     **/


    for(layer=lastInputLayer; layer >= lastInputLayer;
	layer--) {
      /** Calculate the errors for layer number LAYER **/

/*       printf("Input %d calc errors for layer %d\n", vec, layer); */

      calcErrors(layer, vec);
    }
  }
}

void calcErrors( int layer, int vec)
{
  /* Calculate the errors using the VECth set of input and activations
   * for layer LAYER. */
  
  /*** Local Variables ***/  
  int unit;
  preCellInfo_t		*preCellInfo;
  int			cellsInLayer;
  Real			**wts;
  int			*outputs;
  int			numOutputs;
  int			*wtsIndex;

  Real			*errors;
  double		sum, bitsum;
  int 			conn, opcell, offset, srcoffset, thiscell;
  Real			wt, error;
  Real			*actns, actn, deriv;
  ActFn			actfn;
  
  /* Work out how many cells there are in this layer.
   * For each cell in this layer:
   *  sum =  sum_k w_kj d_k
   *  sum *= g'(Aj)
   */

  cellsInLayer = layerInfo[layer].nPreCellInfo;
  preCellInfo = layerInfo[layer].preCellInfo;

  /* get the activation function for cells in this layer */
  actfn = layerInfo[layer].actfn;

  /*
  printf("Layer %d cells uses %s activation functions\n", layer,
	 (actfn == Linearfn? "LINEAR" : "TANH" ));
	 */
  
/*   printf("There are %d cells in layer %d\n", cellsInLayer, layer); */

  /* All outputs are relative to the start of the layer, so to get
   * the absolute cell number, we need the offset. */
    
  srcoffset = actInfo.startLayer[layer]; /* Offset for this layer. */
  offset = actInfo.startLayer[layer+1];
  
  for(unit=0; unit<cellsInLayer; unit++) {
    thiscell = unit + srcoffset;
    
    /*printf("Calc error for cell %d\n", unit); */

    numOutputs = preCellInfo[unit].nOutputs;
    /*printf("Unit %d has %d outputs \n", unit, numOutputs); */
    wts = preCellInfo[unit].wts;
    wtsIndex = preCellInfo[unit].wtsIndex;	
    outputs = preCellInfo[unit].outputs;
    errors = allActns.errors[vec];
    actns = allActns.allActs[vec]; /* Array of current activations for
				    * cells in this layer */
				    
				      
    sum = 0.0;
  
    for(conn=0; conn<numOutputs; conn++) {
      wt = *wts[conn];
      opcell = outputs[conn] + offset;
      error = errors[opcell];
      bitsum = error * wt;
      sum += bitsum;

      /*
      printf("Weight index %d Val %lf * Error from output cell %d %lf\n",
	     wtsIndex[conn], wt, opcell, error);
	     */
/*       printf("Bitsum %lf Sum %lf\n", bitsum, sum); */
    }

    /* Pass through derivative function */
    /* Note that this assumes that the cell has a tanh() activation
     * function.  This should be tested by looking at the layerInfo
     * structures.
     */
       
    actn = actns[thiscell];

    /* We need the derivative of the activation function: choose between
     * the possible activation functions */
    
    switch(actfn) {
      
    case Linearfn:
      /* linear unit. Derivative of f(x) = x, is f'(x) = 1.0 */
      deriv = 1.0;
      break;
    case Tanhfn:
      deriv = dtanh(actn);
      break;
    default:
      printf("%s: Error - unknown activation function \n", __FUNCTION__);
      break;

    }
    
      
    /*
    printf("Activation of this cell is %lf - after deriv %lf\n",
	   actn, deriv); */

    sum *= deriv;
    
    /* Store the result */

    errors[thiscell] = sum;
    
  }
    
}


double dtanh( double x)
{
  /* return the derivative of tanh at point x.
   * Calculated using defn of d/dx tanh(x) + 1.0 / cosh^2(x)
   * as this has better saturation properties (see Jim!)
   */
  double c;
  c = cosh(x);
  return (1.0 /(c*c));
}

void testdtanhx()
{
  /* test procedure for dtanh()  - working fine.*/
  double x,y;
  FILE	*fp;
  
  fp = fopen( "dtanh", "w");
  if (! fp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, "dtanh");
    exit(-1);
  }
  
  for(x=-3.0; x<3.0; x+= 0.05) {
    y = dtanh(x);
    fprintf(fp, "%lf %lf\n", x, y);
  }

  fclose(fp);
}

void createPartials()
{
  /* create the whole array of weight changes.
   * The changes are stored in a global called dw
   */

  /*** Local Variables ***/
  char string[80];
  int vec, numVecs;
  int	weight, preCell, postCell, numWeights;
  Real	*outputs, *errors;
  
  /* Clear the dw array. */
  setArray(dw, 0.0);

  numVecs = allActns.num;
  numWeights = weightInfo.numWts;
  
  for(vec=0; vec<numVecs; vec++) {
    outputs = allActns.allOps[vec];
    errors = allActns.errors[vec];
/*     printf("Creating onedw for activations %d\n", vec); */

    setArray(onedw, 0.0);
    for(weight=0; weight<numWeights; weight++) {
      preCell = weightInfo.preCell[weight];
      postCell = weightInfo.postCell[weight];

      /*
      printf("weight %d connects cell %d to cell %d\n",
	     weight, preCell, postCell);
	     */
      onedw.data[weight] = outputs[preCell] * errors[postCell];

      /*
      printf("output %lf * error %lf = %lf\n",
	     outputs[preCell], errors[postCell],
	     onedw.data[weight]);
	     */
      
    }

    /* save this array */
#ifdef dumpArrays
    sprintf(string, "onedw.%d", vec);
    writeArray(onedw, string);
#endif
    
    /* add this result to the dw partials vector. */
    addArrayInPlace(dw, onedw);
  }

  /* dw now created. */
#ifdef dumpArrays
  writeArray(dw, "dw");
#endif

}

  
/*************************** Version Log ****************************/
/*
 * $Log: disperrors.c,v $
 * Revision 1.6  1996/01/16  01:27:35  stephene
 * now have the code in place so that weight sharing can be done or left
 * out.
 *
 * Revision 1.5  1995/12/15  17:13:45  stephene
 * calcErrors modified so that it takes account of the activation
 * function for cells when taking the derivative of the function -
 * previously it was assumed that all cells were tanh() cells, which was
 * ok, given that our hidden layer was tanh() cells
 *
 * Revision 1.4  1995/12/11  06:25:31  stephene
 * *** empty log message ***
 *
 * Revision 1.3  1995/12/07  15:42:34  stephene
 * post nips sort out
 *
 * Revision 1.2  1995/11/23  16:19:33  stephene
 * CG now installed
 *
 * Revision 1.1  1995/11/21  23:30:07  stephene
 * Initial revision
 *
 */
