/****************************************************************************
***
*** Time-stamp: <11 Nov 95 18:15:04 stephene>
***
*** dispnet.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision: 1.4 $
*** $Date: 1995/11/13 22:10:41 $
****************************************************************************/


#ifndef _DISPNET_H
#define _DISPNET_H

/* Header file for the disparity network */

/* By default, we will be using doubles for weights and activations */
typedef double Real;


/* Do we need bias for a layer? */
typedef  enum  { Nobias, Bias} BiasType;

/* Possible types of activation function */
typedef  enum {Linearfn, Tanhfn} ActFn;

typedef struct {
  int  numInputs;
  int *inputs;
  Real **ptrInputs;
  Real *wtsStart;
} cellInfo_t;

/* cellInfo_t:	Tells a cell where its input is coming from, and how
 *		many inputs there are.  inputs[i] stores the location
 *		of the input cell in the activations array that is
 *		part of the input to this cell.
 *		ptrInputs stores the pointer to the input cell.
 *		i.e. ptrInput[i] = & (actInfo.op[ inputs[i] ]);
 */

/* There is one cellInfo structure for each cell in the op layer, but
   there is no need to have cellInfo structure for the bias units, as
   these cells do not receive any input. */



/*** Presynaptic cell information ***/

typedef struct {
  int nOutputs;			/* number of connections that this
 				   presynaptic cell makes. Initially
 				   will be zero, and then incremented */
  Real **wts;			/* wts[i] stores a pointer to the ith
 				   fan out weight from this cell*/
  int *outputs;			/* output[i] stores the number of the
 				   postsynaptic cell that this weight
 				   is connected to. */
				/* IS THIS Absolute or relative? */

  Real **ptroutputs;		/* ptrouputs[i] is a pointer to the
  				   activation of the ith postsynaptic
  				   cell */
} preCellInfo_t;

 
typedef struct {
  int tlx;	int tly;	/* Top Left */
  int brx;	int bry;	/* Bottom Right */
} Rect;

   
/* layerInfo_t	Information for a layer. */

typedef struct {
  int       	ncells;		/* this equals nrows * ncols. It
				 * therefore does not include any Bias
				 * weight, which will need to be
				 * explicitly checked for.*/
				   
  int       	nrows;		/* same as the ht */
  int       	ncols;		/* same as the wid */
  cellInfo_t 	*cellInfo;	/* cellInfo[i] stores the input
				   information for cell i in this
				   layer */
  ActFn     	actfn;		/* Which activation function cells in this
				   layer use. */
  BiasType	 bias;		/* Does this layer provide bias for
				   the next layer? */
  preCellInfo_t	*preCellInfo;	/* preCellInfo[i] stores the
  				   presynaptic info for cell i of this
  				   layer */
  int		nPreCellInfo;	/* number of elements in the
  				 * preCellInfo array. This value can
  				 * either be ncells or ncells+1,
  				 * depending on whether there is a
  				 * bias cell in this layer or not. */
  
} layerInfo_t;



/* weightInfo_t:	Information about the weight structure */
typedef struct {
  Real	*data;		/* Actual weight vector */
  int	nextFreeWeight;	/* Index to next free weight element, as
			 * weights are being allocated by the function
			 *  nextFreeWeight() */			   
  int	maxIndex;	/* Maximum index into the weight vector */
  int	numWts;		/* Number of weights allocated.
			 * ie, last weight stored in data[numWts-1] */
} weightInfo_t;



typedef struct {
  Real *actn;
  Real *op;
  int  *startLayer;
  int  *biasIndex;
  int  size;
} activationInfo_t;

/* ActivationInfo:	Stores the activation level of the cells 
 * 
 * actn			actn[i] stores the activation level of the ith cell
 *
 * op			op[i] stores the output of the ith cell.
 *
 * startLayer 		startLayer[i] stores the index to the start of the
 *			cells for the ith layer
 * 
 * biasIndex		biasIndex[i] stores the index to the ith bias unit.
 *			If there is no bias unit, then this value should
 *			equal 99999 to try and cause an error.
 *
 * size			size of the activation array, ie. actn[0 to size-1].
 * */





/* netInfo_t:	General information about the network.
 *
 * nLayers -  number of layers in the network. These will be labelled
 *            from 0 to the Nlayers-1. Layer 0 is the input layer.  */

typedef struct {
  int nLayers;
} netInfo_t;



/*** Function definitions ***/
void calcAllActivations();
void calcActivation(int layer);
void clearUpMemory();

/************************* Global Variables *************************/
weightInfo_t	weightInfo;
netInfo_t	netInfo;
activationInfo_t actInfo;
layerInfo_t	*layerInfo;
/************************* Global Variables *************************/



#endif


/*************************** Version Log ****************************/
/* $Log: dispnet.h,v $
 * Revision 1.4  1995/11/13  22:10:41  stephene
 * New structure preCellInfo_t to store pre synaptic cell Information
 *
 * Revision 1.3  1995/11/12  23:06:27  stephene
 * Daily change
 *
 * Revision 1.2  1995/11/10  15:29:10  stephene
 * Daily update - next major step is to present net with input and
 * calculate activations.
 *
 * Will also need to initialise weights.
 *
 * Revision 1.1  1995/11/09  20:59:27  stephene
 * Initial revision
 *
 *
 */
