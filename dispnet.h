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
*** $Revision: 1.2 $
*** $Date: 1995/11/10 15:29:10 $
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
  Real *wtsStart;
} cellInfo_t;

/* cellInfo_t:	Tells a cell where its input is coming from, and how
 *		many inputs there are.  inputs[i] stores the location
 *		of the input cell in the activations array that is
 *		part of the input to this cell.
 */

/* layerInfo_t	Information for a layer. */

typedef struct {
  int       	ncells;		/* this equals nrows * ncols */
  int       	nrows;		/* same as the ht */
  int       	ncols;		/* same as the wid */
  cellInfo_t 	*cellInfo;	/* cellInfo[i] stores the input
				   information for cell i in this
				   layer */
  ActFn     	actfn;		/* Which activation function cells in this
				   layer use. */
  BiasType	 bias;		/* Does this layer provide bias for
				   the next layer?*/ 
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
  int  *startLayer;
  int  *biasIndex;
  int  size;
} activationInfo_t;

/* ActivationInfo:	Stores the activation level of the cells 
 * 
 * actn			actn[i] stores the activation level of the ith cell
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
    
 
typedef struct {
  int tlx;	int tly;	/* Top Left */
  int brx;	int bry;	/* Bottom Right */
} Rect;



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
