/****************************************************************************
***
*** Time-stamp: <>
***
*** dispnet.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/09 20:59:27 $
****************************************************************************/


#ifndef _DISPNET_H
#define _DISPNET_H

/* Header file for the disparity network */

/* By default, we will be using doubles for weights and activations */
typedef double Real;


/* Do we need bias for a layer? */
typedef  enum  { nobias, bias} BiasType;

/* Possible types of activation function */
typedef  enum {linearfn, tanhfn} ActFn;

typedef struct {
  int  numInputs;
  int *inputs;
  Real *wtsStart;
} CellInfo;


typedef struct {
  int      ncells;
  int      nrows;
  int      ncols;
  CellInfo *cellInfo;
  ActFn    actfn;
  BiasType bias;
} LayerInfo;


typedef struct {
  Real	*data;
  int	nextfreeweight;
  int	maxindex;
} WeightInfo;



typedef struct {
  Real *actn;
  int  *startLayer;
  int  *biasIndex;
  int  size;
} ActivationInfo;

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
 *
 */
    
 
typedef struct {
  int tlx;
  int tly;
  int brx;
  int bry;
} Rect;

typedef struct {
  int nLayers;
} NetInfo;


/* NetInfo:	General information about the network.
 *
 * nLayers -  number of layers in the network. These will be labelled
 *            from 0 to the Nlayers-1. Layer 0 is the input layer.  */



/************************* Global Variables *************************/
WeightInfo weightInfo;
NetInfo	netInfo;
ActivationInfo actInfo;
LayerInfo *layerInfo;
/************************* Global Variables *************************/

#define AND &&
#define OR ||




#endif


/*************************** Version Log ****************************/
/* $Log: dispnet.h,v $
 * Revision 1.1  1995/11/09  20:59:27  stephene
 * Initial revision
 *
 *
 */
