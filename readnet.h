/****************************************************************************
***
*** Time-stamp: <>
***
*** readnet.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef _READNET_H
#define _READNET_H

/* Function definitions for readnet.c
 * This file can be included by all other source files. */

typedef double Real;


void readNet(char *fname);
int cellArea(int tlx, int tly, int brx, int bry);
void createWeights(int len);
void freeWeights();
Real *nextFreeWeight();



/* Do we need bias for a layer? */
typedef  enum  { nobias, bias} BiasType;

/* Possible types of activation function */
typedef  enum {linearfn, tanhfn} ActFn;

typedef struct {
  int  numInputs;
  Real *inputs;
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

#define AND &&
#define OR ||

#endif


/*************************** Version Log ****************************/
/* $Log$
 *
 */
