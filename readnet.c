/****************************************************************************
***
*** Time-stamp: <09 Nov 95 20:48:32 stephene>
***
*** readnet.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header$";
#endif


/* Code to read the network structure from a file and create the
 * neccessary structures for the activations and weights */


/* TO do: sort out connectivity for the "full" marker
 * look at the inclusion of bias weights */


/* - Include Files - */
#include <stdio.h>
#include "dispnet.h"
#include "readnet.h"
/* - Defines - */
#define MAX_NUM_WEIGHTS 2000

/* - Function Declarations - */


/* - Global Variables - */ 
WeightInfo weightInfo;



/* - Start of Code  - */

void readNet(char *fname)

{
  /* read in the network structure from the file fname. */

  /*** Local Variables ***/
  int layer,cell;
  CellInfo  *cellInfo;
  Real      *cellinputs;
  int       numInputCells, inputnum;
  int       firstInput;
  int       x,y;
  int       inputcell;
  int       unitnum;
  int       numunits;
  int       tlx, tly, brx, bry;
  int       readingConnections;
  int       nextLayer;
  char      afn[30], anybias[30];
  FILE      *structFP;
  int       i, returnval;
  int       nLayers;
  int       layernum, wid, ht;
  char      line[256];
  int       linenum, linelen;
  int       ncols, nrows;
  LayerInfo *layerInfo;


  createWeights(MAX_NUM_WEIGHTS);

  structFP = fopen( fname, "r");
  if (! structFP ) {
    printf("%s: %s could not be opened for reading",
	   __FUNCTION__, fname);
    exit(-1);
  }

  /* initialize line data. */

  linelen = sizeof(line);
  linenum = 0;  /* First line of the file will be line 1 iff linenum
		   initialised to 0 */

  fgets(line, linelen, structFP); linenum++;
  sscanf(line, "nLayers %d\n", &nLayers);


  /* Allocate space for the layer data */

  layerInfo = (LayerInfo*)calloc(nLayers, sizeof(LayerInfo));
  if (! layerInfo) { 
    printf("%s: could not allocate space for layerInfo\n", __FUNCTION__);
    exit(-1);

  }

  /* Read in the size of each layer and store it */

  
  for(i=0; i< nLayers; i++) {
  fgets(line, linelen, structFP); linenum++;
  returnval = sscanf(line, "layer %d %d %d\n", &layernum,  &wid, &ht);
    printf("return value = %d\n", returnval);
    if ( returnval != 3 ) {
      printf("%s: Error  reading layer %d data from %s (line %d): only found %d elements\n", 
	     __FUNCTION__, layernum, fname, linenum,  returnval);
      exit(-1);
    }

  layerInfo[layernum].ncols = wid;
  layerInfo[layernum].nrows = ht;
  layerInfo[layernum].ncells = wid*ht;
  }

  /***  Now read in the details about how the layers are connected ***/

  for( nextLayer=0; nextLayer<nLayers; nextLayer++) {
    /* Read in details for a layer. */
    
    fgets(line, linelen, structFP); linenum++;
    returnval = sscanf(line, "layer %d %s %s", &layernum, afn, anybias);
    if ( returnval != 3 ) {
      printf("%s: Error  reading layer %d afn and bias from %s (line %d)\n",
	     __FUNCTION__, nextLayer, fname, linenum);
      exit(-1);
    }
    
    /*** Convert the activation function string into a corresponding code ***/
    if (!strcmp(afn, "TANH") ) {
      printf("Tanh activation function for layer %d\n", nextLayer);
      layerInfo[nextLayer].actfn = tanhfn;
    }
    else if ( !strcmp( afn, "LINEAR")) {
      printf("Linear activation function for layer %d\n", nextLayer);
      layerInfo[nextLayer].actfn = linearfn;
    }
    
    else {
      /* Unknown activation function */
      printf("Unknown activation function %s on line %d\n", afn, linenum);
      exit(-1);
    }
    
    
    /*** Decide whether we need bias or not. ***/
    if (!strcmp(anybias, "bias")) {
      printf("Bias needed for layer %d\n", nextLayer);
      layerInfo[nextLayer].bias = bias;
    }
    else   if (!strcmp(anybias, "nobias")) {
      printf("No bias needed for layer %d\n", nextLayer);
      layerInfo[nextLayer].bias = nobias;
    }
    else {
      /* Invalid bias value */
      printf("Unknown bias value %s on line %d\n", anybias, linenum);
      exit(-1); 
    }

    /***************************************************/
    /*** Start to allocate the weights to the inputs ***/
    /***************************************************/


    readingConnections = 1;
    if ( nextLayer == 0 ) {
      /* we dont need to read any connections for the input layer */
      printf("no connections needed for layer 0\n");
      readingConnections = 0;
    }

    while ( readingConnections ) {
      fgets(line, linelen, structFP); linenum++;
      if (!strncmp(line, "full",4) ) {
	/* We need full connectivity between the layers */
	printf("Full connectivity from layer %d to layer %d\n", nextLayer-1,
	       nextLayer);
	readingConnections = 0;
      } else {
	/* read in specific connectivity details */
	printf("%s : specific connections given from layer %dto %d\n", 
	       __FUNCTION__, nextLayer-1, nextLayer);

	/* ncols and nrows refer to the dimensionality of the cells in
	   the layer providing input to this one. */

	ncols = layerInfo[nextLayer-1].ncols;
	nrows = layerInfo[nextLayer-1].nrows;
	/* we will need to read in a set of conections for 
	   every weight in this layer */

	numunits = layerInfo[nextLayer].ncells;

	/* Create space for the cellInfo structure */
	

	cellInfo = (CellInfo*)calloc(numunits, sizeof(CellInfo));
	if (! cellInfo) { 
	  printf("%s: could not allocate space for cellInfo\n",
		 __FUNCTION__);
	  exit(-1);
	}
	layerInfo[nextLayer].cellInfo = cellInfo;

	
	/* The first line of unit connections has already been read
           into line, and so we dont need to do a fgets until we want
           to get the info for the second line */

	for (unitnum=0; unitnum < numunits; unitnum++) {
	  if (unitnum != 0 ) {
	    /* now get the next line*/
	    fgets(line, linelen, structFP); linenum++; 
	  }

	  returnval = sscanf(line, "%d %d %d %d", &tlx, &tly, &brx, &bry);
	  if ( returnval != 4 ) {
	    printf("%s: Error  reading connectivity in layer %d (line %d)\n",
		   __FUNCTION__, nextLayer,  linenum);
	    exit(-1);
	  }

	  /* check that box coords are ok */

	  if ( (tlx <0 ) OR (tlx >= ncols)) {
	    printf("%s : Error line %d - tlx (%d) must be in range [%d,%d]\n",
		   __FUNCTION__, linenum, tlx, 0, ncols-1);
	    exit(-1);
	  }

	  if ( (brx <0 ) OR (brx >= ncols)) {
	    printf("%s : Error line %d - brx (%d) must be in range [%d,%d]\n",
		   __FUNCTION__, linenum, brx, 0, ncols-1);
	    exit(-1);
	  }

	  if ( (tly <0 ) OR (tly >= nrows)) {
	    printf("%s : Error line %d - tly (%d) must be in range [%d,%d]\n",
		   __FUNCTION__, linenum, tly, 0, nrows-1);
	    exit(-1);
	  }

	  if ( (bry <0 ) OR (bry >= nrows)) {
	    printf("%s : Error line %d - bry (%d) must be in range [%d,%d]\n",
		   __FUNCTION__, linenum, bry, 0, nrows-1);
	    exit(-1);
	  }
	  /* end checks for box */
	  
	  printf("Layer %d unit %d:  %d %d %d %d\n", nextLayer, unitnum,
		 tlx, tly, brx, bry);

	  numInputCells =  cellArea( tlx, tly, brx, bry);
	  printf("Unit %d receiving %d inputs\n", unitnum, numInputCells);

	  /* Allocate space to store the input cells to each unit */

	  cellinputs = (Real*)calloc(numInputCells, sizeof(Real));
	  if (! cellinputs) { 
	    printf("%s: could not allocate space for cellinputs\n",
		   __FUNCTION__);
	    exit(-1);
	  }
	  cellInfo[unitnum].inputs = cellinputs;
	  cellInfo[unitnum].numInputs = numInputCells;
	  
#define TWODTOONED(x,y,wid) ( (wid)*(y) + (x))
	  inputnum = 0;
	  firstInput =1;
	  for(y=tly; y<=bry; y++) {
	    for(x=tlx; x<=brx; x++) {
	      inputcell = TWODTOONED(x,y,ncols);
	      if (firstInput) {
		firstInput=0;
		cellInfo[unitnum].wtsStart = nextFreeWeight();
	      }
	      else {
		(void)nextFreeWeight();
	      }
	      cellinputs[inputnum] = inputcell;
	      inputnum++;
	    }
	  }
	  
	} /* get next unit */

	readingConnections = 0;
      }
    }

    
  } /* now read in the details for the next layer */

  /* End of reading in net - close up files  */

  /* Show how many weights have been allocated. */
  printf("%d weights have been allocated\n", weightInfo.nextfreeweight);
  fclose(structFP);




  /* Print out the information in a coherent manner. */
  for(layer=0; layer < nLayers; layer++) {
    printf("\n\nLayer %d\n",layer);
    
    printf("Units: %d x %d = %d\n", layerInfo[layer].ncols,
	   layerInfo[layer].nrows, layerInfo[layer].ncells);
    printf("Activation: %d Bias %d\n", 	   layerInfo[layer].actfn,
	   layerInfo[layer].bias);
    if (layer != 0 ) {
      for(cell=0; cell<layerInfo[layer].ncells; cell++) {
	printf("Unit %d: Num inputs %d\n", cell, layerInfo[layer].cellInfo[cell].numInputs);
      }
    }   
  }

  /* This must be moved elsewhere at the end of the day */
  cfree(layerInfo);
}



int cellArea(int tlx, int tly, int brx, int bry)
{
  /* Work out the cell area bounded by the rectangle */
  int wid, ht;
  wid = (brx - tlx) + 1;
  ht  = (bry - tly) + 1;
  return (wid*ht);
}

/********************************************/
/*** Functions for allocating the weights ***/
/********************************************/

void createWeights(int len)
{
  /* Create the weight vector of length len. Initialise other relevant
   * structures */

  Real	*data;

  data = (Real*)calloc(len, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }

  weightInfo.data = data;
  weightInfo.nextfreeweight = 0;
  weightInfo.maxindex = len-1;
}


void freeWeights()
{
  /* Clear up the weights structure */
  cfree(weightInfo.data);
  weightInfo.nextfreeweight = -999;
  weightInfo.maxindex = -999;
}
  

Real *nextFreeWeight()
{
  /* Allocate the next free weight */
  REAL *nextwt;
  nextwt = &(weightInfo.data[weightInfo.nextfreeweight]);
  weightInfo.nextfreeweight++;
  return nextwt;
}







/*************************** Version Log ****************************/
/*
 * $Log$
 */
