/****************************************************************************
***
*** Time-stamp: <10 Nov 95 21:20:33 stephene>
***
*** readnet.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision: 1.2 $
*** $Date: 1995/11/10 15:28:36 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/readnet.c,v 1.2 1995/11/10 15:28:36 stephene Exp stephene $";
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
#define TWODTOONED(x,y,wid) ( (wid)*(y) + (x))

/* - Function Declarations - */


/* - Global Variables - */ 

/* Variables private to this file */
static char      line[256];
static int       linenum, linelen;
static FILE      *structFP;

/*********************/
/* Private functions */
/*********************/

static void readNextLine();
static void createActivationArray();
static void showActivationsArray();
static void checkRect(Rect rect, int linenumd, int layer);
static void connectCells(int sourceLayer, int unitnum,
			 int destLayer, Rect rect);
static void printNet();

/* - Start of Code  - */

void readNet(char *fname)

{
  /* read in the network structure from the file fname. */

  /*** Local Variables ***/

  int      destLayer;
  int      layer,cell;
  CellInfo *cellInfo;
  int      unitnum;
  int      numunits;
  int      tlx, tly, brx, bry;
  Rect     rect;
  int      readingConnections;
  int      nextLayer;
  char     afn[30], anybias[30];
  int      i, returnval;
  int      nLayers;
  int      layernum, wid, ht;
  int      ncols, nrows;



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

  readNextLine();
  sscanf(line, "nLayers %d\n", &nLayers);
  netInfo.nLayers = nLayers;

  /* Allocate space for the layer data */

  layerInfo = (LayerInfo*)calloc(nLayers, sizeof(LayerInfo));
  if (! layerInfo) { 
    printf("%s: could not allocate space for layerInfo\n", __FUNCTION__);
    exit(-1);

  }

  /* Read in the size of each layer and store it */

  
  for(i=0; i< nLayers; i++) {
    /*   fgets(line, linelen, structFP); linenum++; */
    readNextLine();
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
    
/*     fgets(line, linelen, structFP); linenum++; */
    readNextLine();
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
  } /* now read in info for the next layer */
  /**************************************************/
  /*** - Allocate space for the activations array ***/
  /**************************************************/


  createActivationArray();
  showActivationsArray();

  /***************************************************/
  /*** Start to allocate the weights to the inputs ***/
  /***************************************************/

  for(nextLayer=0; nextLayer < nLayers; nextLayer++) {
    readingConnections = 1;
    if ( nextLayer == 0 ) {
      /* we dont need to read any connections for the input layer */
      printf("no connections needed for layer 0\n");
      readingConnections = 0;
    }
      
	

    while ( readingConnections ) {

      /*       fgets(line, linelen, structFP); linenum++; */
      readNextLine();
      sscanf(line, "connections to layer %d", &destLayer);
      if ( destLayer != nextLayer) {
	printf("Error: need connections for layer %d - given %d (line %d)\n",
	       nextLayer, destLayer, linenum);
	exit(-1);
      }

      /*** Create the space for the cellInfo structures ***/

      numunits = layerInfo[nextLayer].ncells;
      
      /* Create space for the cellInfo structure */
      
      
      cellInfo = (CellInfo*)calloc(numunits, sizeof(CellInfo));
      if (! cellInfo) { 
	printf("%s: could not allocate space for cellInfo\n",
	       __FUNCTION__);
	exit(-1);
      }
      layerInfo[nextLayer].cellInfo = cellInfo;



/*       fgets(line, linelen, structFP); linenum++; */
      readNextLine();
      if (!strncmp(line, "full",4) ) {
	/* We need full connectivity between the layers */
	printf("Full connectivity from layer %d to layer %d\n", nextLayer-1,
	       nextLayer);
	  

	ncols = layerInfo[nextLayer-1].ncols;
	nrows = layerInfo[nextLayer-1].nrows;

	/* rect will be the size of all of the input layer */
	rect.tlx = 0;
	rect.tly = 0;
	rect.brx = ncols-1;
	rect.bry = nrows-1;
	numunits = layerInfo[nextLayer].ncells;
	for (unitnum=0; unitnum < numunits; unitnum++) {
	  /* xxxx */
	  connectCells( nextLayer-1,  unitnum,  nextLayer, rect);
	}

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
	
	/* The first line of unit connections has already been read
           into line, and so we dont need to do a fgets until we want
           to get the info for the second line */

	for (unitnum=0; unitnum < numunits; unitnum++) {
	  if (unitnum != 0 ) {
	    /* now get the next line*/
 	    readNextLine(); 
/* 	    fgets(line, linelen, structFP); linenum++;  */
	  }

	  returnval = sscanf(line, "%d %d %d %d", &tlx, &tly, &brx, &bry);
	  if ( returnval != 4 ) {
	    printf("%s: Error  reading connectivity in layer %d (line %d)\n",
		   __FUNCTION__, nextLayer,  linenum);
	    exit(-1);
	  }

	  /* check that box coords are ok */
	  
	  rect.tlx = tlx;
	  rect.tly = tly;

	  rect.brx = brx;
	  rect.bry = bry;

	  checkRect(rect, linenum, nextLayer-1);
	  
	  printf("Layer %d unit %d:  %d %d %d %d\n", nextLayer, unitnum,
		 tlx, tly, brx, bry);
	  connectCells( nextLayer-1,  unitnum,  nextLayer, rect);
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


  printNet();

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
  Real *nextwt;
  nextwt = &(weightInfo.data[weightInfo.nextfreeweight]);
  weightInfo.nextfreeweight++;
  return nextwt;
}




void readNextLine()
{
  /* Read the next line from the file, skipping over blank lines. */
  /* Data is returned in the variable line.
   * The variable linenum gives you the current line number. */

  /* This assumes that you are reading from the file pointer structFP */

  int reading = 1; /* flag for reading */
  while (reading) {
    fgets(line, linelen, structFP);
    linenum++;
    /* If line is not just whitespace, then dont read any more yet. */
    if (! emptyLine(line)) {
      reading=0;
    }
  }
}


int emptyLine(char *str)
{
  /* Return true if the string str is empty - ie only has whitespace in it. */
  /* These are just blank lines */
  /* Or ignore the line if the line starts with a hash mark */
  int returnval=1;
  int looping=1;

  /* Simple test that the comment marker # should then ignore the rest
     of the line */


  if (*str =='#') {
    return 1;
  }
  while (looping) {
    if ((*str == ' ') OR (*str == '\t')) {
      /* ok, just whitespace */
      str++;
      ;
    }
    else if (*str == '\n') {
      looping=0; returnval =1;
    } else {
      /* non whitespace char */
      looping=0; returnval=0;
    }
  }
  return returnval;
}
/********************************************************************/
/******* Activation Array Functions: Create, Show and Delete  *******/
/********************************************************************/

void createActivationArray()
{
  /* Create the activation array and set up the other data structures */
  /*** Local Variables ***/
  int cellnum;
  int	*biasIndex;
  int nLayers;
  int layer;
  int	*startLayer;
  int unitsInLayer;

  nLayers = netInfo.nLayers;

  startLayer = (int*)calloc(nLayers, sizeof(int));
  if (! startLayer) { 
    printf("%s: could not allocate space for startLayer\n", __FUNCTION__);
    exit(-1);
  }
  
  biasIndex = (int*)calloc(nLayers, sizeof(int));
  if (! biasIndex) { 
    printf("%s: could not allocate space for biasIndex\n", __FUNCTION__);
    exit(-1);
  }

  cellnum = 0;
  for(layer=0; layer < nLayers; layer++) {
    startLayer[layer] = cellnum;
    unitsInLayer = layerInfo[layer].ncells;
    cellnum += unitsInLayer;

    if ( layerInfo[layer].bias == bias ) {
      /* We have a bias also to include */
      biasIndex[layer] = cellnum;
      cellnum++;
    } else {
      /* This layer has no bias for the next layer */
      biasIndex[layer] = 9999;
    }

  }

  actInfo.size = cellnum;
  actInfo.startLayer = startLayer;
  actInfo.biasIndex = biasIndex;


}

void showActivationsArray()
{
  /* Print out the details of the activation Information */
  int nLayers, layer;
  nLayers = netInfo.nLayers;
  
  printf("Total number of cells: %d\n", actInfo.size);
  
  for(layer=0; layer< nLayers; layer++) {
    printf("Layer %d: Start %d Bias %d\n", layer, actInfo.startLayer[layer],
	   actInfo.biasIndex[layer]);
  }
}

void freeActivationsArray()
{
  /* Free up all of the arrays allocated for the activation Information */

  cfree(actInfo.actn); 
  cfree(actInfo.startLayer); 
  cfree(actInfo.biasIndex);
}


void checkRect(Rect rect, int linenumd, int layer)
{
  /* Check to see that the coordinates of the box are ok. */

  int tlx, tly, brx, bry;
  int ncols, nrows;
  tlx = rect.tlx;
  tly = rect.tly;
  brx = rect.brx;
  bry = rect.bry;

  ncols = layerInfo[layer].ncols;
  nrows = layerInfo[layer].nrows;

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
}


void connectCells(int sourceLayer, int unitnum, int destLayer, Rect rect)
{
  /* Connect Cells from the area rect in the source layer to unit
     UNITNUM of the destination layer. */
  /* This will also take account of the bias units from the source layer
   * into the destLayer. */

  /*** Local Variables ***/
  int  inputcell;
  int  firstInput;
  int  numInputCells, inputnum;
  int *cellinputs;
  int  tlx, tly, brx, bry;
  int  x,y;
  int  ncols;
  int  needbias;
  int	inputOffset;
  CellInfo *cellInfo;


  ncols = layerInfo[sourceLayer].ncols;
  cellInfo = layerInfo[destLayer].cellInfo;

  inputOffset = actInfo.startLayer[sourceLayer];
  tlx = rect.tlx;
  tly = rect.tly;
  brx = rect.brx;
  bry = rect.bry;

  /* Does this cell also need to receive bias input? */
  if ( layerInfo[sourceLayer].bias == bias) {
    needbias =1;
  } else {
    needbias = 0;
  }
    
  numInputCells =  cellArea( tlx, tly, brx, bry);
  if (needbias) {
    numInputCells++;
  }

  printf("Unit %d receiving %d inputs\n", unitnum, numInputCells);
  /* Allocate space to store the input cells to each unit */
  
  cellinputs = (int*)calloc(numInputCells, sizeof(Real));
  if (! cellinputs) { 
    printf("%s: could not allocate space for cellinputs\n",
	   __FUNCTION__);
    exit(-1);
  }
  cellInfo[unitnum].inputs = cellinputs;
  cellInfo[unitnum].numInputs = numInputCells;
  

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
      cellinputs[inputnum] = (inputOffset + inputcell);
      inputnum++;
    }
  }

  if (needbias) {
    /* Add the input from the bias cell */
    cellinputs[inputnum] = actInfo.biasIndex[sourceLayer];
    (void)nextFreeWeight();
  }


}



void printNet()
{
  /* Print out the network */

  /*** Local Variables ***/
  int   activationcell, x,y, nrows, ncols, cellnum;
  int layer, nopcells;
  int cell, numInputs, ip;
  Real	*wt;
  for(layer=0; layer < netInfo.nLayers; layer++) {
    printf("\n\nLayer %d\n",layer);
    
    printf("Units: %d x %d = %d\n", layerInfo[layer].ncols,
	   layerInfo[layer].nrows, layerInfo[layer].ncells);
    printf("Activation: %d Bias %d\n", 	   layerInfo[layer].actfn,
	   layerInfo[layer].bias);
    if (layer != 0 ) {
      nopcells = layerInfo[layer].ncells;
      for(cell=0; cell<nopcells; cell++) {
	numInputs = layerInfo[layer].cellInfo[cell].numInputs;
	printf("Unit %d: Num inputs %d\n", cell,
		 numInputs);
	wt = layerInfo[layer].cellInfo[cell].wtsStart;
	for(ip=0; ip < numInputs; ip++) {
	  printf("Ip %d weight %d\n", layerInfo[layer].cellInfo[cell].inputs[ip], (int)wt);
	  wt++;
	}
      }
    }
  } /* Next layer */

  activationcell = 0;
  for(layer=0; layer < netInfo.nLayers; layer++) {
    ncols = layerInfo[layer].ncols;
    nrows = layerInfo[layer].nrows;
    printf("\n\nLayer %d [0 0 %d %d]\n\n",layer, ncols-1, nrows-1);
    cellnum =0;
    for(y=0; y< nrows; y++) {
      for(x=0; x<ncols; x++) {
	printf("%2d ", cellnum++);
      }
      printf("\n");
    }
    if (layerInfo[layer].bias == bias) {
      printf("Bias\n");
    }
  } /* next layer */

  printf("** Activations **\n");
  for(layer=0; layer < netInfo.nLayers; layer++) {
    printf("\nLayer %d\n",layer);
    ncols = layerInfo[layer].ncols;
    nrows = layerInfo[layer].nrows;
    for(y=0; y< nrows; y++) {
      for(x=0; x<ncols; x++) {
	printf("%2d ", activationcell++);
      }
      printf("\n");
    }
    if (layerInfo[layer].bias == bias) {
      printf("Bias %d\n", activationcell++);
    }
  } /* next layer */
	

}


	  
/*************************** Version Log ****************************/
/*
 * $Log: readnet.c,v $
 * Revision 1.2  1995/11/10  15:28:36  stephene
 * Daily update - next major step is to present net with input and
 * calculate activations.
 *
 * Will also need to initialise weights.
 *
 * Revision 1.1  1995/11/09  20:48:36  stephene
 * Initial revision
 *
 */
