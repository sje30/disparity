/****************************************************************************
***
*** Time-stamp: <12 Nov 95 22:31:25 stephene>
***
*** dispinputs.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 12 Nov 95
***
*** $Revision: 1.3 $
*** $Date: 1995/11/21 23:31:48 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispinputs.c,v 1.3 1995/11/21 23:31:48 stephene Exp stephene $";
#endif

/* Functions to provide the input to the disparity network */

/* -  Include Files - */
#include <stdio.h>
#include "dispnet.h"
#include "dispglobals.h"
#include "dispvars.h"
#include "dispinputs.h"
#include "rnd.h"

/* - Defines - */
#define INPUTLAYER 0		/* Layer number of the input layer. */
#define NUMEYES 2		/* Two images, left and right. */

/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */


void getInputVector()
{
  /* Store the next inputVector into the activations array.  */
  /* Bias is set by the function setBiases(), which is called
     separately */

  /* For the moment, we can just create random numbers. */
  /*** Local Variables ***/
  int i;
  Real *actn, *op;
  int offset;
  i = layerInfo[ INPUTLAYER ].ncells;

  offset = actInfo.startLayer[INPUTLAYER]; /* Should be zero. */ 
  actn = &(actInfo.actn[offset]);
  op = &(actInfo.op[offset]);
  for(; i-->0; ) {
    *actn = rnd();
    *op = *actn;
    actn++;
    op++;
  }
}


void getNextInputVector(int vecnum)
{
  /* Store the inputVector numbered VECNUM into the activations array.
   * Bias is set by the function setBiases(), which is called
   * separately.
   */
   

  /*** Local Variables ***/
  int	i;
  Real	*actn, *op;
  int	offset;
  Real	*ivdata;		/* Pointer to the inputvector data. */
  
  i = layerInfo[ INPUTLAYER ].ncells;

  ivdata = &(inputs.data[vecnum *inputs.wid]);
  
  offset = actInfo.startLayer[INPUTLAYER]; /* Should be zero. */ 
  actn = &(actInfo.actn[offset]);
  op = &(actInfo.op[offset]);

  for(; i-->0; ) {
    *actn = *ivdata;
    *op = *actn;		/* Assuming identity function here for the
				 * input cells*/
    actn++;
    op++;
    ivdata++;
  }
}



void setBiases()
{
  /* If bias is required for a layer, this function will set the
   * appropriate bias elements in the activation array and the output
   * array to 1.0 .  */
  int layer;
  int lastInputLayer;
  
  lastInputLayer = netInfo.nLayers-1;
  for(layer=0; layer <lastInputLayer; layer++) {
    if ( layerInfo[layer].bias == Bias) {
      /* Set the activation level */

      /* This will need checking  */
      actInfo.actn[actInfo.biasIndex[layer] ] = 1.0;
      actInfo.op[actInfo.biasIndex[layer] ] = 1.0;
    }
  }
}



void clearActivationArray()
{
  /* Reset the activations array, so that all of the elements are 0.0
   * for both the activations and the outputs*/
  
  int i;
  Real *actn;			/* Pointer to the activations array */
  Real *op;			/* Pointer to the output array. */
  i = actInfo.size;
  actn = actInfo.actn;
  op = actInfo.op;

  for(;i-- > 0; ) {
    *actn++ = 0.0;
    *op++ = 0.0;
  }
}

void readInData()
{
  /* Read in the data from the files, along with the shift values. */
  ;
}

void freeInputsAndShifts()
{
  cfree(shifts.data);
  cfree(inputs.data);
}


void testInputVectors()
{
  /* test procedure for reading input vectors and shifts. */


  createInputVectorsAndShifts();
  writeArray(shifts, "shifts.test");
  writeArray(inputs, "inputs.test");

  freeInputsAndShifts();
}



void createInputVectorsAndShifts()
{
  /* Create the input vectors to be supplied to the network, and 
   * the array of shifts.
   * Sat Dec  9 1995.
   * Modifiying routine so that it can read in inputs from a 2d image
   * in multiple rows, rather than assuming input image is one long row.
   */

  /*** Local Variables ***/  
  int i;
  int colstart, rowstart;
  Array leftVector, rightVector;
  Array leftImage, rightImage, shiftsArr;
  int im, elem;
  int vec; /* Current vector being created. */
  int inputCentreRow, inputCentreCol;
  int oneInputSize = inputWid * inputHt;    
  int inputVectorSize =  NUMEYES * inputWid*inputHt; 
  Real *shiftsData, *inputsdata, *data;


  createArray( inputWid, inputHt, &leftVector);
  createArray( inputWid, inputHt, &rightVector);
  
  createArray( totalInputWid, totalInputHt, &leftImage);
  createArray( totalInputWid, totalInputHt, &rightImage);
  createArray( totalInputWid, totalInputHt, &shiftsArr);
  
  /*** The following two arrays are global variables ***/
  createArray( inputVectorSize, numInputVectors, &inputs);
  createArray( numInputVectors, 1, &shifts);
  shiftsData = shifts.data;
  
  readInputFile( image1File, leftImage);
  readInputFile( image2File, rightImage);
  readInputFile( shiftsFile, shiftsArr);

  /* The first input vector is to be read from the top left hand
   *  corner of the left and right image.
   *
   * Consecutive inputs are then taken by moving colstart along the
   * image until it gets to the end of a row.  Once at the end of a
   * row, if more images are needed, then rowstart is incremented to
   * move onto some more data.  */
     
  colstart = 0;
  rowstart = 0;
  inputsdata = inputs.data;
  
  /* Find the centre position of the input vector */
  inputCentreRow = (inputHt -1)/2;
  inputCentreCol = (inputWid-1)/2;
  
#ifdef debug
  printf("%s: Centre of input vector: row %d col %d\n", inputCentreRow,
	 inputCentreCol);
#endif
  
  /* Extract the input vectors */
  for(vec=0; vec< numInputVectors; vec++) {
    
    extractInputVector(leftImage, leftVector, colstart, rowstart);
    extractInputVector(rightImage, rightVector, colstart, rowstart);
	
    /* Get the relevant shift, by taking the shift value that is
     * at the centre of the inputVector.
     */
    /* The relevant shift element is at 
     * column colstart + inputCentreRow,
     * row    rowstart + inputCentreCol
     */

    elem = ( (rowstart+inputCentreRow) * shiftsArr.wid) +
      (colstart + inputCentreCol);
    shifts.data[vec] = shiftsArr.data[elem];

    /* Now store these vectors in the inputs array. */
    for(im=0; im< NUMEYES; im++) {
      switch(im) {
      case 0:
	data = leftVector.data;
	break;
      case 1:
	data = rightVector.data;
	break;
      default:
	printf("Unknown value of im: %d\n", im);
	exit(-1);
	break;
      }

      for(i=oneInputSize; i-->0;) {
	*inputsdata++ = *data++;
      }
    } /* next image. */

    /* move onto the next vector */
    colstart += inputWid + inputSkipX;

    if ( ( (colstart+inputWid) > leftImage.wid) &&
	(vec != (numInputVectors - 1) ) ) /* not last input vector */
      {
	/* We have gone over the edge of the input images, and so we
	 * need to move onto the next row of the image.
	 * Reset colstart and increment rowstart. */
	
	colstart = 0;
	rowstart += inputHt + inputSkipY;
	
	/* Check to see if rowstart is valid, ie. we have not gone over
	   the edge of the images. */
	if ( ( (rowstart+inputHt) > leftImage.ht) &&
	    (vec != (numInputVectors - 1) ) )/* not last input vector */ {
	      printf("%s: Error:  Will be reading beyond the image. Reading InputVector %d, rowstart %d\n", __FUNCTION__, vec, rowstart);
	      exit(-1);
	    }
      }
  } /* next vec */
  
  /*** Have now finished creating the input vectors. We can now
   * free the arrays that are no longer needed.  */
  
  cfree(leftVector.data); cfree(rightVector.data);
  cfree(leftImage.data); cfree(rightImage.data);
  cfree(shiftsArr.data);
  
}


void extractInputVector(Array source, Array dest,
			int colstart, int rowstart)
{
  /* Fill in the DEST array with elements from the SOURCE array,
   * with the top left hand corner of the DEST array aligned
   * to position (colstart, rowstart) in the SOURCE array. */

  /*** Local Variables ***/
  int destwid, destht;
  Real *destdata;
  int x,y;
  int sourcewid, sourceht, lhsource;
  Real *sourcedata;
  Real *sourceindex;

  destwid = dest.wid; destht = dest.ht; destdata = dest.data;
  sourcewid = source.wid; sourceht = source.ht; sourcedata = source.data;
  /* work out index to the top left hand corner of the source array
   * to be copied from. */
  lhsource = (sourcewid * rowstart) + colstart;	

  for(y=0; y<destht; y++) {
    sourceindex = &sourcedata[lhsource + (y*source.wid)];
    for(x=destwid; x-->0;) {
      *destdata++ = *sourceindex++;
    }
  }
}



void readInputFile(char *inputFile, Array arr)
{

  /* Read in the data from one data file, and then store it 
   * in arr. Space has already been allocated for the array. 
   */

  /* Maybe should check that enough values were read in to
   * the array from the file. */

  /*** Local Variables ***/
  int	x,y;
  FILE	*fp;
  int	wid, ht, size;
  Real	*data;
  Real	val;
  
  wid = arr.wid;
  ht = arr.ht;
  data = arr.data;
  size = wid*ht;
	

  fp = fopen( inputFile, "r");
  if (! fp ) {
    printf("%s: %s could not be opened for reading",
	   __FUNCTION__, inputFile);
    exit(-1);
  }

   for(y=0; y< ht; y++) {
     for(x=0; x < wid; x++) {
       fscanf(fp, "%lf", &val);
       *data++ = val;
     }
     fscanf(fp, "\n");		/* end of line terminator? */
   }

  fclose(fp);
}

/*************************** Version Log ****************************/
/*
 * $Log: dispinputs.c,v $
 * Revision 1.3  1995/11/21  23:31:48  stephene
 * About to include CG Code
 *
 * Revision 1.2  1995/11/17  00:04:26  stephene
 * *** empty log message ***
 *
 * Revision 1.1  1995/11/12  23:06:05  stephene
 * Initial revision
 *
 */
