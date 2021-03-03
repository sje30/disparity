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
*** $Revision: 1.11 $
*** $Date: 1998/03/24 17:43:10 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /home/stephen/disparity/dispinputs.c,v 1.11 1998/03/24 17:43:10 stephen Exp $";
#endif

/* Functions to provide the input to the disparity network */

/* -  Include Files - */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dispnet.h"
#include "dispglobals.h"
#include "dispvars.h"
#include "dispinputs.h"
#include "rnd.h"

/* - Defines - */
#define INPUTLAYER 0		/* Layer number of the input layer. */
#define NUMEYES 2		/* Two images, left and right. */


#define RAD_TO_DEG 57.29577951
#define NoOrientation -100 /* This is used to indicate that a region does
			    * not have a valid orientation. */
/* - Function Declarations - */
extern int Globify(char *fname);
Real getOrientation(Array source, FILE *fp);
Real sqr(Real x);
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
  ActFn	actfn;
  
  /* get the activation function for cells in this layer */
  actfn = layerInfo[INPUTLAYER].actfn;

  i = layerInfo[ INPUTLAYER ].ncells;

  ivdata = &(inputs.data[vecnum *inputs.wid]);
  
  offset = actInfo.startLayer[INPUTLAYER]; /* Should be zero. */ 
  actn = &(actInfo.actn[offset]);
  
  op = &(actInfo.op[offset]);

  if (actfn != Linearfn) {
    printf("%s: Error - can only cope with linear input units at the moment - you will have to change this routine if you want to pass the input activation through an activation function\n", __FUNCTION__);
    exit(-1);
  }
    
  /* Fri Dec 15 1995
   * check what the activation function is for these units: */

  for(; i-->0; ) {
    *actn = *ivdata;
    *op = *actn;		/* Assuming identity function here for
				 * the input cells - change here if
				 * you want to code for other
				 * activation functions.  eg *op =
				 * tanh(*actn). See calcErrors() for
				 * an example of how to do this.*/
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
  ActFn			actfn;

  /* modifying 15/12/95 to check for activation function. */
  
  lastInputLayer = netInfo.nLayers-1;
  for(layer=0; layer <lastInputLayer; layer++) {
    if ( layerInfo[layer].bias == Bias) {
      /* Set the activation level */


      actInfo.actn[actInfo.biasIndex[layer] ] = 1.0;
      /* The output of the cell should be the activation passed
       * through the activation function. */
      actfn = layerInfo[layer].actfn;

      /*
	 printf("Layer %d cells uses %s activation functions\n", layer,
	 (actfn == Linearfn? "LINEAR" : "TANH" ));
	 */
      
      switch(actfn) {
	
      case Linearfn:
	/* linear unit. */
	actInfo.op[actInfo.biasIndex[layer] ] = 1.0;
	break;
      case Tanhfn:
	actInfo.op[actInfo.biasIndex[layer] ] = tanh(1.0);
	break;
      default:
	printf("%s: Error - unknown activation function \n", __FUNCTION__);
	break;

      }
      /*	 
      printf("OP of bias cell in layer %d = %lf\n", layer, 	actInfo.op[actInfo.biasIndex[layer] ]); 
      */

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

    if (noshifting) {
      shifts.data[vec] = shiftsArr.data[vec];
    } else {
      /* default behaviour preserved for disparity expts. */
      elem = ( (rowstart+inputCentreRow) * shiftsArr.wid) +
	(colstart + inputCentreCol);
      shifts.data[vec] = shiftsArr.data[elem];
    }

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

void createInputVectorsAndShiftsOneImage()
{
  /* Create the input vectors to be supplied to the network, and 
   * the array of shifts.
   * Sat Dec  9 1995.
   * Modifiying routine so that it can read in inputs from a 2d image
   * in multiple rows, rather than assuming input image is one long row.
   *
   * Sun Jan 14 1996
   * This function is used when we only want one image to be read in.
   */

  /*** Local Variables ***/  
  int i;
  int colstart, rowstart;
  Array leftVector;
  Array leftImage, shiftsArr;
  int im, elem;
  int vec; /* Current vector being created. */
  int inputCentreRow, inputCentreCol;
  int oneInputSize = inputWid * inputHt;    
  int inputVectorSize;
  Real *shiftsData, *inputsdata, *data , sum, sumsq, k, m ;

  /* Orientation stuff */
  char *ornfile = "realorientation.dat";
  FILE	*orfp;
  Real orientation;

  inputVectorSize=  1 * inputWid*inputHt; /* Here we only have one eye */
  createArray( inputWid, inputHt, &leftVector);
  
  createArray( totalInputWid, totalInputHt, &leftImage);
  createArray( totalInputWid, totalInputHt, &shiftsArr);
  
  /*** The following two arrays are global variables ***/
  createArray( inputVectorSize, numInputVectors, &inputs);
  createArray( numInputVectors, 1, &shifts);
  shiftsData = shifts.data;
  
  readInputFile( image1File, leftImage);
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



  orfp = fopen( ornfile, "w");
  if (! orfp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, ornfile);
    exit(-1);
  }

  
  /* Extract the input vectors */
  for(vec=0; vec< numInputVectors; vec++) {
    
    extractInputVector(leftImage, leftVector, colstart, rowstart);
    orientation = getOrientation(leftVector,orfp);
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

    data = leftVector.data;

    /* copy across the data from the image into the inputs. */
    if (normInput) {

      /* normalise input data here */
      sum=0.0 ; sumsq = 0.0 ;
      for (i=0;i<oneInputSize;i++)
	{ sum += data[i] ; sumsq += (data[i]*data[i]) ; } ;
      m = sum/oneInputSize ;
      k = 1.0/(sqrt(sumsq-(sum*sum/oneInputSize))+0.000000001) ;

      /* copy across the data from the image into the inputs. */
      for(i=oneInputSize; i-->0;) {
	*inputsdata++ = k * (*data++ - m)  ;
      }
    }
    else {
      /* don't bother normalising input */
      
      for(i=oneInputSize; i-->0;) {
	*inputsdata++ = *data++;
      }
    }
    
  
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
    
  }
  /*** Have now finished creating the input vectors. We can now
   * free the arrays that are no longer needed.  */
  
  cfree(leftVector.data);
  cfree(leftImage.data); 
  cfree(shiftsArr.data);


  fclose(orfp);
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



Real getOrientation(Array source, FILE *fp)
{
  /* Fill in the DEST array with elements from the SOURCE array,
   * with the top left hand corner of the DEST array aligned
   * to position (colstart, rowstart) in the SOURCE array. */

  /*** Local Variables ***/
  int x,y;
  int sourcewid, sourceht;
  Real *sourcedata;
  Real *sourceindex;

  

  int		i;
  int		numPoints;

  Real		sum00 = 0.0,	/* sumij = x^i, y^j moment. */
  		sum01 = 0.0,
  		sum10 = 0.0,
  		sum11 = 0.0,
  		sum20 = 0.0,
  		sum02 = 0.0;
  Real	xBar, yBar;
  Real	numerator, denominator, orientation;
  Real	xd,yd;
  Real	a,b,c, eMin, eMax, sq1, sinTerm, cosTerm, roundness;
  Real  pixelthresh=0.2;	/* min. value for a pixel to be accepted
				 * in calculation. */

  sourcewid = source.wid; sourceht = source.ht; sourcedata = source.data;
  numPoints = 0;
  
  for (y=0; y < sourceht; y++) {
    for (x=0; x < sourcewid; x++) {

      if ( *sourcedata++ > pixelthresh) {
	numPoints++;
	sum00 += 1;  	/* *B(x,y) */
	sum01 += (1*y);  	/* *B(x,y) */
	sum10 += (1*x);  	/* *B(x,y) */
      }
    }
  }

  xBar = ( sum10 / sum00);
  yBar = ( sum01 / sum00);

  /* Now pass over data again to calculate the 2nd order moments, since we
   * now know the Centre of Gravity of the region. */

  sourcedata = source.data;

  for (y=0; y < sourceht; y++) {
    for (x=0; x < sourcewid; x++) {

      if ( *sourcedata++ > pixelthresh) {
	
	xd = x-xBar;
	yd = y-yBar;
	sum20 += (sqr(xd));  	/* *B(x,y) */
	sum11 += (xd)*(yd);  /* *B(x,y) */
	sum02 += (sqr(yd));  	/* *B(x,y) */
      }
    }
  }
  

  numerator = 2.0 * sum11;
  denominator = sum20 - sum02;


  /* If the numerator and denominator are both zero, then we have a region
   * with no meaningful orientation.  Such objects are circles and squares. */

  if ( (numerator == 0 ) && ( denominator == 0 )) {
    /* no orientation . */
    orientation = NoOrientation;
    roundness = 1;
  }
  else {
    /* orientation is valid. */
    orientation = (atan2(numerator, denominator)/ 2.0);

    /* orientation angle is 2*theta, so we need to divide by 2. */

    /* convert orientation from radians to degrees. */
    orientation *= RAD_TO_DEG;


    a = sum20;
    b = 2.0 * sum11;
    c = sum02;

    sq1 = sqrt( sqr(b) + sqr(a-c));
    sinTerm = b/sq1;
    cosTerm = (a-c)/sq1;
    
    eMin = (0.5*(a+c)) - (0.5*(a-c)*cosTerm) - (0.5*b*sinTerm);
    eMax = (0.5*(a+c)) + (0.5*(a-c)*cosTerm) + (0.5*b*sinTerm);
    
    roundness = eMin / eMax;
  }



#ifdef printMoments
  if (logFP) {
    fprintf(logFP, "Moments\ns00\t%lf\ns01\t%lf\ns10\t%lf\n",
	    sum00, sum01, sum10);
    fprintf(logFP, "Moments\ns20\t%lf\ns11\t%lf\ns02\t%lf\n",
	    sum20, sum11, sum02);
  }
#endif

#ifdef printMoments
  if (logFP)
    fprintf(logFP, "a = %lf b= %lf c = %lf sq1 = %lf\nsinTerm = %lf cosTerm = %lf \neMin = %lf eMax = %lf roundness = %lf\n",
	    a, b, c, sq1, sinTerm, cosTerm, eMin, eMax, roundness);
#endif

  /*
  stats.xBar = (int)xBar;
  stats.yBar = (int)yBar;
  stats.orientation = orientation;
  stats.roundness = roundness;
  */

  fprintf(fp, "%f %f %d\n", orientation, roundness, numPoints);

  return orientation;
  
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
  char	globnm[255];
  wid = arr.wid;
  ht = arr.ht;
  data = arr.data;
  size = wid*ht;
	
  
  if (inputFile[0] == '~') {
    strcpy(globnm, inputFile);
    Globify(globnm);
    inputFile = globnm;
    printf("%s: filename expanded to %s\n", __FUNCTION__, inputFile);
  }

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


Real sqr(Real x)
/* Post - Return the square of x. */
{
  return (x*x);
}
