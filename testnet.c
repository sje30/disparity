/****************************************************************************
***
*** Time-stamp: <12 Nov 95 23:06:40 stephene>
***
*** testnet.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision: 1.5 $
*** $Date: 1995/11/21 23:34:02 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/testnet.c,v 1.5 1995/11/21 23:34:02 stephene Exp stephene $";
#endif


/* Test file for testing the network */

/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "dispnet.h"
#include "dispvars.h"
#include "readnet.h"
#include "dispwts.h"
#include "dispglobals.h"
#include "disperrors.h"
#include "testnet.h"
#include "cg_williams_module.h"
#include "bp_check_deriv.h"

/* - Defines - */

/* - Function Declarations - */
extern int	yylex();
double evalMeritFunction();
static void printParams();


/* - Global Variables - */ 



/* - Start of Code  - */
main(int argc, char *argv[])
{

/*  normalmain(argc, argv); */
  netmain(argc, argv);
  

  /*   testMasks2(); */
  /*
     getParams(argv[1]);
     testInputVectors();
     */
  
} /* end of main() */

normalmain(int argc, char *argv[])
{
  /*** Local Variables ***/
  int vecnum;
  extern int	yylex();
  double u,v;
  double u1, v1;	/* u1 = 1.0/U and v1 = 1.0/V */

/*  testArrayDist();   exit(-1); */
  
  if (argc !=2 ) {
    printf("Usage %s <parameter file> \n", argv[0]);
    exit(-1);
  }
  getParams(argv[1]);

  readNet(netFile);
/*   calcAllActivations(); */

  initWtsRnd();
  writeWts("wts0.wts");

  printPreCellInfo();

  createInputVectorsAndShifts();

  writeArray(shifts, "shifts.test");
  writeArray(inputs, "inputs.test");

  /* Check that the size of the input vectors and the number of input
   *  cells is the same. */
  if ( layerInfo[0].ncells != inputs.wid) {
    printf("%s: Error - number of cells in input layer (%d) does not match
the dimensionality of the input vectors (%d)\n",
	   __FUNCTION__, layerInfo[0].ncells, inputs.wid);
    exit(-1);
  }



  setUpGnuplot();
  createZs();


  /* Weights must have been allocated before this routine is called,
   * as it must know the number of weights that were allocated. */
  createdw();
  
  createMasks();
  writeMask(uMask, "umask");
  writeMask(vMask, "vmask");
  writeMask(uMaskD, "umaskD");
  writeMask(vMaskD, "vmaskD");


  for(vecnum=0; vecnum< numInputVectors; vecnum++) {
    /* 
    if ((vecnum % 10 )== 0) {
      printf("Using input vector %d\n", vecnum);
    }
    */
    clearActivationArray();
    setBiases();
    getNextInputVector(vecnum);
    
    calcAllActivations();
    showAllActivations();
/*    showGnuplot(); */

    storeActivations(vecnum);
  }
  getZ();
  writeArray(z, "z");

  /* Create the zbar and ztilde arrays by convolution. */

  double_convolve_wrap(z, uMask, ztilde);
  double_convolve_wrap(z, vMask, zbar);

  writeArray(zbar, "zbar");
  writeArray(ztilde, "ztilde");


  u = arrayDist(z, ztilde);
  v = arrayDist(z, zbar);
  
  printf("U %lf\t V %lf\n", u, v);

  /* zbar min z is simpley zbar - z */
  /* ztilde min z is ztilde - z */
  
  subArray(zbar, z, zbarminz);
  subArray(ztilde, z, ztildeminz);

  writeArray(zbarminz, "zbarminz");
  writeArray(ztildeminz, "ztildeminz");

  /* now create dU/dX and dV/dX by convolution. */
  double_convolve_wrap(ztildeminz,
		       uMaskD,
		       dudx);


  double_convolve_wrap(zbarminz,
		       vMaskD,
		       dvdx);

  writeArray(dudx, "dudx");
  writeArray(dvdx, "dvdx");	

  v1 = 1.0 / v;
  u1 = 1.0 / u;

  printf("1/u %lf 1/v %lf \n", u1, v1);

  multArrayInPlace( dudx, u1);
  multArrayInPlace( dvdx, v1);

  subArray(dvdx, dudx, da);

  writeArray(da, "da");
  
  /* da = dF/dx_a  = 1/v dV/Dx_a - 1/u dU/dx_a */


  printAllActns("allacts");
  storeTopLayerErrors(da);

  propagateErrors();

  printAllActns("allacts2");
  createPartials();
  clearUpMemory();

  /*     iteration(); */
} 



void setUpNetwork(char *paramFile)
{
  /* Set up the network. Read in the parameters from the paramFile and then
   * set up all data structures, read inputs and so on.
   */


  getParams(paramFile);

  readNet(netFile);


  initWtsRnd();
  writeWts("wts0.wts");

  printPreCellInfo();

  createInputVectorsAndShifts();

  writeArray(shifts, "shifts.test");
  writeArray(inputs, "inputs.test");
  
  /* Check that the size of the input vectors and the number of input
   *  cells is the same. */
  if ( layerInfo[0].ncells != inputs.wid) {
    printf("%s: Error - number of cells in input layer (%d) does not match
the dimensionality of the input vectors (%d)\n",
	   __FUNCTION__, layerInfo[0].ncells, inputs.wid);
    exit(-1);
  }


/*   setUpGnuplot(); */
  createZs();


  /* Weights must have been allocated before this routine is called,
   * as it must know the number of weights that were allocated. */
  createdw();

  createMasks();
  writeMask(uMask, "umask");
  writeMask(vMask, "vmask");
  writeMask(uMaskD, "umaskD");
  writeMask(vMaskD, "vmaskD");

  printf("%s: end of function\n", __FUNCTION__);
}

  
netmain(int argc, char *argv[])
{
  /*** Local Variables ***/
  double	f;
  int		fin;
  int		maxWtIndex;
  Real		tol = 1.0;
  int		maxLineSearches= 10;
  
  
  if (argc !=2 ) {
    printf("Usage %s <parameter file> \n", argv[0]);
    exit(-1);
  }

  setUpNetwork( argv[1]);

  maxWtIndex = weightInfo.numWts - 1;
/*  
  f = evalFn(weightInfo.data);
  printf("F has been evaluated to %lf\n", f);
*/


  if (checker) {
    /** Call the function checker. **/
    
    fin = bp_check_func_deriv(weightInfo.data,
			      0,
			      maxWtIndex,
			      evalFn,
			      evalPartials,
			      tol,
			      maxLineSearches);

/*    fin = bp_check_func_deriv(weightInfo.data -1,
			      1,
			      maxWtIndex+1,
			      evalFn,
			      evalPartials,
			      tol,
			      maxLineSearches); */

  }
  else {
    /* Call conjugate gradient */
    fin = cg_williams(weightInfo.data,
		      0,
		      maxWtIndex,
		      evalFn,
		      evalPartials,
		      finishedFn);

  }

  /*** End of training ***/
  clearUpMemory();

  /*     iteration(); */
} 




void calcMeritAndPartials()
{

  /* This function does all the hard work - pass all the inputs
   * through the network; determine U,V and F, and then create the
   * partials vector.
   *
   * The result of the merit function, f, is stored in netInfo.f
   *
   * The partials will be stored in the global array dw.
   *
   */

  /*** Local Variables ***/
  double u,v,f;
  double u1, v1;	/* u1 = 1.0/U and v1 = 1.0/V */
  int vecnum;
  
  /* Evaluate the merit function */

  
  for(vecnum=0; vecnum< numInputVectors; vecnum++) {
    /*
       
    if ((vecnum % 10 )== 0) {
      printf("Using input vector %d\n", vecnum);
    }
    */
    clearActivationArray();
    setBiases();
    getNextInputVector(vecnum);
    
    calcAllActivations();
    showAllActivations();
/*    showGnuplot(); */

    storeActivations(vecnum);
  }
  getZ();
  writeArray(z, "z");

  /* Create the zbar and ztilde arrays by convolution. */

  double_convolve_wrap(z, uMask, ztilde);
  double_convolve_wrap(z, vMask, zbar);

  writeArray(zbar, "zbar");
  writeArray(ztilde, "ztilde");

  /******************************/
  /*** - Calculate U, V and F ***/
  /******************************/
  
  u = arrayDist(z, ztilde);
  v = arrayDist(z, zbar);

/*   printf("Before halving: u %lf v %lf \n", u,v); */
  u *= 0.5;
  v *= 0.5;
/*   printf("After halving: u %lf v %lf \n", u,v); */
  f = log(v/u);

  /* Store f so that it can be retrieved by the other functions.*/
  netInfo.f = f;
  
  printf("U %lf\t V %lf F %lf\n", u, v, f);

  /* zbar min z is simpley zbar - z */
  /* ztilde min z is ztilde - z */
  
  subArray(zbar, z, zbarminz);
  subArray(ztilde, z, ztildeminz);

  writeArray(zbarminz, "zbarminz");
  writeArray(ztildeminz, "ztildeminz");

  /* now create dU/dX and dV/dX by convolution. */
  double_convolve_wrap(ztildeminz,
		       uMaskD,
		       dudx);


  double_convolve_wrap(zbarminz,
		       vMaskD,
		       dvdx);

  writeArray(dudx, "dudx");
  writeArray(dvdx, "dvdx");	

  v1 = 1.0 / v;
  u1 = 1.0 / u;

  printf("1/u %lf 1/v %lf \n", u1, v1);

  multArrayInPlace( dudx, u1);
  multArrayInPlace( dvdx, v1);

  subArray(dvdx, dudx, da);

  writeArray(da, "da");
  
  /* da = dF/dx_a  = 1/v dV/Dx_a - 1/u dU/dx_a */


  printAllActns("allacts");
  storeTopLayerErrors(da);

  propagateErrors();

  printAllActns("allacts2");
  createPartials();

}




double evalFn2(Real *wts)
{
  /*
   * Evaluation function.
   * Given the weights, it evaluates the merit function.
   */

  /*** Local Variables ***/
  int		numWeights;
  double	rval;
  printf("Simple eval fn\n");
  
  rval =1.0;
  return rval;

}

double evalFn(Real *wts)
{
  /*
   * Evaluation function.
   * Given the weights, it evaluates the merit function.
   */

  /*** Local Variables ***/
  int		numWeights;
  double	rval;
  
  numWeights = weightInfo.numWts;
  /* Copy the weights passed into my global variable. */
  copyVec(weightInfo.data, wts, numWeights);

  /* Do the hard work */
  calcMeritAndPartials();

  rval = netInfo.f;

  if (cgmax) {
    /* negate if we are maximising function */
    rval *= -1.0;
  }

  return rval;

}

void evalPartials(Real *wts, Real *derivs)
{
  /* Function to be passed to Conjugate gradient method.
   *
   * First argument is the weights that are to be used.
   * Second argument will store on return the value of the
   * derivatives for these weights.
   */

  /*** Local Variables ***/
  int numWeights;

  printf("eval partials\n");

  /* return; */

  numWeights = weightInfo.numWts;
  /* Copy the weights passed into my global variable. */
  copyVec(weightInfo.data, wts, numWeights);

  /* do the hard work again */
  calcMeritAndPartials();

  /* Now we need to return the partial derivatives. At the moment,
   * these are stored in the array dw.  */

  copyVec(derivs, dw.data, numWeights);

  if (cgmax) {
    /* negate vector if we are maximising function */
    vecNegate(derivs, numWeights);
  }

}


char finishedFn(Real *wts, int iteration)
{
  /* Decide whether we have finished? */
  char rval;
  if (iteration > 20) {
    rval =1;
  }
  else {
    rval = 0;
  }

  return rval;
}
    
void vecNegate(Real *vec, int len)
{
  /* Negate each element of the vector of length LEN. */
  for(; len-->0; ) {
    *vec *= -1.0;
    vec++;
  }
}

  
     
void getParams(char *fname) { /* Get the lex defaults from the file
   FNAME */
  setParamDefaults();
  freopen( fname, "r", stdin);
  
  printf("About to parse %s\n", fname);
  while ( yylex() ) ;
  
  printf("End of lex\n");
  
  printParams();
}
     
void setParamDefaults()
{
  /* Set up the defaults for the parameters. */
  strcpy(netFile, "default.net");
  strcpy(inputFile, "inputs");

  ulambda = 1.0;
  vlambda = 1.0;


  outputWid = 1000;
  outputHt = 1;

  totalInputWid = 7000;
  totalInputHt = 5;

  strcpy(image1File, "leftimage");
  strcpy(image2File, "rightimage");
  strcpy(shiftsFile, "shifts");


  inputWid = 5;
  inputHt = 1;
  inputSkip = 2;
  numInputVectors = 1000;

  cgmax = 1;			/* By default, we will maximise function. */
  checker = 1;			/* Check code rather than use CG to learn. */
}

void printParams()
{
  /* Print the parameters after the parameters file has been loaded */
  printf("*** System Parameters ***\n");
  printf("netFile %s\n", netFile);
  printf("inputFile %s\n", inputFile);
  printf("ulambda = %lf\n", ulambda);
  printf("vlambda = %lf\n", vlambda);
  printf("outputWid = %d\n", outputWid);
  printf("outputHt = %d\n", outputHt);

  printf("totalInputWid = %d\n", totalInputWid);
  printf("totalInputHt = %d\n", totalInputHt);

  printf("image1File = %s\n", image1File);
  printf("image2File = %s\n", image2File); 
  printf("shiftsFile = %s\n", shiftsFile);

  printf("numInputVectors = %d\n", numInputVectors);
  printf("inputSkip = %d\n", inputSkip);
  printf("inputHt = %d\n", inputHt);
  printf("inputWid = %d\n", inputWid);

  printf("cgmax = %d\n", cgmax);
  printf("checker = %d\n", checker); 
}


void copyVec(Real *wdest, Real *wsrc, int numWeights)
{
  /* Copy data from wsrc to wdest.  numWeights elements are copied
   * across.
   */
  for(;numWeights-->0;) {
    *wdest = *wsrc;
    wdest++; wsrc++;
  }
}

  
/*************************** Version Log ****************************/
/*
 * $Log: testnet.c,v $
 * Revision 1.5  1995/11/21  23:34:02  stephene
 * About to include CG Code
 *
 * Revision 1.4  1995/11/21  02:32:27  stephene
 * Update - moving towards a merit function
 *
 * Revision 1.3  1995/11/17  00:07:35  stephene
 * # daily update
 *
 * Revision 1.2  1995/11/10  22:05:10  stephene
 * about to make a snapshot
 *
 * Revision 1.1  1995/11/09  20:49:36  stephene
 * Initial revision
 *
 */

