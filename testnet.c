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
*** $Revision: 1.11 $
*** $Date: 1995/12/11 06:27:37 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/testnet.c,v 1.11 1995/12/11 06:27:37 stephene Exp stephene $";
#endif


/* Test file for testing the network */

/* -  Include Files - */
#include <stdio.h>
#include <math.h>
#include "dispnet.h"
#include "dispvars.h"
#include "readnet.h"
#include "dispwts.h"
#include "dispmasks.h"
#include "dispglobals.h"
#include "disperrors.h"
#include "testnet.h"
#include "cg_williams_module.h"
#include "bp_check_deriv.h"

/* - Defines - */
#define dumpCorrn

/* - Function Declarations - */
void createLayerOpFiles();
void checkNetSize();
extern int	yylex();
double evalMeritFunction();
static void printParams();
void checkNetPerformance();

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
#ifdef dumpArrays
  writeArray(z, "z");
#endif
  
  /* Create the zbar and ztilde arrays by convolution. */

  double_convolve_wrap(z, uMask, ztilde);
  double_convolve_wrap(z, vMask, zbar);
#ifdef dumpArrays
  writeArray(zbar, "zbar");
  writeArray(ztilde, "ztilde");
#endif

  u = arrayDist(z, ztilde);
  v = arrayDist(z, zbar);
  
  printf("U %lf\t V %lf\n", u, v);

  /* zbar min z is simpley zbar - z */
  /* ztilde min z is ztilde - z */
  
  subArray(zbar, z, zbarminz);
  subArray(ztilde, z, ztildeminz);

#ifdef dumpArrays
  writeArray(zbarminz, "zbarminz");
  writeArray(ztildeminz, "ztildeminz");
#endif
  
  /* now create dU/dX and dV/dX by convolution. */
  double_convolve_wrap(ztildeminz,
		       uMaskD,
		       dudx);


  double_convolve_wrap(zbarminz,
		       vMaskD,
		       dvdx);

#ifdef dumpArrays  
  writeArray(dudx, "dudx");
  writeArray(dvdx, "dvdx");	
#endif
  
  v1 = 1.0 / v;
  u1 = 1.0 / u;

/*   printf("1/u %lf 1/v %lf \n", u1, v1); */

  multArrayInPlace( dudx, u1);
  multArrayInPlace( dvdx, v1);

  subArray(dvdx, dudx, da);
#ifdef dumpArrays
  writeArray(da, "da");
#endif
  /* da = dF/dx_a  = 1/v dV/Dx_a - 1/u dU/dx_a */


#ifdef dumpArrays
  printAllActns("allacts");
#endif
  
  storeTopLayerErrors(da);

  propagateErrors();


#ifdef dumpArrays
  printAllActns("allacts2");
#endif
  
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

/*   printPreCellInfo();  */


  if (oneImage) {
    createInputVectorsAndShiftsOneImage();
  } else {
    createInputVectorsAndShifts();
  }


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

  if (usegnuplot) {
    setUpGnuplot();
  }
  
  createZs();


  /* Weights must have been allocated before this routine is called,
   * as it must know the number of weights that were allocated. */
  createdw();

  /************************/
  /*** Create the masks ***/
  /************************/

  /* See if we want to convert half lifes to lambda */
  if (useHalf) {
    ulambda = half2lambda(uhalf);
    vlambda = half2lambda(vhalf);
  }

  printf("Ulambda %lf\tVlambda %lf\n", ulambda, vlambda);
  
    
  /* At this point, we need to decide whether to set up 2d masks or 1d
   * masks.  This will be done according to the value of outputHt.
   */

  if (outputHt == 1) {
    /* 1d Masks needed */
    createMasks();
  }
  else {
    /* 2d masks needed */
    createMasks2();
  }
    
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


  /* Check the network configuration */
  checkNetSize();
  
  /* Decide if we want to do learning, or just test the performance of
   * the network on the input images */

  if ( !doLearning) {
    checkNetPerformance();
    
    exit(-1);
  }

  /* Ok, we must want to do learning... */
  opfp = fopen( "dispoutputs", "w");
  if (! opfp ) {
    printf("%s: %s could not be opened for writing",
	   __FUNCTION__, "dispoutputs");
    exit(-1);
  }

  
  maxWtIndex = weightInfo.numWts - 1;
/*  
  f = evalFn(weightInfo.data);
  printf("F has been evaluated to %lf\n", f); */


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


  /* Close the output file printing diagnostic information */
  fclose(opfp);
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
/*     showAllActivations(); */
    if (usegnuplot) {
      createLayerOpFiles();
      showGnuplot();
    }

    storeActivations(vecnum);
  }
  getZ();
#ifdef dumpArrays
  writeArray(z, "z");
#endif
  
  /* Create the zbar and ztilde arrays by convolution. */

  double_convolve_wrap(z, uMask, ztilde);
  double_convolve_wrap(z, vMask, zbar);

#ifdef dumpArrays
  writeArray(zbar, "zbar");
  writeArray(ztilde, "ztilde");
#endif
  
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

#ifdef dumpArrays
  writeArray(zbarminz, "zbarminz");
  writeArray(ztildeminz, "ztildeminz");
#endif
  
  /* now create dU/dX and dV/dX by convolution. */
  double_convolve_wrap(ztildeminz,
		       uMaskD,
		       dudx);


  double_convolve_wrap(zbarminz,
		       vMaskD,
		       dvdx);
#ifdef dumpArrays
  writeArray(dudx, "dudx");
  writeArray(dvdx, "dvdx");	
#endif
  
  v1 = 1.0 / v;
  u1 = 1.0 / u;

  printf("1/u %lf 1/v %lf \n", u1, v1);

  multArrayInPlace( dudx, u1);
  multArrayInPlace( dvdx, v1);

  subArray(dvdx, dudx, da);
#ifdef dumpArrays
  writeArray(da, "da");
#endif
  
  /* da = dF/dx_a  = 1/v dV/Dx_a - 1/u dU/dx_a */

#ifdef dumpArrays
  printAllActns("allacts");
#endif
  storeTopLayerErrors(da);

  propagateErrors();


  /*
   printAllActns("allacts2");
   */


  createPartials();

}

void createLayerOpFiles()
{
  /* Create the output files for gnuplot to plot */
  /* Will create files layer<n>.op and layer<n>.act
   * for all layers
   */

  /*** Local Variables ***/

  FILE	*actfp;
  FILE	*opfp;  
  int layer;
  char opfile[80], actfile[80];
  int ncells, cellnum, i;
  
  for(layer=0; layer< netInfo.nLayers; layer++) {
    sprintf(opfile, "layer%d.op", layer);
    sprintf(actfile, "layer%d.act", layer);


    opfp = fopen( opfile, "w");
    if (! opfp ) {
      printf("%s: %s could not be opened for writing",
	     __FUNCTION__, opfile);
      exit(-1);
    }

    actfp = fopen( actfile, "w");
    if (! actfp ) {
      printf("%s: %s could not be opened for writing",
	     __FUNCTION__, actfile);
      exit(-1);
    }


    /* print the state of cells in layer LAYER */
    ncells = layerInfo[layer].ncells;
    cellnum = actInfo.startLayer[layer];
    printf("layer %d has %d cells, starting at location %d\n",
	   layer, ncells,  cellnum);
    for(i=0; i< ncells; i++) {
      fprintf(opfp, "%f\n", actInfo.op[cellnum]);
      fprintf(actfp, "%f\n", actInfo.actn[cellnum]);
      cellnum++;
    }
    if (layerInfo[layer].bias == Bias ) {
      cellnum = actInfo.biasIndex[layer];
      printf("bias cell at loc %d\n", cellnum);
      fprintf(opfp, "%f\n", actInfo.op[cellnum]);
      fprintf(actfp, "%f\n", actInfo.actn[cellnum]);
    }      


    fclose(actfp);
    fclose(opfp);
  } /* next layer */
    
  exit(-1);
  
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


void checkNetSize()
{
  /* Check to see if the network is set up ok. There should be a
  ** certain relationshp between the number of input vectors, the
  ** number of cells  in the output layer of the net, and the number
  ** of cells in the virtual output array
  */

  int oplayer, opwid, opht;
  int zwid, zht;
  int nVirtualBlocks;
  
  oplayer = netInfo.nLayers-1;
  opwid = layerInfo[oplayer].ncols;
  opht =  layerInfo[oplayer].nrows;
  printf("real net: op layer is layer %d, %d wide by %d high\n",
	 oplayer, opwid, opht);

  zht = z.ht; zwid = z.wid;
  printf("Virtual array is %d wide by %d high\n", zwid, zht);


  nVirtualBlocks = (zht / opht ) * (z.wid / opwid);
  printf("I think there are %d virtual blocks\n", nVirtualBlocks);

  
  if ( z.ht % opht ) {
    printf("%s: error - z.ht (%d)  must be exactly divisible by oplayer.ht (%d)
\n", __FUNCTION__, z.ht, opht);
    exit(-1);
  }
  
  if ( z.wid % opwid ) {
    printf("%s: error - z.wid (%d) must be exactly divisible by oplayer.wid (%d) \n", __FUNCTION__, z.wid, opwid);
    exit(-1);
  }

  if (nVirtualBlocks != numInputVectors) {
    printf("%s: Error: I think nVirtualBlocks (%d) should equal numInputVectors (%d)\n", __FUNCTION__, nVirtualBlocks, numInputVectors);
    exit(-1);
  }
  
}
  
char finishedFn(Real *wts, int iteration)
{
  /* Decide whether we have finished? */
  char rval;
  Real corrn;

  /* added Sat Dec  9 1995
   * dump out correlations into a file called corrn.dat.
   * This is slightly naughty, as it never closes this file,
   * just keeps flushing it out...
   */
  

#ifdef dumpCorrn
  static FILE *corrnfp = NULL;
#endif

  /* Print out the correlation for the fun of it, although at the
   * moment it is not being used to decide whether to stop running the
   * program.
   */
  /* Decide if we want to output the correlation */
  if (compcorrn) {
    corrn = Rvec_correlate(z.data, shifts.data, 0, numInputVectors );
    printf("Correlation %lf\n", corrn);

#ifdef dumpCorrn
    if (corrnfp == NULL ) {
      /* open up corrn file */
      corrnfp = fopen( "corrn.dat", "w");
      if (! corrnfp ) {
	printf("%s: %s could not be opened for writing",
	       __FUNCTION__, "corrn.dat");
	exit(-1);
      }
    }

    fprintf(corrnfp, "%lf\n", corrn);
    fflush(corrnfp);
    
#endif
  }
  /* maxiterations is a global parameter. */
  
  if (iteration > maxiterations) {
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
  inputSkipX = 2;
  inputSkipY = 0;
  
  numInputVectors = 1000;

  cgmax = 1;			/* By default, we will maximise function. */
  checker = 1;			/* Check code rather than use CG to learn. */

  maxiterations = 100;

  useHalf = 0;
  uhalf = 32;
  vhalf = 320;

  strcpy(initWts, "none");
  strcpy(results, "results");
  doLearning = 1;
  oneImage = 0;
  compcorrn = 1;
  usegnuplot = 0;
}

void printParams()
{
  /* Print the parameters after the parameters file has been loaded */
  printf("*** System Parameters ***\n");
  printf("netFile %s\n", netFile);
  printf("inputFile %s\n", inputFile);

  printf("useHalf = %d\n", useHalf);
  printf("uhalf = %d\n", uhalf);
  printf("vhalf = %d\n", vhalf);
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
  printf("inputSkipX = %d\n", inputSkipX);
  printf("inputSkipY = %d\n", inputSkipY); 
  printf("inputHt = %d\n", inputHt);
  printf("inputWid = %d\n", inputWid);

  printf("cgmax = %d\n", cgmax);
  printf("checker = %d\n", checker);

  printf("maxiterations = %d\n", maxiterations);

  printf("initWts = %s\n", initWts);
  printf("results = %s\n", results); 
  printf("doLearning = %d\n", doLearning);
  printf("oneImage = %d\n", oneImage);
  printf("compcorrn = %d\n", compcorrn);
  printf("usegnuplot = %d\n", usegnuplot); 
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


void checkNetPerformance()
{
  /* Just read in the initial weights and run the input images through
   * the network to see how the network performs.
   *
   * We dont do any learning here.
   */


  /*** Local Variables ***/
  int vecnum;
  double corrn;
  
  printf("Just checking the network's performance\n");
  
  printf("reading in wts from file %s\n", initWts);
  readWts(initWts);
  
  /* Debug: write out the weights again just to check that the wts were
   * read in ok.  */
  writeWts("checkwts");

  /* Now run each input vector through the network... */
  /* This chunk of code was copied from  calcMeritAndPartials(). */
     
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
/*     showAllActivations(); */
/*    showGnuplot(); */

    storeActivations(vecnum);
  }
  getZ();

  /* Finally, now we have the outputs in z, we can dump them out to the file.
   */
  
  writeArray(z, results);

  /* Get the correlation */
  corrn = Rvec_correlate(z.data, shifts.data, 0, numInputVectors );
  printf("Correlation %lf\n", corrn);


  
}

  
/*************************** Version Log ****************************/
/*
 * $Log: testnet.c,v $
 * Revision 1.11  1995/12/11  06:27:37  stephene
 * *** empty log message ***
 *
 * Revision 1.10  1995/12/10  17:55:45  stephene
 * Looking at 2d networks now.
 *
 * Revision 1.9  1995/12/09  16:42:01  stephene
 * *** empty log message ***
 *
 * Revision 1.8  1995/12/08  00:18:00  stephene
 * Inclusion of Rvec_correlate
 *
 * Revision 1.7  1995/12/07  15:43:14  stephene
 * post nips sort out
 *
 * Revision 1.6  1995/11/23  16:20:34  stephene
 * CG now installed
 *
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
