/****************************************************************************
***
***
*** testconvolve.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 17 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header$";
#endif

/* Simple functions to test the convolve function */


/* -  Include Files - */
#include <stdio.h>
#include <math.h>

#include "convolve.h"
#include "gen.h"

/* - Defines - */

/* - Function Declarations - */
void test2d2();
void test2d();
void test1d();


/* - Global Variables - */ 



/* - Start of Code  - */
main(int argc, char *argv[])
{
  printf("hello world\n");
/*   test1d(); */
  test2d2();
} /* end of main() */


void test1d()
{
  double input[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double  mask[] = {-1.0, -2.0, 5.0, 3.0, 4.0};
  double output[6];
  int inputWid = 6;
  int maskwid = 5;
  int maskextent = 2;

  double_convolve1d_wrap(input,  inputWid,
			 mask, maskwid, maskextent,
			 output);

  printDoubleArray(stdout, output, inputWid, 1);
}


void test2d()
{
  double input[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  double  mask25[] = {-1.0, -2.0, 5.0, 3.0, 4.0, 1.0, -1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, -1.0, 1.0, 0.0, 0.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0};	
  
  double  mask9[] = {-1.0, -2.0, 5.0, 3.0, 4.0, 1.0, -1.0, 1.0, 2.0};
  double output[12];
  int inputWid = 4;
  int inputHt = 3;
  int maskwid = 3;
  int maskht = 3;
  int maskextent = 1;
  
  double_convolve2d_wrap(input,  inputWid,  inputHt,
			  mask9,  maskwid, maskht,  maskextent,
			  output);

  printDoubleArray(stdout, output, inputWid, inputHt);
}


void test2d2()
{
  double input[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  double  mask25[] = {-1.0, -2.0, 5.0, 3.0, 4.0, 1.0, -1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, -1.0, 1.0, 0.0, 0.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0};	
  double mask1[] = {1.0};
  double  mask9[] = {0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11};
  double output[12];
  int inputWid = 4;
  int inputHt = 3;
  int maskwid = 1;
  int maskht = 1;
  int maskextent = 0;
  
  double_convolve2d_wrap(input,  inputWid,  inputHt,
			  mask1,  maskwid, maskht,  maskextent,
			  output);

  printDoubleArray(stdout, output, inputWid, inputHt);
}

  
/*************************** Version Log ****************************/
/*
 * $Log$
 */

