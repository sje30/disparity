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
*** $Revision: 1.2 $
*** $Date: 1995/11/21 23:33:49 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/testconvolve.c,v 1.2 1995/11/21 23:33:49 stephene Exp stephene $";
#endif

/* Simple functions to test the convolve function */


/* -  Include Files - */
#include <stdio.h>
#include <math.h>

#include "convolve.h"
/*#include "gen.h" */
#include "rnd.h"

/* - Defines - */

/* - Function Declarations - */
/* Copied printDoubleArray from mygen.c for portability */
void printDoubleArray(FILE *stream, double *arr, int wid, int ht);
void test2d2();
void test2d();
void test1d();


/* - Global Variables - */ 



/* - Start of Code  - */
#ifdef needmain
main(int argc, char *argv[])
{
  printf("hello world\n");
/*   test1d(); */
  test2d2();
} /* end of main() */
#endif

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



void testConvolve1d()
{
  Array a1, output;
  Mask m;
  Real	*data;
  int wid, ht, i;
  
  wid = 9; ht=1;
  
  createRndArray(wid, ht, &a1);
  createArray(wid, ht, &output);

  m.wid = wid;
  m.ht = 1;

  data = (Real*)calloc(wid, sizeof(Real));
  if (! data) { 
    printf("%s: could not allocate space for data\n", __FUNCTION__);
    exit(-1);
  }
  m.maskExtent =4;
  m.data = data;

  for(i=0; i<wid; i++) {
    a1.data[i] = i;
    m.data[i] = 2 *i;
  }

  double_convolve_wrap(a1, m, output);
  
  printf("My result: \n");
  printDoubleArray(stdout, output.data, wid, 1);

}

void printDoubleArray(FILE *stream, double *arr, int wid, int ht)
{
  /* Array is a 2d floating point array of size wid*ht, stored as a 1d
   * array, row by row. The array is printed to the filepointer (which
   * can be stdout)
   */
  
  /*** Local Variables ***/
  int i;
  while (--ht >= 0)  {
    i=wid;
    while (--i >= 0) {
      fprintf(stream, "%lf ", *(arr++));
    }
    fprintf(stream, "\n");
  }
}


/*************************** Version Log ****************************/
/*
 * $Log: testconvolve.c,v $
 * Revision 1.2  1995/11/21  23:33:49  stephene
 * About to include CG Code
 *
 * Revision 1.1  1995/11/17  21:03:47  stephene
 * Initial revision
 *
 */

