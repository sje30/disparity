/****************************************************************************
***
*** Time-stamp: <10 Nov 95 21:58:13 stephene>
***
*** dispvars.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 10 Nov 95
***
*** $Revision: 1.7 $
*** $Date: 1995/12/11 06:26:04 $
****************************************************************************/


#ifndef _DISPVARS_H
#define _DISPVARS_H


/* Global Disparity Variables used by the lex file. */
char netFile[100];		/* Name of the file that stores the network
				 * structure information.
				 * Default: default.net
				 */
char inputFile[100];		/* Name of the file where inputs can
				 * be read in from.
				 */

double ulambda;			/* Lambda for short range mean U */
double vlambda;			/* Lambda for long range mean V */
				/* The width of the means for U and V
                                 * can either be set using lambda, or
                                 * by providing the half lives
                                 * directly. See useHalf. uhalf and vhalf. */


int outputWid; 			/* Width of the output cells. See also
				 * outputHt. */

int outputHt;			/* Height of the output cells. For a
				 * 1d arrangement of output cells,
				 * this should be set to 1.  If this
				 * value is not one, it is assumed
				 * that the network is 2d.  */

int totalInputWid;		/* Total width of Jim's input files.
				 * Default to 7000. */
int totalInputHt;		/* Total height of Jim's input files.
				 * Default to 5. */
char  image1File[100]; 		/* File for left image */
char  image2File[100]; 		/* File for right image */
char  shiftsFile[100];		/* File for shifts data */
int numInputVectors;		/* number of input vectors to make
				   from the total input files. */

int inputSkipX;			/* How many columns to skip between
				 * successive input vectors. Default 2.
				 */
int inputSkipY; 		/* How many rows to skip between input
				 * vectors from a 2d image. Default 0.
				 */

int inputHt; 			/* Height of each input vector from
				   one image */
int inputWid; 			/* Width of each input vector from one im. */

int cgmax;			/* If non zero, then maximise merit
				 * function in conjugate gradient
				 * method. Default is 1.
				 */
int checker;			/* Do we use conjugate gradient or
				 * just check the derivatives?  Non
				 * zero means that we check.
				 * Default: 1.
				 */

int maxiterations;		/* Maximum number of iterations for CG
				 * to perform.  Default is 100. */

int useHalf;		        /* If nonzero, use uhalf and vhalf to
				 * specify the half lives for U and
				 * V.  Otherwise, use ulambda and
				 * vlambda.  Default is 0. */
int uhalf, vhalf; 		/* Half lives of U and V. If useHalf
				 * is non zero, then these values are
				 * used to specify the size of U and
				 * V. */

char  initWts[100];		/* Name of file storing the initial
				 * weights of the network if any.  If
				 * value = "none" (default value),
				 * then rnd initial weights are
				 * created. */

char  results[100];		/* Where to put the results of
				 * checking the nets
				 * performanace. Default "results". */


int doLearning;			/* If non zero (default), then do
				 * learning. Otherwise, just read in
				 * initial weights and run input
				 * images through the network, and
				 * calculate the correlation. */


#endif


/*************************** Version Log ****************************/
/* $Log: dispvars.h,v $
 * Revision 1.7  1995/12/11  06:26:04  stephene
 * new params: initWts, doLearning, results
 *
 * Revision 1.6  1995/12/10  17:55:00  stephene
 * skipInput changed to skipInputX and skipInputY
 *
 * Revision 1.5  1995/12/08  20:17:49  stephene
 * Adding the code to choose between creating 1d and 2d masks in
 * setUpNetwork
 * Adding also the new params useHalf, uhalf and vhalf
 *
 * Revision 1.4  1995/12/07  15:42:57  stephene
 * post nips sort out
 *
 * Revision 1.3  1995/11/23  16:20:04  stephene
 * CG now installed
 *
 * Revision 1.2  1995/11/21  23:33:01  stephene
 * About to include CG Code
 *
 * Revision 1.1  1995/11/10  21:59:00  stephene
 * Initial revision
 *
 *
 */
