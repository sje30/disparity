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
*** $Revision: 1.10 $
*** $Date: 1997/06/12 14:05:53 $
****************************************************************************/


#ifndef _DISPVARS_H
#define _DISPVARS_H


/* Global Disparity Variables used by the lex file. */
extern char netFile[100];		/* Name of the file that stores the network
				 * structure information.
				 * Default: default.net
				 */
extern char inputFile[100];		/* Name of the file where inputs can
				 * be read in from.
				 */

extern double ulambda;			/* Lambda for short range mean U */
extern double vlambda;			/* Lambda for long range mean V */
				/* The width of the means for U and V
                                 * can either be set using lambda, or
                                 * by providing the half lives
                                 * directly. See useHalf. uhalf and vhalf. */


extern int outputWid; 			/* Width of the output cells. See also
				 * outputHt. */

extern int outputHt;			/* Height of the output cells. For a
				 * 1d arrangement of output cells,
				 * this should be set to 1.  If this
				 * value is not one, it is assumed
				 * that the network is 2d.  */

extern int totalInputWid;		/* Total width of Jim's input files.
				 * Default to 7000. */
extern int totalInputHt;		/* Total height of Jim's input files.
				 * Default to 5. */
extern char  image1File[100]; 		/* File for left image */
extern char  image2File[100]; 		/* File for right image */
extern char  shiftsFile[100];		/* File for shifts data */
extern int numInputVectors;		/* number of input vectors to make
				   from the total input files. */

extern int inputSkipX;			/* How many columns to skip between
				 * successive input vectors. Default 2.
				 */
extern int inputSkipY; 		/* How many rows to skip between input
				 * vectors from a 2d image. Default 0.
				 */

extern int inputHt; 			/* Height of each input vector from
				   one image */
extern int inputWid; 			/* Width of each input vector from one im. */

extern int cgmax;			/* If non zero, then maximise merit
				 * function in conjugate gradient
				 * method. Default is 1.
				 */
extern int checker;			/* Do we use conjugate gradient or
				 * just check the derivatives?  Non
				 * zero means that we check.
				 * Default: 1.
				 */

extern int maxiterations;		/* Maximum number of iterations for CG
				 * to perform.  Default is 100. */

extern int useHalf;		        /* If nonzero, use uhalf and vhalf to
				 * specify the half lives for U and
				 * V.  Otherwise, use ulambda and
				 * vlambda.  Default is 0. */
extern int uhalf, vhalf; 		/* Half lives of U and V. If useHalf
				 * is non zero, then these values are
				 * used to specify the size of U and
				 * V. */

extern char  initWts[100];		/* Name of file storing the initial
				 * weights of the network if any.  If
				 * value = "none" (default value),
				 * then rnd initial weights are
				 * created. */

extern char  results[100];		/* Where to put the results of
				 * checking the nets
				 * performanace. Default "results". */


extern int doLearning;			/* If non zero (default), then do
				 * learning. Otherwise, just read in
				 * initial weights and run input
				 * images through the network, and
				 * calculate the correlation. */

extern int oneImage;			/* Do we want to read in one image or two?
				 * Default 0 - means we have more than one
				 * image file. */
extern int compcorrn;			/* Do we want to compute correlation?
				 * If non zero, correlation is computed.
				 * Default: 1*/
extern int usegnuplot;			/* Should we use gnuplot to display
				 * network activity? If non zero, then
				 * use gnuplot. Default: 0.  */
extern int seed;			/* If non-zero, use this value as the seed.
				 * Otherwise, clock used to set seed.
				 * Default 0.
				 */
extern int normInput;			/* Non-zero if we want to normalize each
				 * image-patch on way in. Default 0.
				 */
extern int noshifting;			/* Non-zero if we don't need to shift the
				 * desired output values when computing
				 * the correlation between actual and
				 * desired output.
				 */
#endif


/*************************** Version Log ****************************/
/* $Log: dispvars.h,v $
 * Revision 1.10  1997/06/12  14:05:53  stephene
 * added seed for controlling random numbers.
 *
 * Revision 1.9  1996/01/16  01:27:40  stephene
 * now have the code in place so that weight sharing can be done or left
 * out.
 *
 * Revision 1.8  1995/12/13  04:04:04  stephene
 * *** empty log message ***
 *
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
