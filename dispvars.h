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
*** $Revision: 1.4 $
*** $Date: 1995/12/07 15:42:57 $
****************************************************************************/


#ifndef _DISPVARS_H
#define _DISPVARS_H


/* Global Disparity Variables used by the lex file. */
char netFile[100];		/* Name of the file that stores the network
				 * structure information.
				 * Default: default.net
				 */
char inputFile[100];		/* Name of the file where inputs can
				 * be read in from
				 */

double ulambda;			/* Lambda for short range mean U */
double vlambda;			/* Lambda for long range mean V */


int outputWid; 			/* Width of the output cells. See also
				 * outputHt. */

int outputHt;			/* Height of the output cells. For a
				 * 1d arrangement of output cells,
				 * this should be set to 1.
				 */

int totalInputWid;		/* Total width of Jim's input files.
				 * Default to 7000. */
int totalInputHt;		/* Total height of Jim's input files.
				 * Default to 5. */
char  image1File[100]; 		/* File for left image */
char  image2File[100]; 		/* File for right image */
char  shiftsFile[100];		/* File for shifts data */
int numInputVectors;		/* number of input vectors to make
				   from the total input files. */

int inputSkip;			/* How many columns to skip between
				   successive input vectors */

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

int useHalf;			 /* If nonzero, use uhalf and vhalf to
				  * specify the half lives for U and
				  * V.  Otherwise, use ulambda and
				  * vlambda.  Default is 0. */
int uhalf, vhalf; 		/* Half lives of U and V */

#endif


/*************************** Version Log ****************************/
/* $Log: dispvars.h,v $
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
