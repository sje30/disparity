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
*** $Revision: 1.1 $
*** $Date: 1995/11/10 21:59:00 $
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

#endif


/*************************** Version Log ****************************/
/* $Log: dispvars.h,v $
 * Revision 1.1  1995/11/10  21:59:00  stephene
 * Initial revision
 *
 *
 */
