/****************************************************************************
***
*** Time-stamp: <>
***
*** testnet.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 09 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header$";
#endif


/* Test file for testing the network */

/* -  Include Files - */
#include "dispnet.h"
#include "readnet.h"
/* - Defines - */

/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */
main(int argc, char *argv[])
{
  readNet("3layer.net");
} /* end of main() */


/*************************** Version Log ****************************/
/*
 * $Log$
 */

