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
*** $Revision: 1.1 $
*** $Date: 1995/11/09 20:49:36 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/testnet.c,v 1.1 1995/11/09 20:49:36 stephene Exp stephene $";
#endif


/* Test file for testing the network */

/* -  Include Files - */
#include "dispnet.h"
#include "dispvars.h"
#include "readnet.h"
#include <stdio.h>

/* - Defines - */

/* - Function Declarations - */
static void setParamDefaults();
static void printParams();


/* - Global Variables - */ 



/* - Start of Code  - */
main(int argc, char *argv[])
{

  extern int	yylex();

  if (argc !=2 ) {
    printf("Usage %s <parameter file> \n", argv[0]);
    exit(-1);
  }

  setParamDefaults();
  freopen( argv[1], "r", stdin);
  
  printf("About to parse %s\n", argv[1]);
  while ( yylex() ) ;

  printf("End of lex\n");

  printParams();

  readNet(netFile);
} /* end of main() */


void setParamDefaults()
{
  /* Set up the defaults for the parameters. */
  strcpy(netFile, "default.net");
  strcpy(inputFile, "inputs");
}

void printParams()
{
  /* Print the parameters after the parameters file has been loaded */
  printf("*** System Parameters ***\n");
  printf("netFile %s\n", netFile);
  printf("inputFile %s\n", inputFile);
}

     
/*************************** Version Log ****************************/
/*
 * $Log: testnet.c,v $
 * Revision 1.1  1995/11/09  20:49:36  stephene
 * Initial revision
 *
 */

