/****************************************************************************
***
***
*** $RCSfile$
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 15 Jan 96
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef _DISPDEFINES_H
#define _DISPDEFINES_H

/* These are some defines that are quite important, and will probably
** want to be changed by users. */

/* Change #define to #undef if you want to undefine the statements,
** rather than deleting them from this file.
*/

#define MAX_NUM_WEIGHTS 50000	/* The maximum number of weighta that
   				** will be allocated by the
   				** program. If you exceed this value,
				** you will have to recompile the program.
				*/

#undef dumpArrays		/* Dump out some important files, such
				** as z, zbar, ztilde, da, allacts
				** - this will slow things down a bit though
				*/
#define seedrndtime()  srand48((long)time(NULL))
#define seedrnd(x) srand48((long)x)
#endif


/*************************** Version Log ****************************/
/* $Log$
 *
 */
