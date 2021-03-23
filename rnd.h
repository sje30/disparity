/* rnd.h	-- Stephen Eglen
 *
 * $Header: /rsuna/home2/stephene/Clib/rnd.h,v 1.2 1994/10/17 21:53:08 stephene Exp stephene $
 * $Log: rnd.h,v $
 * Revision 1.2  1994/10/17  21:53:08  stephene
 * about to update to use random and srandom instead of rand and srand.
 *
 * Revision 1.1  1994/09/27  14:56:30  stephene
 * Initial revision
 *
 */



#ifndef _RND_H
#define _RND_H

/* Random number generator. */
/* Usage : call rnd() to get a random number as a float between 0 and 1.
 * Seeding it:
 * Call seedrnd(x) to seed it with X.
 * Call seedrndtime() to seed it with the system clock */


/*** Include Files ***/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h> /* see man -s 2 time */
#include <sys/types.h>

/*** Function Declarations ***/
time_t time(time_t *tloc);

void srand48(long s);
double drand48(void);

#define rnd() (float)(drand48())
#define seedrnd(x) srand48((long)x)
#define seedrndtime()  srand48((long)time(NULL))

/* Double version */
#define drnd() (float)(drand48())

#endif
