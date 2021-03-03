/****************************************************************************
***
***
*** dispcorrn.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 17 Nov 95
***
*** $Revision: 1.1 $
*** $Date: 1995/11/17 12:07:30 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/dispcorrn.c,v 1.1 1995/11/17 12:07:30 stephene Exp $";
#endif


/*****************************************************************************/
/* Determine the correlation between the shift vectors and the cell outputs. */
/*****************************************************************************/

/* -  Include Files - */

/* - Defines - */

/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */


Real
Rvec_correlate(Real *x, Real *y, int imin, int imax)
{
	/* *x and *y are indexed from imin to imax */
        /* modified from NUM REC code */
        
        int             j;
        Real    yt,xt,t,df;
        Real    syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
        Real            den;
        int             n;
        Real            prob, r, z;
        
        n = imax-imin+1;

        for (j=1;j<=n;j++) 
                {
                ax += x[j];
                ay += y[j];
                }
                
        ax /= n;
        ay /= n;
        for (j=1;j<=n;j++) {
                xt=x[j]-ax;
                yt=y[j]-ay;
                sxx += xt*xt;
                syy += yt*yt;
                sxy += xt*yt;
        }
        den = JSQRT(sxx*syy);
        r =  sxy / den;
        
        return((Real)r);
}


/*************************** Version Log ****************************/
/*
 * $Log: dispcorrn.c,v $
 * Revision 1.1  1995/11/17  12:07:30  stephene
 * Initial revision
 *
 */

