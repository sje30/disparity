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
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header$";
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
 * $Log$
 */

