/* 
 * This proc does normal bp and also checks that the func deriv really
 * *is the derivative of the function by taking small steps and
 * assuming *that surface is locally planar. It can also be used for
 * normal BP by *making step larger (actually step-size adjusts itself
 * *automatically).
 */

/*
 * Changes made: Stephen Eglen, Wed Nov 22 1995
 *
 * Occurences of "%Lf" have been changed to "%lf" - wasnt working with "L".
 *
 * All arrays indexed to go from imin to imax, instead of assuming
 * that they go from 1 to len.
 *
 * Definition of Rvec_dot put at top of file, so that the Rvec_dot
 * returns a double rather than an int.
 *
 * Definitions of APPLY_FUNCION changed - dont need to use * in front
 * of function name I think.
 */

 
#define ANSI /* #define MAC */

#include <math.h>
#include <stdio.h>


#include "macros.h"
#ifdef jimscodenotused

#include <typedefs.h>

#include <vec_mat.h>
#include <jutil.h>
#include <interrupt.h>
#include <random.h>

#endif

extern double
Rvec_dot(double *v0,double *v1,int imin,int imax);

#include "bp_check_deriv.h"
/* #include <cg_williams.h> */

#ifdef MAC
#include <interrupt.h>
#include <graph_plot.h>
#define JGRAPHICS 1
#endif

#ifndef MAC
int
plot();
#endif

int				CG_TRY_NUM;
REAL			BP_CHECK_STEP_FACTOR = 1e-6;
char			STEPPED_TO_MIN;

/* FIXED_STEP_FACTOR is normally TRUE to do a fixed small step. */
/* If set to false, it seems to take the ratio away from 1.0 */

char				FIXED_STEP_FACTOR = TRUE;
int					*TRAINING_ITERATION;

int
bp_check_func_deriv(
					REAL *w_orig, 				/* INITIAL VALUES OF PARAMETERS */
					int imin, 
					int imax, 						/* MIN AND MAX INDICES OF PARAMS IN W_ORIG */
					REAL (*func)(REAL *w),
					void (*dfunc)(REAL *w,REAL *g),
					REAL tol,
					int max_num_line_searches
					)
{
    /*--------------------- 
    TRAINING PARAMETERS 
    -----------------------*/
    REAL                eta;
    REAL                alpha;
    REAL		kappa; 
   	
    /*---------- 
    VECTORS 
    -----------*/
    REAL                *w0, *w1;       		/* CURRENT AND NEW WEIGHT VEC */
    REAL                *g0, *g1;       		/* DIRECTION OF STEEPEST DESCENT */
    REAL                *w_try, *g_try; 	/* w_try USED TO EST 2ND DERIV */
    REAL                *s;             			/* SEARCH DIRECTION */
    REAL                *w_change;
    REAL		*temp_vec;
    REAL		*min_vec;
	
    /* MISC */
    REAL            	E0, E1;
    REAL               	temp;
    REAL                *ptr;
    REAL		E_last, E_lastb1, E_lastb2; /* used in print_profile */
    REAL		E_try;
    REAL		dE_dw, temp_eta, temp_alpha;
    REAL		ratio=1.0;
    
    int                	success=TRUE;
    int                	len;
    int                	i, iter=0;
    int                	itemp, num_contiguous_fails;
    int 		counter=0, fin=FALSE;
    
    char		print_profile=FALSE;
    char		too_linear;
    
    REAL z;
    
    len = imax-imin+1;
/*
    w0  = Rvec_create(1,len);    
    w_try = Rvec_create(1,len);
    w_change = Rvec_create(1,len);
*/
    

    w0  = Rvec_create(imin, imax);
    w_try = Rvec_create(imin, imax);
    w_change = Rvec_create(imin, imax);

/*     g0 = Rvec_create(1,len); */
/*     s  = Rvec_create(1,len); */


    g0 = Rvec_create(imin, imax);
    s  = Rvec_create(imin, imax);
    
    
/*     temp_vec  	= Rvec_create(1,len);*/
    temp_vec  	= Rvec_create(imin, imax);
    min_vec 	= Rvec_create(-2,2);

    /* COPY ORIG WTS INTO W0 */
    Rvec_copy(w0,w_orig,imin,imax);
    
    TRAINING_ITERATION = &iter;
    CG_TRY_NUM=1;
    
    /*APPLY_FUNCTION2(dfunc,w0,g0);
      Rvec_plot(g0,imin,imax);
      kappa 	= Rvec_dot(g0,g0,1,len); 
      E0 		= APPLY_FUNCTION1(func,w0);
      printf("Initial E=%.3lf\n",E0); 
      plot();*/
    /* sje mod two lines above in last printf */
    
    while ( NOT(finished(w0,iter)) )
      {
        check_interrupt();	++iter;
        printf("success=%d\n",success);
        CG_TRY_NUM = (success ? iter : CG_TRY_NUM);
        
        if (success)
	  {
	    STEPPED_TO_MIN = TRUE; 
	    APPLY_FUNCTION2(dfunc,w0,g0);
	    E0 = APPLY_FUNCTION1(func,w0);
	    STEPPED_TO_MIN = FALSE;  
	    set_graph_title("Weight partials vector");
	    Rvec_plot(g0,imin,imax);

	    /* This assumes we are going from 1 to len. */
	    /*Rvec_negate(s,g0,1,len); */

	    Rvec_negate(s,g0,imin, imax);

	    /* kappa = Rvec_dot(s,s,1, len); */
	    kappa = 0.0;
	    kappa = (double)Rvec_dot(s,s,imin, imax);
	    z = (double)Rvec_dot(s,s,imin, imax);
/*	    printf("Kappa is %lf z is %lf\n", kappa, z); */
	    num_contiguous_fails=0;
	  } 
		
	printf("\nLine search number = %d E=%.6lf kappa=%lG\n",iter,E0,kappa); /* sje mod */
        /*plot();*/
        
	/* find step size that pushes slope of scale of stepsize just
	 * outside linear region */
	too_linear = TRUE;
	while ( too_linear )
	  {
	    eta = BP_CHECK_STEP_FACTOR / JSQRT(kappa); 
/*	    Rvec_linear_comb(w_try,1.0,w0,eta,s,1,len); */
	    Rvec_linear_comb(w_try,1.0,w0,eta,s,imin, imax);
	    
	    E_try = APPLY_FUNCTION1(func,w_try);
/* 	    Rvec_scale(w_change,s,eta,1,len); */

	    Rvec_scale(w_change,s,eta,imin, imax);

	    
/* 	    ratio = (E0-E_try) / (Rvec_dot(w_change,s,1,len));*/
	    ratio = (E0-E_try) / (Rvec_dot(w_change,s,imin, imax));
	    
	    if NOT(FIXED_STEP_FACTOR)
	      {
		if (ratio < 0.1) 	/* stepped across valley floor! E increased */ 
		  BP_CHECK_STEP_FACTOR /= 2;
		else
		  if ( ratio > 0.66 AND ratio < 1.5 ) /* almost linear region of surface */
		    {
		      BP_CHECK_STEP_FACTOR *= 2;
		      too_linear = TRUE;
		    }
		  else 
		    too_linear = FALSE;
		/* slope is +ve and outside range 0.66->1.5, ie non-linear  */
     			
	      }  /* FIXED_STEP_FACTOR  */
     	
	    /* BP_CHECK_STEP_FACTOR = JMAX(BP_CHECK_STEP_FACTOR,0.001); */
	    printf("BP_CHECK_STEP_FACTOR=%lG\neta=%lG\nActual change in E = %e\nPredicted change in E = %e\nRatio = %lf\n", /* sje mod */
		   BP_CHECK_STEP_FACTOR,
		   eta,
		   E0-E_try,	
		   /* 		   Rvec_dot(w_change,s,1,len),  */
		   Rvec_dot(w_change,s,imin, imax),
		   ratio);
	    if (FIXED_STEP_FACTOR) too_linear=FALSE;
	  } /* too_linear  */
     	
	if (print_profile OR !success)	
	  {
	    alpha = eta;
	    printf("alpha = %.6lG\n",alpha); /* sje mod */
	    itemp = 2;
	    E_last=0.0;
	    for (i=(-itemp);i<=itemp;i++)
	      {
            	temp_alpha = (2.0/itemp)*i*alpha;
            	
            	if 		(i EQ 0) temp=E0; 
            	else if (i EQ 1) temp=E_try;
            	else 
		  {
/* 		    Rvec_linear_comb(temp_vec,1.0,w0,temp_alpha,s,1,len); */
		    Rvec_linear_comb(temp_vec,1.0,w0,temp_alpha,s,imin, imax);
		    temp = APPLY_FUNCTION1(func,temp_vec);
		  }
            	
            	printf("E = %.6lG   	Alpha_factor = %.6lG\n",temp,temp_alpha/alpha); /* sje mod */
            	min_vec[i] = temp;
            	E_lastb2 = E_lastb1;
            	E_lastb1 = E_last;
            	E_last = temp;
	      }
            set_graph_title("Cross section of cost function");
            Rvec_plot(min_vec,-2,2);
            wait();
	  }
        
        success = (E_try<=E0 ? TRUE : FALSE );
        	
        if NOT(success)
	  {
            printf("WARNING: E_try>E0: Value of func after stepping eta=%.3e in search dir = %le\n", /* sje mod */
		   eta,APPLY_FUNCTION1(func,w_try));
            ++num_contiguous_fails;
	  }
      	else
	  {
	    E0 = E_try;
	    ptr = w0; w0 = w_try; w_try = ptr;
	  }
        
        if (num_contiguous_fails>4)
	  success = TRUE;
      } /* iter */
    
    printf("bp_for_IMAX terminating: iter=%d fin=%d\n",iter,fin);
     
    Rvec_copy(w_orig,w0,imin,imax);		
    
    Rvec_destroy(w0,1,len);
    Rvec_destroy(w_try,1,len);
    Rvec_destroy(g0,1,len);
    Rvec_destroy(s,1,len);
    Rvec_destroy(temp_vec,1,len);
    
    return(fin);  /* FINISH IS CODE FOR REASON FOR TERMINATING */
}

int
finished(REAL *wts, int iter)
{
	int	fin = 0;
	if (iter>100)
		fin = 1;
		
	return(fin);
}


/* extra code Wed Nov 22 1995 SJE */
set_graph_title(char *str)
{
  ; /* do nothing */
}

Rvec_plot(REAL *g0, int imin,int imax)
{
  ; /* do nothing */
}

/* Taken from jims/MISC/BP_CG/MAC/GENERAL/VEC_MAT/vec_mat.c */

int
Rvec_scale(REAL *v0, REAL *v1, REAL fac,int imin,int imax)
{
    	register int counter = imax-imin+1;
    	register REAL rg_fac, *rg_v0, *rg_v1;
	check_length(imin,imax);
	v0 += imin;
	v1 += imin;
	
	rg_fac = fac;
    rg_v0 = v0;
    rg_v1 = v1;
    
    while (counter--)
        *v0++ = *v1++ * fac;
}


check_interrupt()
{
  ; /* do nothing */
}
#ifndef MAC


int
check_length(int imin, int imax)
{
	int len;
	len = imax-imin+1;
	
	if (DEBUG)
	/*if (imax-imin<=0)
		printf("check_length: Warning: Length %d imin %d imax %d\n",len,imin,imax);
	*/
	return(len);
}


#endif

/****************************************************************************
***
***
*** bp_check_deriv.c
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 23 Nov 95
***
*** $Revision: 1.2 $
*** $Date: 1995/12/08 00:16:23 $
****************************************************************************/


#ifndef lint
static char *rcsid = "$Header: /rsuna/home2/stephene/disparity/bp_check_deriv.c,v 1.2 1995/12/08 00:16:23 stephene Exp $";
#endif



/* -  Include Files - */

/* - Defines - */

/* - Function Declarations - */


/* - Global Variables - */ 



/* - Start of Code  - */


/*************************** Version Log ****************************/
/*
 * $Log: bp_check_deriv.c,v $
 * Revision 1.2  1995/12/08  00:16:23  stephene
 * *** empty log message ***
 *
 * Revision 1.1  1995/11/23  16:30:22  stephene
 * Initial revision
 *
 */

