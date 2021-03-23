
/* extra code needed to use Jim's code Wed Nov 22 1995 */
#ifndef REAL
#define REAL double
#endif

/* #include "cg_williams_module.h" */


double *
Rvec_create();

int
bp_check_func_deriv(
		    REAL *w_orig,	/* INITIAL VALUES OF PARAMETERS */
		    int imin, 
		    int imax,		/* MIN AND MAX INDICES OF PARAMS IN W_ORIG */
		    REAL (*func)(REAL *w),
		    void (*dfunc)(REAL *w,REAL *g),
		    REAL tol,
		    int max_num_line_searches
		    );

int
finished(REAL *wts, int iter);


#ifdef MAC
void
set_graph_title(char *str);
void
plot(void);
#endif
