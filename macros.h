#ifndef _macros_
#define _macros_

/* #define MAC */

#define __TCL_DEBUG__ 1

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define FALSE                  0
#define TRUE                   1

#define VERBOSE 		TRUE
/*#define DEGREES*/
#define RADIANS
#define TRACE FALSE

#define    type_PTR							0
#define 	type_REAL						1
#define 	type_INT							2
#define 	type_CHAR 						3
#define 	type_PTR_REAL 				5
#define 	type_PTR_C_CU 				6
#define 	type_PTR_C_OBJECT 		7
#define 	type_PTR_C_UNIT			8

#define 	OR					||
#define	AND					&&

#define jswitch(a) 		thestring=(a); if (FALSE) ;
#define jcase(b)			else if ( SAME_STRING(thestring,b) )
#define jdefault			else

#define DEBUG	FALSE

/* SET UP FLOATING POINT TYPES */

#define 	FLOAT 				FALSE
#define 	DOUBLE 			TRUE
#define	LONG_DOUBLE 	FALSE

#if (FLOAT==TRUE)
#define  REAL 							float
#define  CS_REAL_READ			"%f "
#define  CS_REAL_PRINT			"%.3f "
#define    CS_REAL_FILE_PRINT 		"%.10f "
#define  jf								f
#define	REAL_MAX				FLT_MAX
#define	REAL_MIN					FLT_MIN
#define	REAL_DELTA				FLT_DELTA
/* REAL_DELTA smallest fraction that can be represented */
#endif

#if (DOUBLE==TRUE)
#define 	REAL 						double
#define	CS_REAL_READ			"%lf "
#define	CS_REAL_PRINT		"%.3lf "
#define    CS_REAL_FILE_PRINT 		"%.10lf "
#define    jf								lf
#define	REAL_MAX				DBL_MAX
#define	REAL_MIN					DBL_MIN
#define	REAL_DELTA				DBL_DELTA
#endif

#if (LONG_DOUBLE==TRUE)
#define 	REAL							long double
#define 	CS_REAL_READ 		"%Lf "
#define    CS_REAL_PRINT 		"%.3Lf "
/* CS_REAL_FILE_PRINT with LG gives a precision of 19 sig places */
#define    CS_REAL_FILE_PRINT 	" %.30LG "
#define	jf								Lf
#define	REAL_MAX				LDBL_MAX
#define	REAL_MIN					LDBL_MIN
#define	REAL_DELTA				LDBL_DELTA
#endif

/*****************************************/
/* FILE I/O MACROS */

/* these are macros, so fp is incremented in code */
#define REAL_READ(fp,r) 	fscanf(fp,CS_REAL_READ,&r)
#define INT_READ(fp,i) 		fscanf(fp,"%d",&i)

/*****************************************/

#define MAX_STRING_LENGTH 				512

#ifdef jimdefapply
#define APPLY_FUNCTION1(f,x)  			(*f)(x)
#define APPLY_FUNCTION2(f,x1,x2)  	(*f)((x1),(x2))
#else
#define APPLY_FUNCTION1(f,x)  			f((x))
#define APPLY_FUNCTION2(f,x1,x2)  	f((x1),(x2))
#endif

#define SAME_STRING(s1,s2)			(strcmp((s1),(s2))==0)
#define NOT_SAME_STRING(s1,s2)	( ( ! SAME_STRING( s1, s2 ) ) )


#define EOS '\0'
#define EOL '\n'

#define FAIL 					(-1.0)
#define SUCCESS 				(1.0)

#define NOT(a)                 ((a) ? FALSE : TRUE)
#define EQ						==

#define 	CHECK_FALSE(bool)	if (bool) \
	printf("CHECK_FALSE: WARNING. BOOL IS TRUE AND SHOULD BE FALSE\n")

#define 	CHECK_TRUE(bool)	if NOT(bool) \
	printf("CHECK_TRUE: WARNING. BOOL IS FALSE AND SHOULD BE TRUE\n")
	
#define RESCALE_VALUE(val,old_rmin, old_rmax, rmin,rmax) \
										( ((val)-(old_rmin))*((rmax)-(rmin))/((old_rmax)-(old_rmin)) + (rmin) )

#define CONVERT_INDICES(v,old_imin,old_imax,new_imin,new_imax) \
	if ( (old_imax) - (old_imin) == (new_imax) - (new_imin) ) \
		v += (old_imin) - (new_imin); \
	else 	\
		{printf("ERROR: Indices are incompatible\n"); exit(-1);}

#define ERROR_NULL(ptr)		if (ptr==NULL) \
											{printf("ERROR_NULL: Pointer==NULL\n"); exit(0);}

#define ERROR_TRUE(bool)		if (bool==TRUE) \
											{printf("ERROR_TRUE: bool==TRUE\n"); exit(0);}

#define ERROR_FALSE(bool)		if (bool==FALSE) \
											{printf("ERROR_FALSE: bool==FALSE\n"); exit(0);}
											
#define ERROR_ZERO(bool)		if (bool==0) \
											{printf("ERROR_ZERO: num==0\n"); exit(0);}

#define WARN_NULL(ptr)		if (ptr==NULL) \
											{printf("\n\n\nWARN_NULL: Pointer==NULL\n\n\n"); }
#define WARN_TRUE(bool)		if (bool==TRUE) \
											{printf("\n\n\nWARN_TRUE: bool==TRUE\n\n\n"); }

#define WARN_FALSE(bool)		if (bool==FALSE) \
											{printf("\n\n\nWARN_FALSE: bool==FALSE\n\n\n"); }
											
#define WARN_ZERO(bool)		if (bool==0) \
											{printf("\n\n\nWARN_ZERO: num==0\n\n\n"); }
										
/***************************************************************
*
* 		MATH MACROS
*
****************************************************************/

#include <float.h>

/* USE J BECUASE SOME MACROS ARE DEFINED IN THINK C LIBRARIES */
#define ZERO 					0.0
#define ONE						1.0
#define JINFINITY    			1.0e7
/* defined in float.h */
#define JTINY 					FLT_MIN



/* sd = (2*sd^2) */
#define ODD_GAUSS(x,odd_sd)  exp( (double) - JDIV(JSQR(x),(odd_sd)) )

#define JEXP(x)					exp((double)(x))
#define JLOG(x)					log((double)(x))

#define IRANGE(imin,imax)	( (imax)!=(imin) && (imax)>(imin) ? (imax) - (imin) + 1 : 0 )

#define IS_EVEN(n)				((n)%2==0)
#define IS_ODD(n)				((n)%2!=0)

#define MEAN2(a,b)			( (a+b)/2.0 )

#define JSIGN(x) 				( (x)>0 ? 1.0 : ( (x)<0 ? -1.0 : 0.0 ) )
#define JSQR(x)               ( ((REAL)x)*(x) )
#define JSQRT(x)				( (x)==0.0 ? 0.0 : sqrt((double)(x))  )
#define JPOW(n,e)			(n)==0.0 ? 0.0 : pow((double)(n),(double)(e))	

/* SAFE DIVIDE */
/* NB This may multiply evaluate b. */
#define JDIV(a,b) 	 (b)==0.0 ? \
		( (a)==0.0 ? ZERO : print_error("ERROR:JDIV divide by zero") ) \
		 	:		(a)/(b)
		 	
/*				( (b)==ZERO ? ((a)==ZERO ? ZERO : JINFINITY ) : (a)/(b) )*/

/* INVERSE TRIG FUNCS */
/* acos returns NaN if a=1.0 which is evaluated by C as 0.0,
    but by pop as garbage */

#define     PI              3.14159265
#define	e					2.71828

#define     ONERAD      57.2958

#define RAD2DEG(r)		(ONERAD*r)
#define DEG2RAD(d)		(d/ONERAD)

#ifdef RADIANS
#define JACOS(a)     	acos( (double) (a) )
#define JASIN(a)			asin( (double) (a) )
#define JATAN(a)  		( atan( (double) (a)) )
/* #define JACOS(a)	fabs(a)<1.0 ? acos( (double) (a)) : 0.0 ) */
/*#define JASIN(a)		( fabs(a)<1.0 ? asin( (double) (a)) : PI )*/

#define JTAN(a)			tan( (double) (a) )
#define JSIN(a)			sin( (double) (a) )
#define JCOS(a)			cos( (double) (a) )
#endif

#ifdef DEGREES
#define JACOS(a)     	( ONERAD * ( fabs(a)<1.0 ? acos((double)(a)) : 0.0 ) )
#define JASIN(a)      	( ONERAD * ( fabs(a)<1.0 ? asin((double)(a)) : PI ) )
#define JATAN(a)    		( ONERAD * ( atan((double)(a)) ) )

#define JTAN(a)			tan( (double) DEG2RAD(a) )
#define JSIN(a)			sin( (double) DEG2RAD(a) )
#define JCOS(a)			cos( (double) DEG2RAD(a) )
#endif

#define JMAX(x,y) 								((x)>(y)?(x):(y))
#define JMIN(x,y) 								((x)<(y)?(x):(y))

#define BETWEEN(x,xmin,xmax)    	( ((x) < (xmax)) && ((x) > (xmin)) )
#define INCLUDES(x,xmin,xmax)    	( ((x) <= (xmax)) && ((x) >= (xmin)) )

#define JABS(x)					( (x)>=0.0 ? (x) : -(x) )
#define JROUND(x) ( ceil((double)x)-(x) < (x)-floor((double)x) ? ceil((double)x) : floor((double)x) )

/***************************************************************************/
/*
* PRINT MACROS, nbg with think c pre_processor
*/

#ifndef jfile
#define jprintf							printf
#endif

#define FILE_RPRINT(fp,r)		( r>=0 ? fprintf(fp,"%s"," "), fprintf(fp,CS_REAL_PRINT,r) : fprintf(fp,CS_REAL_PRINT,r) )

#define RPRINT(r)					( r>=0 ? jprintf(" "), jprintf(CS_REAL_PRINT,r) : jprintf(CS_REAL_PRINT,r) )
#define SPRINT(s)					jprintf("%s ",s);
#define IPRINT(i)						( i>=0 ? jprintf(" %d ",i) : jprintf("%d ",i) )


#define NEWLINE					jprintf("\n")
#define FILE_NEWLINE(fp)	fprintf(fp,"\n")

/* #define SPR(str)				jprintf("str=%s\n",str) */
/*#define SSPR(str)				jprintf(str); */

#define PPRINT(str,val) if (val>=0) {jprintf(" "); jprintf(str,val);} else jprintf(str,val)

#if (TRACE==TRUE)
#define TRACE_SSPR(str) jprintf("%s\n",str);
#else
#define TRACE_SSPR(str) /* do nothing */
#endif

#define IPR1(x1) jprintf("x1=%d\n",x1)
#define IPR2(x1,x2) jprintf("x1=%d x2=%d\n",x1,x2)
#define IPR3(one,two,three) jprintf("one=%d two=%d three=%d\n",one,two,three)
#define IPR4(one,two,three,four) \
    jprintf("one=%d two=%d three=%d four=%d\n",one,two,three,four)

#define FPR1(one) \
    jprintf("one=%f\n",one)
#define FPR2(one,two) \
    jprintf("one=%.8f two=%.8f\n",one,two)
#define FPR3(one,two,three) \
    jprintf("one=%.3f two=%.3f three=%f \n",one,two,three);
#define FPR4(one,two,three,four) \
    jprintf("one=%.3f two=%.3f three=%.3f four=%f\n",one,two,three,four)
#define FPR5(one,two,three,four,five) \
    jprintf("one=%.3f two=%.3f three=%.3f four=%.3f five=%f\n",one,two,three,four,five)
#define FPR6(one,two,three,four,five,six) \
    jprintf("one=%.3f two=%.3f three=%.3f four=%.3f five=%.3f six=%f\n",one,two,three,four,five,six)
#define FPR7(one,two,three,four,five,six,seven) \
    jprintf("one=%.3f two=%.3f three=%.3f four=%.3f five=%.3f six=%.3f seven=%f\n", one,two,three,four,five,six,seven)
#define FPR8(one,two,three,four,five,six,seven,eight) \
    jprintf("one=%.3f two=%.3f three=%.3f four=%.3f five=%.3f six=%.3f seven=%.3f eight=%f\n", one,two,three,four,five,six,seven,eight)

#define ILOOP							for ((i)=0;(i)<(n);(i)++)
#define LOOP(i,imin,imax)       for ((i)=(imin);(i)<(imax);(i)++)
#define NEGLOOP(i,imax)         for ((i)=(-imax);(i)<=(imax);(i)++)

#define PRINT_CONST(b)          jprintf("b=%d\n",b);
#define PRINT_FLAG(b)           if (b!=0) jprintf("b\n");



#endif
