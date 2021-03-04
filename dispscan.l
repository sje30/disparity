%{
/****************************************************************************
***
*** dispscan.l
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 10 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/
%}

%{
/* Lex code for disparity program. */
/* Load in the parameter variables */
#include "dispvars.h"
%}
%{
void static	skipcmnts();
%}

	/* Comments can go into this section as long as the line begins 
	 * with whitespace. */

	/* These are definitions for later use.  
	 * S - White space.
	 * W - Word - anything delimited by space. 
	 * F - Float - floating point number ? indicates minus is optional 
	 * I - Int - integer, pos or negative */

S	[ \t\n]*
W	[^ \t\n]+
F	[-]?[0-9\.]+
I	[-]?[0-9]+
C	^"\#"[^\n]*
%%
"#"[^\n]*		; /*printf("skip comment to end of line\n");*/
"/*"			skipcmnts();
"netFile"{S}{W}		sscanf(yytext,"netFile%s", netFile);
"inputFile"{S}{W}	sscanf(yytext,"inputFile%s", &inputFile);

"useHalf"{S}{I}		sscanf(yytext,"useHalf %d", &useHalf);
"ulambda"{S}{F}		sscanf(yytext,"ulambda %lf", &ulambda);
"vlambda"{S}{F}		sscanf(yytext,"vlambda %lf", &vlambda);
"uhalf"{S}{I}		sscanf(yytext,"uhalf %d", &uhalf);
"vhalf"{S}{I}		sscanf(yytext,"vhalf %d", &vhalf);

"outputWid"{S}{I}	sscanf(yytext,"outputWid %d", &outputWid);
"outputHt"{S}{I}	sscanf(yytext,"outputHt %d", &outputHt);
"totalInputWid"{S}{I}	sscanf(yytext,"totalInputWid %d", &totalInputWid);
"totalInputHt"{S}{I}	sscanf(yytext,"totalInputHt %d", &totalInputHt);

"image1File"{S}{W}	sscanf(yytext,"image1File %s", image1File);
"image2File"{S}{W}	sscanf(yytext,"image2File %s", image2File);
"shiftsFile"{S}{W}	sscanf(yytext,"shiftsFile %s", shiftsFile);

"inputWid"{S}{I}	sscanf(yytext,"inputWid %d", &inputWid);
"inputHt"{S}{I}	sscanf(yytext,"inputHt %d", &inputHt);
"inputSkipX"{S}{I}	sscanf(yytext,"inputSkipX %d", &inputSkipX);
"inputSkipY"{S}{I}	sscanf(yytext,"inputSkipY %d", &inputSkipY);

"numInputVectors"{S}{I}	sscanf(yytext,"numInputVectors %d", &numInputVectors);
"cgmax"{S}{I}		sscanf(yytext,"cgmax %d", &cgmax);
"checker"{S}{I}		sscanf(yytext,"checker %d", &checker);
"maxiterations"{S}{I}	sscanf(yytext,"maxiterations %d", &maxiterations);

"initWts"{S}{W}		sscanf(yytext,"initWts %s", initWts);
"results"{S}{W}		sscanf(yytext,"results %s", results);
"doLearning"{S}{I}	sscanf(yytext,"doLearning %d", &doLearning);

"oneImage"{S}{I}	sscanf(yytext,"oneImage %d", &oneImage);
"compcorrn"{S}{I}	sscanf(yytext,"compcorrn %d", &compcorrn);
"usegnuplot"{S}{I}	sscanf(yytext,"usegnuplot %d", &usegnuplot);
"seed"{S}{I}		sscanf(yytext,"seed %d", &seed);
"normInput"{S}{I}	sscanf(yytext,"normInput %d", &normInput);
"noshifting"{S}{I}	sscanf(yytext,"noshifting %d", &noshifting);
"\n"			; /* Remove Newlines */

%%
void static
skipcmnts()
{
for (;;) {
    while ( input() != '*' )
    	;
    if ( input () != '/' ) {
      unput( yytext[yyleng-1]);
    } else return;
}
}


/* routine used by lex when end of input is reached. */
int yywrap()
{
  return(1);
}



/*************************** Version Log ****************************/
/* $Log: storkscan.l,v $
 * Revision 1.2  1995/09/15  12:18:10  stephene
 * updating lex to handle indnwaves()
 *
 *
 */

