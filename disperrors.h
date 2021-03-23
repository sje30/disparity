/****************************************************************************
***
***
*** disperrors.h
*** 
*** Stephen Eglen
*** COGS, University of Sussex.
***
*** Created 21 Nov 95
***
*** $Revision$
*** $Date$
****************************************************************************/


#ifndef _DISPERRORS_H
#define _DISPERRORS_H

void	storeTopLayerErrors(Array da);
void	propagateErrors();
void	calcErrors( int layer, int vec);
double	dtanh( double x);
void	createPartials();
#endif


/*************************** Version Log ****************************/
/* $Log$
 *
 */
