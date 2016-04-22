/*** iowef_conv.h ************************************************************/
/* Distanzfunktion eines Faltungscodes berechnen                             */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*****************************************************************************/

#ifndef _iowef_conv_h
#define _iowef_conv_h

#include "dipoly.h"

DIpolyP iowef_conv(long Cl, long n, long k, long *G, long *R,
                      long dmax, long jmax,long punctanz, long *punct);

#endif

/*** EOF *********************************************************************/
