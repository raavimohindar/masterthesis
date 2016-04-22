/*** iowef_block_conv.h ******************************************************/
/* IOWEF eines Faltungscodes berechnen                                       */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*                        03.11.2000                                         */
/*****************************************************************************/

#ifndef _iowef_block_conv_h
#define _iowef_block_conv_h

#include "dipoly.h"

DIpolyP iowef_block_conv(long Cl, long n, long k, long *G,
    long *R, long htrunc, long wtrunc, long wmax, long punctanz, long *punct);

#endif

/*** EOF *********************************************************************/
