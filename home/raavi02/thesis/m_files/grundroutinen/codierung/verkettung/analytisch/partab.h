/*** partab.h*****************************************************************/
/* Berechnung der Paritaet bzw. MOD2 Addition ueber Tabellen                 */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*****************************************************************************/

#ifndef _partab_h
#define _partab_h

/* Tabelle mit 8-Bit Paritaet */
extern unsigned char partab[];

/* 16-Bit Paritaet berechnen */
int     parity16(int x);

/* 32-Bit Paritaet berechnen */
long    parity32(long x);

#endif

/*** EOF **********************************************************************/
