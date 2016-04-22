/*** iowef_conv.c ************************************************************/
/* Distanzfunktion eines Faltungscodes berechnen                             */
/* Quelle: Friedrichs:                                                       */
/*         Kanalcodierung                                                    */
/*         Springer 1994                                                     */
/*         S. 264ff.                                                         */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*                        05.11.2000                                         */
/*****************************************************************************/

#include <malloc.h>
#include <math.h>
#include <stdio.h>
 
#include "iowef_conv.h" 
#include "dipoly.h"
#include "partab.h"


/* #define DEBUG_PRINTFS */

DIpolyP iowef_conv(long Cl, long n, long k, long *G, long *R, long dmax, long jmax, long punctanz, long *punct)
{
  /* Vektoren zur Berechnung */
  DIpolyVec *a;
  DIpolyVec *b;
  DIpolyMtx *S;

  DIpolyVec dx = NULL;
  DIpolyVec lastdx = NULL;
  DIpolyP res = NULL;
  DIpolyP adx = NULL;

  long    Dpot;
  long    Ipot;

  long    Idata;

  long    z;                    /* alle Zustaende durchzaehlen */
  long    startz;               /* Startzustand */
  long    endz;                 /* Endzustand */

  long    ncount;
  long    kcount;


  long    p;                    /* Zähler für Punktierungszwecke */
  long    punctstate;           /* Zähler für den Punktierungszustand */

  long    minDpot;              /* minimale D-Potenz bei der Berechnung */

  long    vecdim;

  long    j;                    /* Zaehler fuer Matrizenmultiplikationen */
  long    ROK;

  /* Vektoren zur Berechnung erzeugen */
  vecdim = (1 << (Cl - k)) - 1;

  a = (DIpolyVec *) malloc(punctanz * sizeof(DIpolyVec *));
  b = (DIpolyVec *) malloc(punctanz * sizeof(DIpolyVec *));
  S = (DIpolyMtx *) malloc(punctanz * sizeof(DIpolyMtx *));
  for (p = 0; p < punctanz; p++)
  {
    a[p] = DIpolyVecMalloc(vecdim);
    b[p] = DIpolyVecMalloc(vecdim);
    S[p] = DIpolyMtxMalloc(vecdim);
  }  


  ROK = 0;
  for (ncount = 0; ncount < n; ncount++)
  {
    if (R[ncount] != 0)
    {
      if ((ROK != 0) && (ROK != R[ncount]))
        fprintf(stderr, "unterschiedliche Rueckkopplungen! Verwende %o\n",
                ROK);
      else
        ROK = R[ncount];
    }
  }

  if ((ROK != 0) && (k != 1))
  {
    fprintf(stderr, "Rekursive Codes mit k!=1 nicht implementiert\n");
  }


  /* Vektoren a,b und Matrix C für alle Punktierungszustände berechnen */
  for (p = 0; p < punctanz; p++)
  {
  for (startz = 0; startz < (1 << (Cl - k)); startz++)
  {
    for (Idata = 0; Idata < (1 << k); Idata++)
    {
      /* Gewicht der Eingangsdaten bestimmen */
      /* alle k Eingangsbits durchzaehlen */
      Ipot = 0;
      for (kcount = 0; kcount < k; kcount++)
        if ((Idata >> kcount) & 0x01)
          Ipot++;

      /* Inhalt des Zustandsregisters */
      /* von links mit neuen Daten aufüllen */
      z = startz | ((Idata ^ parity16(startz & ROK)) << (Cl - k));

      /* Maskierung Endzustand */
      endz = (z >> k) & ((1 << (Cl - k)) - 1);  

      Dpot = 0;
      for (ncount = 0; ncount < n; ncount++)  /* alle Ausgaenge durchzaehlen */
      {
        if ((ROK != 0) && (R[ncount] == 0) && (G[ncount] == 1))
          /* Sonderfall: Systematische Bits eines rekursiven Codes */
          Dpot += (punct[p] >> ncount) & 0x01 & parity16(G[ncount] & Idata);
        else
          /* Normalfall: Maskierung der Zustandsbits */
          Dpot += (punct[p] >> ncount) & 0x01 & parity16(G[ncount] & z);
      }


      if ((startz != 0) && (endz != 0))
        S[p][endz - 1][startz - 1] = DIpolyAppendTo(S[p][endz - 1][startz - 1], 1,Dpot, Ipot);

      if ((startz == 0) && (endz!=0))
        b[p][endz - 1] = DIpolyAppendTo(b[p][endz - 1], 1, Dpot, Ipot);


      if ((endz == 0) && (startz!=0))
          a[p][startz - 1] = DIpolyAppendTo(a[p][startz - 1], 1, Dpot, Ipot);


#ifdef DEBUG_PRINTFS
      printf("\n");
      printf("Dpot=%ld\n", Dpot);
      printf("Ipot=%ld\n", Ipot);
      printf("z=%ld\n", z);
      printf("Idata %ld\n", Idata);
      printf("startz %ld\n", startz);
      printf("endz %ld\n", endz);
#endif
    }
  }
  }


#ifdef DEBUG_PRINTFS
  for (p=0;p<punctanz;p++)
  {
    printf("a:\n");
    DIpolyVecPrint(a[p], vecdim);
    printf("\n");

    printf("b:\n");
    DIpolyVecPrint(b[p], vecdim);
    printf("\n");

    printf("S:\n");
    DIpolyMtxPrint(S[p], vecdim);
    printf("\n");
  }
#endif


  
  /* Punktierung beginnt mit dem  ersten Punktierungsmuster (Index [0]) */
  punctstate=0; 

  /* Codespektrum berechnen */
  adx = DIpolyRVecMultCVec(a[(punctstate+1)%punctanz], b[punctstate], vecdim);
  res = DIpolyAddTo(res, adx);
  DIpolyFree(adx);

  lastdx = DIpolyVecCopy(b[punctstate], vecdim);


  punctstate=(punctstate+1)%punctanz;
  
  j = 1;
  while (j < jmax)
  {
    j++;
    dx = DIpolyMtxMultCVecdmax(S[punctstate], lastdx, vecdim,dmax);


    minDpot = DIpolyVecMinDpot(dx, vecdim); 
    /* falls der Vektor dx leer ist, liefert DIpolyVecMinDpot LONG_MAX */
#ifdef DEBUG_PRINTFS
    printf("minDpot=%ld\n", minDpot);
#endif
    /* Abbruch, falls für Distanzen<=dmin 
       keine Änderung mehr zu erwarten ist!
    */
    if (minDpot > dmax)
      break;

    DIpolyVecSimplify(dx, vecdim);
    
    punctstate=(punctstate+1)%punctanz;
    adx = DIpolyRVecMultCVec(a[punctstate], dx, vecdim);

    res = DIpolyAddTo(res, adx);
    DIpolyFree(adx);

    DIpolySimplify(res);

    DIpolyVecFree(lastdx, vecdim);
    lastdx = DIpolyVecCopy(dx, vecdim);
    DIpolyVecFree(dx, vecdim);
  }


  DIpolyVecFree(lastdx, vecdim);

  for (p = 0; p < punctanz; p++)
  {
    DIpolyVecFree(a[p], vecdim);
    DIpolyVecFree(b[p], vecdim);
    DIpolyMtxFree(S[p], vecdim);
  }  


  DIpolySimplify(res);

  return res;
}

/*** EOF *********************************************************************/
