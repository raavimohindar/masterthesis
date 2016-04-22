/*** iowef_block_conv.c ******************************************************/
/* IOWEF eines Faltungscodes berechnen                                       */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*                        03.11.2000                                         */
/*****************************************************************************/

#include <malloc.h>
#include <math.h>
#include <stdio.h>

#include "iowef_block_conv.h"
#include "dipoly.h"
#include "partab.h"

/*
#define DEBUG_PRINTFS
#define PRINT_PROGRESS
*/

DIpolyP
iowef_block_conv(long Cl, long n, long k, long *G, long *R,
           long htrunc, long wtrunc, long wmax, long punctanz, long *punct)
{
  /* Vektoren zur Berechnung */
  DIpolyVec a;
  DIpolyVec b;
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

  long    vecdim;

  long    p;
  long    punctstate;

  long    w;                    /* Zaehler fuer die Matrizenmultiplikationen */
  long    ROK;


  /* Vektoren zur Berechnung erzeugen */
  vecdim = (1 << (Cl - k));

  a = DIpolyVecMalloc(vecdim);
  b = DIpolyVecMalloc(vecdim);
  S = (DIpolyMtx *) malloc(punctanz * sizeof(DIpolyMtx *));
  for (p = 0; p < punctanz; p++)
    S[p] = DIpolyMtxMalloc(vecdim);

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
        z = startz | ((Idata ^ parity16(startz & ROK)) << (Cl - k));

        endz = (z >> k) & ((1 << (Cl - k)) - 1);  /* Maskierung Endzustand */

        Dpot = 0;
        for (ncount = 0; ncount < n; ncount++) /* alle Ausgaenge durchzaehlen */
        {                        
          if ((ROK != 0) && (R[ncount] == 0) && (G[ncount] == 1))
             /* Sonderfall: systematische Bits eines rekursiven Codes */
            Dpot += (punct[p] >> ncount) & 0x01 & parity16(Idata); 
          else
             /* Normalfall: Maskierung der Zustandsbits */
            Dpot += (punct[p] >> ncount) & 0x01 & parity16(G[ncount] & z);
        }

        S[p][endz][startz] = DIpolyAppendTo(S[p][endz][startz], 1, Dpot, Ipot);
      }
    }
  }

#ifdef DEBUG_PRINTFS
  printf("Matrix OK\n");
  fflush(stdout);
#endif

  /* Vektor für Start- und Endzustand */
  b[0] = DIpolyAppendTo(b[0], 1, 0, 0);
  a[0] = DIpolyAppendTo(a[0], 1, 0, 0);


#ifdef DEBUG_PRINTFS
  /* Debug-Ausgabe der Matrizen und Vektoren */
  printf("a:\n");
  DIpolyVecPrint(a, vecdim);
  printf("b:\n");
  DIpolyVecPrint(b, vecdim);
  for (p = 0; p < punctanz; p++)
  {
    printf("punct pattern %ld\n", p);
    printf("S:\n");
    DIpolyMtxPrint(S[p], vecdim);
  }
#endif

  /* Zustand der Punktierung */
  punctstate = 0;              

  /* Codespektrum berechnen */
  lastdx = DIpolyVecCopy(b, vecdim);

  for (w = 1; w <= (wmax / k); w++)
  {
    dx = DIpolyMtxMultCVecdmaximax(S[punctstate], lastdx, vecdim,
                                   htrunc, wtrunc);

    punctstate = (punctstate + 1) % punctanz;

#ifdef PRINT_PROGRESS
    printf("Trellissegment (Nutzdaten) Nr.:%ld\n", w);
#endif

    DIpolyVecSimplify(dx, vecdim);

    DIpolyVecFree(lastdx, vecdim);
    lastdx = DIpolyVecCopy(dx, vecdim);
    DIpolyVecFree(dx, vecdim);
  }

  /* Terminierung */
  {
    /* Anzahl der Trellissegmente zur Terminierung bestimmen*/
    long    termtrellisseg;
    termtrellisseg = (long) ceil(((double) (Cl - k)) / ((double) k));

    /* Die Terminierung leistet keinen Beitrag zur
       Erhöhung des Informationsgewichtes....
       Daher werden Die I-Polinome hier auf Eins gesetzt.
    */

/*
   Wenn folgender Abschnitt auskommentiert ist, werden
   die Info-Gewichte der Terminierung in der IOWEF mitgezählt.
     Dieses kann für verkettete Codes sinvoll sein!
   
   Ist der folgende Abschnitt aktiv, so
   werden die I-Terme (Info-Gewichte) zu 1 gesetzt 
   bzw. die I-Potenzen auf 0 festgelegt, so daß die Tail-Bits 
   keinen Beitrag zu den Info-Gewicht liefern.
     Dieses ist für reine Faltungscodes (ohne Verkettung)
     sinvoll.
*/
/*
    for (p = 0; p < punctanz; p++)
      for (startz = 0; startz < vecdim; startz++)
        for (endz = 0; endz < vecdim; endz++)
          if (S[p][endz][startz] != NULL)
            DIpolyIOne(S[p][endz][startz]);
*/

    /* Terminierung durchführen */
    for (w = 1; w <= termtrellisseg; w++)
    {
      dx = DIpolyMtxMultCVecdmax(S[punctstate], lastdx, vecdim, htrunc);
      punctstate = (punctstate + 1) % punctanz;

#ifdef PRINT_PROGRESS
      printf("Trellissegment (Terminierung) Nr.:%ld\n", w);
#endif
      DIpolyVecSimplify(dx, vecdim);

      DIpolyVecFree(lastdx, vecdim);
      lastdx = DIpolyVecCopy(dx, vecdim);
      DIpolyVecFree(dx, vecdim);
    }
  }
  res = DIpolyRVecMultCVecdmax(a, lastdx, vecdim, htrunc);

  DIpolyVecFree(lastdx, vecdim);

  DIpolyVecFree(a, vecdim);
  DIpolyVecFree(b, vecdim);
  for (p = 0; p < punctanz; p++)
    DIpolyMtxFree(S[p], vecdim);

  DIpolySimplify(res);

  return res;
}

/*** EOF *********************************************************************/
