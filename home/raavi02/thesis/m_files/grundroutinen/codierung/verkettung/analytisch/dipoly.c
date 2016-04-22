/*** dipoly.c ****************************************************************/
/* Routinen zur Berechnung von Polynomen, Vektoren und Matrizen              */
/* in D und I                                                                */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*****************************************************************************/

#include <limits.h>
#include <malloc.h>
#include <stdio.h>

#include "dipoly.h"

/*** Polynomroutinen *********************************************************/

/* Polynom-Summand an die Liste anhaengen */
DIpolyP
DIpolyAppendTo(DIpolyP next, double Fak, long Dpot, long Ipot)
{
  DIpolyP newDIpoly = (DIpolyP) malloc((size_t) (sizeof(DIpoly)));

  newDIpoly->Fak = Fak;
  newDIpoly->Dpot = Dpot;
  newDIpoly->Ipot = Ipot;

  newDIpoly->next = next;

  return newDIpoly;
}

/* verkettete Liste mit Polynomen freigeben */
void
DIpolyFree(DIpolyP list)
{
  DIpolyP next;

  while (list != NULL)
  {
    next = list->next;
    free(list);
    list = next;
  }
}

/* Debug-Ausgabe fuer ein Polynom */
void
DIpolyPrint(DIpolyP poly)
{
  if (poly != NULL)
    while (poly != NULL)
    {
      printf("%fD^%ldI^%ld ", poly->Fak, poly->Dpot, poly->Ipot);
      poly = poly->next;
    }
  else
    printf("0");
}

/* Terme zusammenfassen */
void
DIpolySimplify(DIpolyP poly)
{
  DIpolyP lastsearchp;
  DIpolyP searchp;

  while (poly != NULL)
  {
    lastsearchp = poly;
    searchp = poly->next;

    while (searchp != NULL)
    {
      if ((searchp->Dpot == poly->Dpot) && (searchp->Ipot == poly->Ipot))
      {
        poly->Fak += searchp->Fak;
        lastsearchp->next = searchp->next;
        free(searchp);
        searchp = lastsearchp->next;
      }
      else
      {
        lastsearchp = searchp;
        searchp = searchp->next;
      }
    }
    poly = poly->next;
  }
}

/* Alle Ipot=0 setzen */
void
DIpolyIOne(DIpolyP poly)
{
  while (poly != NULL)
  {
    poly->Ipot = 0;
    poly = poly->next;
  }
}

/* Ableitung bezueglich Ipot und Ipot=0 setzen */
void
DIpolyDeriv(DIpolyP poly)
{
  while (poly != NULL)
  {
    if (poly->Ipot > 0)
      poly->Fak = poly->Fak * poly->Ipot;
    poly->Ipot = 0;
    poly = poly->next;
  }
}

/* Polynom kopieren */
DIpolyP
DIpolyCopy(DIpolyP src)
{
  DIpolyP dest = NULL;

  while (src != NULL)
  {
    dest = DIpolyAppendTo(dest, src->Fak, src->Dpot, src->Ipot);
    src = src->next;
  }

  return dest;
}

/* zwei Polynome addieren */
DIpolyP
DIpolyAdd(DIpolyP sum1, DIpolyP sum2)
{
  DIpolyP sum = NULL;

  while (sum1 != NULL)
  {
    sum = DIpolyAppendTo(sum, sum1->Fak, sum1->Dpot, sum1->Ipot);
    sum1 = sum1->next;
  }

  while (sum2 != NULL)
  {
    sum = DIpolyAppendTo(sum, sum2->Fak, sum2->Dpot, sum2->Ipot);
    sum2 = sum2->next;
  }

  return sum;
}

/* ein Polynom zu einem anderen hinzuaddieren */
DIpolyP
DIpolyAddTo(DIpolyP sum1, DIpolyP sum2)
{
  while (sum2 != NULL)
  {
    sum1 = DIpolyAppendTo(sum1, sum2->Fak, sum2->Dpot, sum2->Ipot);
    sum2 = sum2->next;
  }

  return sum1;
}

/* zwei Polynome miteinander multiplizieren */
DIpolyP
DIpolyMult(DIpolyP fak1, DIpolyP fak2)
{
  DIpolyP res = NULL;
  DIpolyP fak2loop;

  while (fak1 != NULL)
  {
    fak2loop = fak2;
    while (fak2loop != NULL)
    {
      res = DIpolyAppendTo(res, fak1->Fak * fak2loop->Fak,
                           fak1->Dpot + fak2loop->Dpot,
                           fak1->Ipot + fak2loop->Ipot);
      fak2loop = fak2loop->next;
    }
    fak1 = fak1->next;
  }

  return res;
}

/* zwei Polynome miteinander multiplizieren */
DIpolyP
DIpolyMultdmax(DIpolyP fak1, DIpolyP fak2, long dmax)
{
  DIpolyP res = NULL;
  DIpolyP fak2loop;
  long    newdpot;

  while (fak1 != NULL)
  {
    fak2loop = fak2;
    while (fak2loop != NULL)
    {
      newdpot = fak1->Dpot + fak2loop->Dpot;
      if (newdpot <= dmax)
        res = DIpolyAppendTo(res, fak1->Fak * fak2loop->Fak,
                             newdpot,
                             fak1->Ipot + fak2loop->Ipot);
      fak2loop = fak2loop->next;
    }
    fak1 = fak1->next;
  }

  return res;
}

/* zwei Polynome miteinander multiplizieren */
DIpolyP
DIpolyMultdmaximax(DIpolyP fak1, DIpolyP fak2, long dmax, long imax)
{
  DIpolyP res = NULL;
  DIpolyP fak2loop;
  long    newdpot;
  long    newipot;

  while (fak1 != NULL)
  {
    fak2loop = fak2;
    while (fak2loop != NULL)
    {
      newdpot = fak1->Dpot + fak2loop->Dpot;
      newipot = fak1->Ipot + fak2loop->Ipot;
      if ((newdpot <= dmax) && (newipot <= imax))
        res = DIpolyAppendTo(res, fak1->Fak * fak2loop->Fak,
                             newdpot,
                             newipot);
      fak2loop = fak2loop->next;
    }
    fak1 = fak1->next;
  }

  return res;
}

/* Min Dpot suchen */
long
DIpolyMinDpot(DIpolyP poly)
{
  long    min = LONG_MAX;

  while (poly != NULL)
  {
    if (poly->Dpot < min)
      min = poly->Dpot;

    poly = poly->next;
  }

  return min;
}

/* Max Dpot suchen */
long
DIpolyMaxDpot(DIpolyP poly)
{
  long    max = 0;

  while (poly != NULL)
  {
    if (poly->Dpot > max)
      max = poly->Dpot;

    poly = poly->next;
  }

  return max;
}

/* Max Ipot suchen */
long
DIpolyMaxIpot(DIpolyP poly)
{
  long    max = 0;

  while (poly != NULL)
  {
    if (poly->Ipot > max)
      max = poly->Ipot;

    poly = poly->next;
  }

  return max;
}

/*** Vektoren mit DI-Polynomen ***********************************************/

/* Vektor mit DI-Polynomen erzeugen und initialisieren */
DIpolyVec
DIpolyVecMalloc(long n)
{
  long    count;
  DIpolyVec newvec;

  newvec = (DIpolyVec) malloc((size_t) ((n) * sizeof(DIpolyP)));

  for (count = 0; count < n; count++)
    newvec[count] = NULL;

  return newvec;
}

/* Vektor mit DI-Polynomen freigeben */
void
DIpolyVecFree(DIpolyVec vec, long n)
{
  long    count;

  for (count = 0; count < n; count++)
    DIpolyFree(vec[count]);

  free(vec);
}

/* Vektor mit Einsen erzeugen */
DIpolyVec
DIpolyVecOnes(long vecdim)
{
  DIpolyVec dest;
  long    n;

  dest = DIpolyVecMalloc(vecdim);
  for (n = 0; n < vecdim; n++)
    dest[n] = DIpolyAppendTo(dest[n], 1., 0, 0);

  return dest;
}

/* Vektor mit Polynomen kopieren */
DIpolyVec
DIpolyVecCopy(DIpolyVec src, long vecdim)
{
  DIpolyVec dest;
  long    n;

  dest = DIpolyVecMalloc(vecdim);

  for (n = 0; n < vecdim; n++)
    dest[n] = DIpolyCopy(src[n]);

  return dest;
}

/* Debug-Ausgabe eines Vektors mit DI-Polynomen */
void
DIpolyVecPrint(DIpolyVec vec, long n)
{
  long    i;

  if (vec != NULL)
  {
    for (i = 0; i < n; i++)
    {
      printf("vec[%ld]=", i);
      DIpolyPrint(vec[i]);
      printf("\n");
    }
  }
}

/* Die Polynome eines Vektors vereinfachen */
void
DIpolyVecSimplify(DIpolyVec vec, long n)
{
  long    i;

  if (vec != NULL)
  {
    for (i = 0; i < n; i++)
    {
      DIpolySimplify(vec[i]);
    }
  }
}

/* Skalarmultiplikation - Zeilenvektor mal Spaltenvektor */
DIpolyP
DIpolyRVecMultCVec(DIpolyVec rvec, DIpolyVec cvec, long vecdim)
{
  DIpolyP res = NULL;
  DIpolyP mult = NULL;

  long    count;

  for (count = 0; count < vecdim; count++)
  {
    mult = DIpolyMult(rvec[count], cvec[count]);
    res = DIpolyAddTo(res, mult);
    DIpolyFree(mult);
  }

  return res;
}

/* Skalarmultiplikation - Zeilenvektor mal Spaltenvektor */
DIpolyP
DIpolyRVecMultCVecdmax(DIpolyVec rvec, DIpolyVec cvec, long vecdim, long dmax)
{
  DIpolyP res = NULL;
  DIpolyP mult = NULL;

  long    count;

  for (count = 0; count < vecdim; count++)
  {
    mult = DIpolyMultdmax(rvec[count], cvec[count], dmax);
    res = DIpolyAddTo(res, mult);
    DIpolyFree(mult);
  }

  return res;
}

/* Minimale Dpot im Vektor suchen */
long
DIpolyVecMinDpot(DIpolyVec vec, long dim)
{
  long    min = LONG_MAX;
  long    zeile;

  for (zeile = 0; zeile < dim; zeile++)
  {
    long    newmin = DIpolyMinDpot(vec[zeile]);

    if (newmin < min)
      min = newmin;
  }

  return min;
}

/*** Matrix mit DI-Polynomen *************************************************/

/* Matrix mit DI-Polynomen erzeugen und initialisieren */
DIpolyMtx
DIpolyMtxMalloc(long n)
{
  long    count;

  DIpolyMtx newmtx;

  newmtx = (DIpolyMtx) malloc((size_t) ((n) * sizeof(DIpolyVec)));

  for (count = 0; count < n; count++)
    newmtx[count] = DIpolyVecMalloc(n);

  return newmtx;
}

/* Matrix mit DI-Polynomen freigeben */
void
DIpolyMtxFree(DIpolyMtx mtx, long n)
{
  long    count;

  for (count = 0; count < n; count++)
    DIpolyVecFree(mtx[count], n);

  free(mtx);
}

/* Debug-Ausgabe einer Matrix mit DI-Polynomen */
void
DIpolyMtxPrint(DIpolyMtx mtx, long n)
{
  long    i;
  long    ii;

  if (mtx != NULL)
  {
    for (ii = 0; ii < n; ii++)
      for (i = 0; i < n; i++)
      {
        printf("mtx[%ld][%ld]=", ii, i);
        DIpolyPrint(mtx[ii][i]);
        printf("\n");
      }
  }
}

/* Multiplikation einer Matrix mit einem Spaltenvektor */
DIpolyVec
DIpolyMtxMultCVec(DIpolyMtx Mtx, DIpolyVec vec, long dim)
{
  DIpolyVec res;
  DIpolyP mult;

  long    zeile;
  long    spalte;

  res = DIpolyVecMalloc(dim);
  for (zeile = 0; zeile < dim; zeile++)
  {
    for (spalte = 0; spalte < dim; spalte++)
    {
      mult = DIpolyMult(Mtx[zeile][spalte], vec[spalte]);
      res[zeile] = DIpolyAddTo(res[zeile], mult);
      DIpolyFree(mult);
    }
  }

  return res;
}

/* Multiplikation einer Matrix mit einem Spaltenvektor */
DIpolyVec
DIpolyMtxMultCVecdmax(DIpolyMtx Mtx, DIpolyVec vec, long dim, long dmax)
{
  DIpolyVec res;
  DIpolyP mult;

  long    zeile;
  long    spalte;

  res = DIpolyVecMalloc(dim);
  for (zeile = 0; zeile < dim; zeile++)
  {
    for (spalte = 0; spalte < dim; spalte++)
    {
      mult = DIpolyMultdmax(Mtx[zeile][spalte], vec[spalte], dmax);
      res[zeile] = DIpolyAddTo(res[zeile], mult);
      DIpolyFree(mult);
    }
  }

  return res;
}

/* Multiplikation einer Matrix mit einem Spaltenvektor */
DIpolyVec
DIpolyMtxMultCVecdmaximax(DIpolyMtx Mtx, DIpolyVec vec, long dim,
                          long dmax, long imax)
{
  DIpolyVec res;
  DIpolyP mult;

  long    zeile;
  long    spalte;

  res = DIpolyVecMalloc(dim);
  for (zeile = 0; zeile < dim; zeile++)
  {
    for (spalte = 0; spalte < dim; spalte++)
    {
      mult = DIpolyMultdmaximax(Mtx[zeile][spalte], vec[spalte], dmax, imax);
      res[zeile] = DIpolyAddTo(res[zeile], mult);
      DIpolyFree(mult);
    }
  }

  return res;
}

/*** EOF *********************************************************************/
