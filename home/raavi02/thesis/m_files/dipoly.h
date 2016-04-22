/*** dipoly.h ****************************************************************/
/* Routinen zur Berechnung von Polynomen, Vektoren und Matrizen              */
/* in D und I                                                                */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*****************************************************************************/

#ifndef _dipoly_h
#define _dipoly_h

/* Verkettete Liste zur Darstellung der Polynome in I und D */
struct DIpoly
{
  double  Fak;                  /* Vorfaktor */
  long    Dpot;                 /* Potenz des D-Anteils */
  long    Ipot;                 /* Potenz des I-Anteils */

  struct DIpoly *next;          /* Zeiger auf naechsten Summanden */
};

typedef struct DIpoly DIpoly;
typedef DIpoly *DIpolyP;
typedef DIpolyP *DIpolyVec;
typedef DIpolyVec *DIpolyMtx;

/* Polynomroutinen */

  /* Polynom-Summand an die Liste anhaengen */
DIpolyP DIpolyAppendTo(DIpolyP next, double Fak, long Dpot, long Ipot);

  /* verkettete Liste mit Polynomen freigeben */
void    DIpolyFree(DIpolyP list);

  /* Debug-Ausgabe fuer ein Polynom */
void    DIpolyPrint(DIpolyP poly);

  /* Terme zusammenfassen */
void    DIpolySimplify(DIpolyP poly);

  /* Alle Ipot=0 setzen */
void    DIpolyIOne(DIpolyP poly);

  /* Ableitung bezueglich Ipot und Ipot=0 setzen */
void    DIpolyDeriv(DIpolyP poly);

  /* Polynom kopieren */
DIpolyP DIpolyCopy(DIpolyP src);

  /* zwei Polynome addieren */
DIpolyP DIpolyAdd(DIpolyP sum1, DIpolyP sum2);

  /* ein Polynom zu einem anderen hinzuaddieren */
DIpolyP DIpolyAddTo(DIpolyP sum1, DIpolyP sum2);

  /* zwei Polynome miteinander multiplizieren */
DIpolyP DIpolyMult(DIpolyP fak1, DIpolyP fak2);

  /* zwei Polynome miteinander multiplizieren */
DIpolyP DIpolyMultdmax(DIpolyP fak1, DIpolyP fak2, long dmax);

  /* zwei Polynome miteinander multiplizieren */
DIpolyP DIpolyMultdmaximax(DIpolyP fak1, DIpolyP fak2, long dmax, long imax);

  /* Min Dpot suchen */
long    DIpolyMinDpot(DIpolyP poly);

  /* Max Dpot suchen */
long    DIpolyMaxDpot(DIpolyP poly);

  /* Max Ipot suchen */
long    DIpolyMaxIpot(DIpolyP poly);

/* Vektoren mit DI-Polynomen */

  /* Vektor mit DI-Polynomen erzeugen und initialisieren */
DIpolyVec DIpolyVecMalloc(long n);

  /* Vektor mit DI-Polynomen freigeben */
void    DIpolyVecFree(DIpolyVec vec, long n);

  /* Vektor mit Einsen erzeugen */
DIpolyVec DIpolyVecOnes(long vecdim);

  /* Vektor mit Polynomen kopieren */
DIpolyVec DIpolyVecCopy(DIpolyVec src, long vecdim);

  /* Debug-Ausgabe eines Vektors mit DI-Polynomen */
void    DIpolyVecPrint(DIpolyVec vec, long n);

  /* Die Polynome eines Vektors vereinfachen */
void    DIpolyVecSimplify(DIpolyVec vec, long n);

  /* Skalarmultiplikation - Zeilenvektor mal Spaltenvektor */
DIpolyP DIpolyRVecMultCVec(DIpolyVec rvec, DIpolyVec cvec, long vecdim);

  /* Skalarmultiplikation - Zeilenvektor mal Spaltenvektor */
DIpolyP DIpolyRVecMultCVecdmax(DIpolyVec rvec, DIpolyVec cvec,
                               long vecdim, long dmax);

  /* Minimale Dpot im Vektor suchen */
long    DIpolyVecMinDpot(DIpolyVec vec, long dim);

/* Matrizen mit DI-Polynomen */

  /* Matrix mit DI-Polynomen erzeugen und initialisieren */
DIpolyMtx DIpolyMtxMalloc(long n);

  /* Matrix mit DI-Polynomen freigeben */
void    DIpolyMtxFree(DIpolyMtx mtx, long n);

  /* Debug-Ausgabe einer Matrix mit DI-Polynomen */
void    DIpolyMtxPrint(DIpolyMtx mtx, long n);

  /* Multiplikation einer Matrix mit einem Spaltenvektor */
DIpolyVec DIpolyMtxMultCVec(DIpolyMtx Mtx, DIpolyVec vec, long dim);

  /* Multiplikation einer Matrix mit einem Spaltenvektor */
DIpolyVec DIpolyMtxMultCVecdmax(DIpolyMtx Mtx, DIpolyVec vec, long dim,
                                long dmax);

  /* Multiplikation einer Matrix mit einem Spaltenvektor */
DIpolyVec DIpolyMtxMultCVecdmaximax(DIpolyMtx Mtx, DIpolyVec vec, long dim,
                                    long dmax, long imax);

#endif

/*** EOF *********************************************************************/
