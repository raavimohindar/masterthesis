/*** iowef_conv_mex.c ********************************************************/
/* Distanzfunktion eines Faltungscodes berechnen                             */
/* MEX-Intereface                                                            */
/*                                                                           */
/* AUTOR: Juergen Rinas,  31.05.1999                                         */
/*                        03.11.2000                                         */
/*****************************************************************************/

#include <limits.h>
#include <math.h>

#include "dipoly.h"
#include "mex.h"
#include "iowef_conv.h"

void
mexFunction(int nlhs, mxArray * plhs[],
            int nrhs, const mxArray * prhs[])
{
  long    Cl;
  long    n;                    /* Anzahl der Generatorpolynome */
  long    k;                    /* Anzahl der Eingaenge */
  long   *G;                    /* Vektor mit den Generatorpolynomen */
  double *Gin;
  long   *R;                    /* Vektor mit den Generatorpolynomen */
  double *Rin;
  long    dmax;                 /* Abbruch der spektralen Distanz */
  long    jmax;                 /* maximale Schritte */

  long    punctanz;
  long   *punct;
  double *punctin;


  long    i;                    /* Zaehlvariable fuer Kopierfunktionen */

  /* Anzahl der Argumente OK? */
  if (nrhs != 8)
    mexErrMsgTxt("iowef_conv_mex requires eight arguments.");
  if (nlhs != 1)
    mexErrMsgTxt("distfunc_conv_mex requires one output argument");

  if ((mxGetM(prhs[0]) != 1)
      || (mxGetN(prhs[0]) != 1)
      || (mxIsComplex(prhs[0]))
      || (!mxIsDouble(prhs[0]))
      || (mxIsSparse(prhs[0])))
    mexErrMsgTxt("distfunc_conv_mex argument 1 must be a scalar");
  Cl = (long) *mxGetPr(prhs[0]);

  if ((mxGetM(prhs[1]) != 1)
      || (mxGetN(prhs[1]) != 1)
      || (mxIsComplex(prhs[1]))
      || (!mxIsDouble(prhs[1]))
      || (mxIsSparse(prhs[1])))
    mexErrMsgTxt("distfunc_conv_mex argument 2 must be a scalar");
  n = (long) *mxGetPr(prhs[1]);

  if ((mxGetM(prhs[2]) != 1)
      || (mxGetN(prhs[2]) != 1)
      || (mxIsComplex(prhs[2]))
      || (!mxIsDouble(prhs[2]))
      || (mxIsSparse(prhs[2])))
    mexErrMsgTxt("distfunc_conv_mex argument 3 must be a scalar");
  k = (long) *mxGetPr(prhs[2]);

  if ((mxGetM(prhs[3]) != n)
      || (mxGetN(prhs[3]) != 1)
      || (mxIsComplex(prhs[3]))
      || (!mxIsDouble(prhs[3]))
      || (mxIsSparse(prhs[3])))
    mexErrMsgTxt("distfunc_conv_mex argument 4 must be a col-Vector");
  G = (long *) malloc(n * sizeof(long));
  Gin = mxGetPr(prhs[3]);
  for (i = 0; i < n; i++)
    G[i] = (long) (Gin[i]);

  if ((mxGetM(prhs[4]) != n)
      || (mxGetN(prhs[4]) != 1)
      || (mxIsComplex(prhs[4]))
      || (!mxIsDouble(prhs[4]))
      || (mxIsSparse(prhs[4])))
    mexErrMsgTxt("distfunc_conv_mex argument 5 must be a col-Vector");
  R = (long *) malloc(n * sizeof(long));
  Rin = mxGetPr(prhs[4]);
  for (i = 0; i < n; i++)
  {
    R[i] = (long) (Rin[i]);
  }

  if ((mxGetM(prhs[5]) != 1)
      || (mxGetN(prhs[5]) != 1)
      || (mxIsComplex(prhs[5]))
      || (!mxIsDouble(prhs[5]))
      || (mxIsSparse(prhs[5])))
    mexErrMsgTxt("distfunc_conv_mex argument 6 must be a scalar");
  if (mxIsInf(*mxGetPr(prhs[5])))
    dmax = LONG_MAX;
  else
    dmax = (long) *mxGetPr(prhs[5]);

  if ((mxGetM(prhs[6]) != 1)
      || (mxGetN(prhs[6]) != 1)
      || (mxIsComplex(prhs[6]))
      || (!mxIsDouble(prhs[6]))
      || (mxIsSparse(prhs[6])))
    mexErrMsgTxt("distfunc_conv_mex argument 7 must be a scalar");
  if (mxIsInf(*mxGetPr(prhs[6])))
    jmax = LONG_MAX;
  else
    jmax = (long) *mxGetPr(prhs[6]);

  punctanz = mxGetM(prhs[7]);
  if ((mxGetN(prhs[7]) != 1)
      || (mxIsComplex(prhs[7]))
      || (!mxIsDouble(prhs[7]))
      || (mxIsSparse(prhs[7])))
    mexErrMsgTxt("iowef_conv_mex argument 6 must be a col-Vector");
  punct = (long *) malloc(punctanz * sizeof(long));
  punctin = mxGetPr(prhs[7]);
  for (i = 0; i < punctanz; i++)
    punct[i] = (long) (punctin[i]);


  if (nlhs == 1)
  {
    DIpolyP res;
    DIpolyP poly;
    double *retmtx;
    long imax;


    res=iowef_conv(Cl, n, k, G, R, dmax, jmax,punctanz,punct);
   
    imax = DIpolyMaxIpot(res);


    /* unendlich-Grosse Matrix als Rueckgabe nicht realisierbar */
    /* Begrenzung auf berechnete D-Potenzen */
    if (dmax == LONG_MAX)
      dmax = DIpolyMaxDpot(res);


    plhs[0] = mxCreateDoubleMatrix(imax + 1, dmax + 1, mxREAL);
    retmtx = mxGetPr(plhs[0]);

    poly = res;
    while (poly != NULL)
    {
      if ((poly->Dpot < dmax + 1) && (poly->Ipot < imax + 1))
        retmtx[poly->Dpot * (imax + 1) + poly->Ipot] = poly->Fak;
      poly = poly->next;
    }
    DIpolyFree(res);
  }

  free(punct);
  free(G);
  free(R);
}

/*** EOF *********************************************************************/
