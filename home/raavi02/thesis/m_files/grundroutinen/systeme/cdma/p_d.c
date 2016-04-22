/*
% --------------------------------------------------------------------------------
%
%    function P_d = p_d(d,mu);
%
%    calculating the pairwise error probability for a fading channel
%
% --------------------------------------------------------------------------------
%
%
%
% input parameter       d : scalar, hamming distance to be considered
%		      mu: vector, mu=sqrt(gamma/(1+gamma)), includes SNR per channel
%
%                       length(k) must be length(n) !!!!!
%
%
% output 		      P_d: vector of pairwise error probabilities
%
%
% used mex-files:       p_d.mexlx, p_d.mexsol
%
% ---------------------------------------------------------------------------------
%
%  Author:		Volker Kühn (University of Bremen)
%  Project:             ---
%  Date:		27-Aug-1999
%  Last modified:  	27-Aug-1999
%  Matlab Version:	5.3
%
%  Copyright (c) Department of Telecommunications / University of Bremen 1999
%  This software is property of University of Bremen. Unauthorized 
%  duplication and disclosure to third parties is forbidden.	
% ---------------------------------------------------------------------------------
%
*/
#include <math.h>
#include "mex.h"


/*-------------------------------------------------------------------------------------
 * start MEX-Function
 *------------------------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

/*-------------------------------------------------------------------------------------
 * in/out declarations:
 *------------------------------------------------------------------------------------- */

  double     input_d;             /* real input scalar */
  double    *input_mu;            /* real input vector */
  double     prec;                /* real input scalar, wanted precision */

  double    *output;	          /* output vector */
  double    *overflow;            /* output vector */
 
 
  
/*-------------------------------------------------------------------------------------
 * locale declarations:
 *------------------------------------------------------------------------------------- */


  int		lauf,output_len,j,reduce;
  
  double   	i,n,k,l,tmp,tmp1,k_tmp,mu1,mu2;
  


/*--------------------------- Input Declaration --------------------------------*/

  input_d  	= mxGetScalar(prhs[0]);
  input_mu	= mxGetPr(prhs[1]);
  prec  	= mxGetScalar(prhs[2]);

  output_len    = mxGetN(prhs[1]) * mxGetM(prhs[1]);

  
/*-------------------------- Output Declaration ----------------------------*/

  plhs[0]     = mxCreateDoubleMatrix(output_len,1,mxCOMPLEX);
  output      = mxGetPr(plhs[0]);
  plhs[1]     = mxCreateDoubleMatrix(output_len*4,1,mxCOMPLEX);
  overflow    = mxGetPr(plhs[1]);


/*------------------------------ COMPUTATION -----------------------------*/


  for (lauf=0; lauf<output_len; lauf++) {            /* loop over all mu's */
    mu1 = (1.0 - input_mu[lauf]) / 2.0;
    mu2 = (1.0 + input_mu[lauf]) / 2.0;
    overflow[lauf]              = 0.0;               /*  */
    overflow[lauf+output_len]   = 0.0;
    overflow[lauf+2*output_len] = 0.0;
    overflow[lauf+3*output_len] = 0.0;

    tmp = mu1;                                       /* Initialization for j=0 */
    j = 1;
    while ((j<(int)input_d) & (tmp>1.0e-300)) {
         output[lauf] = tmp;
         tmp *= mu1;
         j++;
    }
    if ((j==(int)input_d) & (tmp>1.0e-300)) 
       output[lauf] = tmp;                           /* store last value if full loop */
    else 
       output[lauf] = 0.0;

    tmp = 1.0;
    j   = (int)input_d - 1;
    while ((j>0) & (tmp/output[lauf]>prec)) {        /* sum over j=1:d-1 */
        n    = input_d - 1.0 + (double)j;
        k    = (double)j;
        l    = n - k;

        reduce = 0;

        tmp   = 1.0;
        k_tmp = k;
        i     = l + 1.0;
        while ((i<=n) & (tmp>1.0e-300)) {
            tmp *= i * mu1 * mu2;
            i   += 1.0;
            if (k_tmp > 1.0){
		        tmp   /= k_tmp;
		        k_tmp -= 1.0;
            }
            if (tmp>1.0e200) {
                tmp *= 1.0e-100;
                reduce++;
	    }
	}
	/*printf("Summe: j=%d\ttmp=%e\tmu1*mu2=%e\treduce=%d\n",j,tmp,mu1*mu2,reduce);*/

        if ((tmp > 1.0e300) | (tmp<1.0e-300)) {
		   if (tmp < 1.0e-300) overflow[lauf] += 1.0;
		   else     overflow[lauf+output_len] += 1.0;
        }
        else {
            tmp1 = tmp;
            i    = k;
            while ((i<input_d) & (tmp1>1.0e-300) & (tmp1<1.0e300)) { 
               tmp   = tmp1;
               tmp1 *= mu1;
               i    += 1.0;
               if ((tmp1<1.0e-200) & (reduce > 0)) {      /* correction of reduction */
                  tmp1 *= 1.0e100;
                  reduce--;
	       }
            }
	    /*printf("Ende: input_d=%f\tj=k=%d\ttmp=%e\t,output=%e\treduce=%d\n",input_d,j,tmp,output[lauf],reduce);*/

            while ((reduce>0) & (tmp1<1.0e200)) {
                  tmp1 *= 1.0e100;
                  reduce--;
            }
            if ((tmp1>1.0e-300) & (reduce==0)) {
      	        tmp = tmp1;                               /* store last value if full loop */
            }
            else {
		if (tmp1<1.0e-300) overflow[lauf+2*output_len] += 1.0;
		else               overflow[lauf+3*output_len] += 1.0;
		/*printf("Abbruch mu1^d: j=%d, i=%f, tmp=%e\n",j,i,tmp);*/
            }
        }
        output[lauf] += tmp;
        j--;

    } /* END OF LOOP j */

  }  /* END OF LOOP lauf */

}/* END OF FUNTION*/


