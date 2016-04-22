/*
% --------------------------------------------------------------------------------
%
%    function n_over_k = bin_coef(n,k);
%
%    computation of binomial coefficients
%
% --------------------------------------------------------------------------------
%
%
%
% input parameter       n: vector
%			k: vector
%
%                       length(k) must be length(n) !!!!!
%
%
% output 		n_over_k: vector
%                                 n_over_k(i)  =  ( n(i) )
%                                                 ( k(i) )
%
% used mex-files:       bin_coef.mexlx, bin_coef.mexsol
%
% ---------------------------------------------------------------------------------
%
%  Author:		Heiko Schmidt (University of Bremen)
%  Project:             ---
%  Date:		25-Aug-1999
%  Last modified:  	30-Aug-1999
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

  double    *input_n;             /* real input vector */
  double    *input_k;             /* imag input vector */

  double    *output;	             /* output vector */

 
 
  
/*-------------------------------------------------------------------------------------
 * locale declarations:
 *------------------------------------------------------------------------------------- */


  int		lauf,input_len,i;
  
  double   	n,k,l,fn,fk,fnmk,tmp,rest,k_tmp,l_tmp;
  


/*--------------------------- Input Declaration --------------------------------*/

  input_n  	= mxGetPr(prhs[0]);
  input_k 	= mxGetPr(prhs[1]);

  input_len     = mxGetN(prhs[0]) * mxGetM(prhs[0]);

  
/*-------------------------- Output Declaration ----------------------------*/

  plhs[0]     = mxCreateDoubleMatrix(input_len,1,mxCOMPLEX);
  output      = mxGetPr(plhs[0]);

/*------------------------------ COMPUTATION -----------------------------*/


  for(lauf=0;lauf<input_len;lauf++){   /* loop of length(n)  */
      
      n    = input_n[lauf];
      k    = input_k[lauf];
      l    = n-k;

      tmp  = 1.0;
      if (l>=k) {
	  k_tmp = k;
	  for(i=l+1;i<=n;i++){
	      tmp = tmp * (double)i;
	      if (k_tmp > 1){
		  tmp = tmp / k_tmp;
		  k_tmp = k_tmp - 1;
	      }
	  }
      }
      else{
	  l_tmp = l;
      	  for(i=k+1;i<=n;i++){
	      tmp = tmp * (double)i;
              if (l_tmp > 1){
		  tmp = tmp / l_tmp;
		  l_tmp = l_tmp - 1;
	      }
	  }
      }

      output[lauf] = tmp;
  } /* END OF LOOP */

 }

/* END OF FUNTION*/


