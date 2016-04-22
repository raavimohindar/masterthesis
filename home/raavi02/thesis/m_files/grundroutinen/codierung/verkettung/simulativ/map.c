/* Function file: SS-MAP-Algorithm for Convolutional Codes (RSC and NSC)
 * Copyright 2000 Volker Kuehn
 * for academic use only
 *
 * function [L_info,L_code]=map(signal, code);
 *
 * -------------------------------------------------------------------------------------------------
 * PROJECT:   Channel-Coding-Algorithmn
 * Submitted by: University of Bremen
 *  
 * -------------------------------------------------------------------------------------------------
 * DESCRIPTION:
 * -------------------------------------------------------------------------------------------------
 * 1. SS-MAP-Decoding of Convolutional Codes supplying LLR's of information bits and codebits
 * 2. Suitable for RSC and NSC Convolutional Codes 
 * -------------------------------------------------------------------------------------------------
 * INPUT:
 * -------------------------------------------------------------------------------------------------
 * signal:    struct containing the following fields
 *                signal.sig        : input signal containing data
 *                signal.L_a        : a priori information for each information bit
 *                signal.last_state : represents last state of trellis if code.term==1
 *
 * code   struct containing code parameters
 *                code.trellis_out  : describes code words corresponding to specific state transitions
 *                code.trellis_next : describes succesive state dependent on preceeding state and input
 *                code.num_state    : number of states in trellis
 *                code.block_len    : number of code words per frame
 *                code.word_len     : number of bits per code word
 *                code.term         : terminated trellis diagram: term == 1
 *                                           otherwise: term~=1
 * -------------------------------------------------------------------------------------------------
 * OUTPUT:
 * -------------------------------------------------------------------------------------------------
 * L_info         LLR's for decoded information bits
 * L_code          LLR's of coded bits (optional)
 * -------------------------------------------------------------------------------------------------
 * FUNCTIONS:
 * -------------------------------------------------------------------------------------------------
 * - function ALLOCProb (included in this file) allocates memory for 2-dim. fields and returns 
 *   pointer on this field
 * -------------------------------------------------------------------------------------------------
 * REMARKS:
 * -------------------------------------------------------------------------------------------------
 * - For further description
 *   see Berrou, C.; Glavieux, A.; Thitimajshima, P.
 *       Near Shannon Limit Error-Correction Coding and Decoding: Turbo_codes(1)
 *       ICC-93, Geneva, pp.1064-1070, May93 
 * -------------------------------------------------------------------------------------------------
 * IMPORTANT PARAMETERS:
 * -------------------------------------------------------------------------------------------------
 * - none
 * * -------------------------------------------------------------------------------------------------
 * Author: Volker Kuehn 06.12.00
 *         last changes from Volker Kuehn 03.04.01, 29.08.01
 * * -------------------------------------------------------------------------------------------------
 */

#include <math.h>
#include "mex.h"
#include "matrix.h"
 
/*-------------------------------------------------------------------------------------
 * S T A R T    F U N C T I O N 
 *------------------------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mxArray *field_array_ptr;
  double *sig;           /* Eingangsvektor */
  double *Z;             /* a-priori information */
  double *L_code;         /* Zeiger auf Decodierergebnis */
  double *L_info;        /* Zeiger auf Decodierergebnis fuer Infobits*/
  double *trellis_out;   /* array containing outputs of trellis branches */
  double *trellis_next; /* array containing successiv states */
  int block_len;         /* number of codewords / frame */
  int n;                 /* codeword length */
  int num_state;         /* number of states in trellis */
  int known_end;         /* terminated trellis ? */
  int last_state;        /* indicates last state of trellis if known_end==1 */
  int i,j,k,m;

/*-------------------------------------------------------------------------------------
 * PROBABILITY_DECLARATIONS:
 *------------------------------------------------------------------------------------- */
  double    *alpha;
  double    *gamma0;    /* Array of channel transition probabilities AND a priori information for d_k=0 */
  double    *gamma1;    /* Array of channel transition probabilities AND a priori information for d_k=0 */
  double    *lambda0;
  double    *lambda1;
  double    *beta;
  double    *norm;
  double    *metric;
  double     dummy;
  int        num_word;
  int        out_info0;
  int        out_info1;
  double     L0;
  double     L1;


if (nrhs!=2) {
   mexErrMsgTxt("Es werden 2 Eingabeargumente erwartet!");
}


/*-------------------------------------------------------------------------------------
 * INPUT_DECLARATIONS:
 *------------------------------------------------------------------------------------- */
  field_array_ptr = mxGetField(prhs[0],0,"sig");
  sig = mxGetPr(field_array_ptr);

  field_array_ptr = mxGetField(prhs[0],0,"L_a");
  Z = mxGetPr(field_array_ptr);

  field_array_ptr = mxGetField(prhs[0],0,"last_state");
  last_state = (int)mxGetScalar(field_array_ptr);


  field_array_ptr = mxGetField(prhs[1],0,"trellis_out");
  trellis_out = mxGetPr(field_array_ptr);

  field_array_ptr = mxGetField(prhs[1],0,"trellis_next");
  trellis_next = mxGetPr(field_array_ptr);

  field_array_ptr = mxGetField(prhs[1],0,"num_state");
  num_state = (int)mxGetScalar(field_array_ptr);

  field_array_ptr = mxGetField(prhs[1],0,"block_len");
  block_len = (int)mxGetScalar(field_array_ptr);

  field_array_ptr = mxGetField(prhs[1],0,"word_len");
  n = (int)mxGetScalar(field_array_ptr);

  field_array_ptr = mxGetField(prhs[1],0,"term");
  known_end = (int)mxGetScalar(field_array_ptr);



/*-------------------------------------------------------------------------------------
 * OUTPUT_DECLARATIONS
 *------------------------------------------------------------------------------------- */
  plhs[0] = mxCreateDoubleMatrix(block_len,1,mxREAL);           /* L_info */
  L_info  = mxGetPr(plhs[0]);

  if (nlhs==2) {
     plhs[1] = mxCreateDoubleMatrix(block_len*n,1,mxREAL);      /* L_code */
     L_code  = mxGetPr(plhs[1]);
  }


/*-------------------------------------------------------------------------------------
 *  Initialize probabilities
 *------------------------------------------------------------------------------------- */
 
  num_word = 1<<n;

  if ((metric = calloc(num_word, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer metric!\n");
     exit(1);
  }


  if ((gamma0 = calloc(num_state*block_len, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer gamma0!\n");
     exit(1);
  }

 if ((gamma1 = calloc(num_state*block_len, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer gamma1!\n");
     exit(1);
  }


  if ((lambda0 = calloc(num_state*block_len, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer lambda0!\n");
     exit(1);
  }


  if ((lambda1 = calloc(num_state*block_len, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer lambda1!\n");
     exit(1);
  }


  if ((alpha = calloc(num_state*(block_len+1), sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer alpha!\n");
     exit(1);
  }

  if ((beta = calloc(num_state*block_len, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer beta!\n");
     exit(1);
  }


  if ((norm = calloc(block_len, sizeof(double))) == NULL) {
     fprintf(stderr, "Nicht genuegend Speicher fuer norm!\n");
     exit(1);
  }


 

  for (m = 0; m < num_state; m++) 
    for (k = 0; k <= block_len; k++) 
        *(alpha+m*(block_len+1)+k) = 0.0;

  *alpha = 1.0;


/*-------------------------------------------------------------------------------------
 * End Initialization
 *------------------------------------------------------------------------------------- */


/*-------------------------------------------------------------------------------------
 * RUN_OUT_CODE:
 *------------------------------------------------------------------------------------- */

  for (k = 0; k < block_len; k++) {
      norm[k] = 0.0;

/*------------------------------------------------------------------------
 * Calculating channel transition probabilities for current time instant
 *------------------------------------------------------------------------ */
      for (i=0; i<num_word; i++) {
          m = i;
          metric[i] = 1.0;
          for (j=0; j<n; j++) {
               metric[i] /= (1 + exp( (m & 1) ? *(sig+j) : -*(sig+j) ) );
               m >>= 1;
          }
      }

/*------------------------------------------------------------------------
 * loop over all states
 *------------------------------------------------------------------------ */
      for (m=0; m<num_state; m++) {
        *(gamma0+m*block_len+k) = metric[(int)trellis_out[m]] / (1 + exp(-(*Z)));
        dummy = *(gamma0+m*block_len+k) * *(alpha+m*(block_len+1)+k);
        *(alpha+((int)trellis_next[m])*(block_len+1)+k+1) += dummy;
        norm[k] += dummy;
 
        *(gamma1+m*block_len+k) = metric[(int)trellis_out[m+num_state]] / (1 + exp(*Z));
        dummy = *(gamma1+m*block_len+k) * *(alpha+m*(block_len+1)+k);
        *(alpha+((int)trellis_next[m+num_state])*(block_len+1)+k+1) += dummy;
        norm[k] += dummy;
      }

      for (m = 0; m < num_state; m++) {
        *(alpha+m*(block_len+1)+k+1) /= norm[k];
      }
      sig += n;
      Z++;

  } /* end for (k ... */



/*---------------------------------------------------------------------------------------------- 
 * End of data block has been reached, starting backward recursion
 *----------------------------------------------------------------------------------------------*/

  if (known_end == 1) {
     for (m = 1; m < num_state; m++) {
       *(beta+m*block_len+block_len-1) = 0.0;
     }
     *(beta+last_state*block_len+block_len-1) = 1.0;
  }
  else {
     for (m = 0; m < num_state; m++) {
         *(beta+m*block_len+block_len-1) = *(alpha+m*(block_len+1)+block_len);
     }
  }


/*--------------------------------------------------------------------------------------------
 * Calculating backwards recursion (beta)
 *--------------------------------------------------------------------------------------------*/
  for (k = block_len-1; k > 0; k--) {
      for (m = 0; m < num_state; m++) {
          *(beta+m*block_len+k-1) = *(gamma0+m*block_len+k)* *(beta+((int)trellis_next[m])*block_len+k)            / norm[k];
          *(beta+m*block_len+k-1) += *(gamma1+m*block_len+k)* *(beta+((int)trellis_next[m+num_state])*block_len+k) / norm[k];
      }
  }


/*---------------------------------------------------------------------------------------------
 * Calculating log-likelihood-ratio for each time instance k
 *---------------------------------------------------------------------------------------------*/
  for (k = 0; k < block_len; k++) {


/*---------------------------------------------------------------------------------------------
 * Calculating log-likelihood-ratio for infobits for each time instant k 
 * 
 * For each infobit the LLR`s to be 1 or 0 are evaluated by the probibilities of the
 * states at time instant k
 * e.g. lambda_j(m,k) or alpha_j(m,k+1) or beta(m,k) 
*---------------------------------------------------------------------------------------------*/

      L0 = 0.0;
      L1 = 0.0;

      for (m = 0; m < num_state; m++) {      
         L0 += *(alpha+m*(block_len+1)+k) * *(gamma0+m*block_len+k) * *(beta+((int)trellis_next[m])*block_len+k);
         L1 += *(alpha+m*(block_len+1)+k) * *(gamma1+m*block_len+k) * *(beta+((int)trellis_next[m+num_state])*block_len+k);
      }
      L_info[k] = L0 / L1;
      if (L_info[k]>1.0e200)
          L_info[k]=1000.0;
      else
      if (L_info[k]<1.0e-200)
          L_info[k]=-1000.0;
      else
          L_info[k] = log( L_info[k] );
      


      if (nlhs==2) {



/*---------------------------------------------------------------------------------------------
 * Calculating log-likelihood-ratio for codebits for each time instant k
 * 
 * For each codebit the LLR`s to be 1 or 0 are evaluated by the probilities of the
 * corresponding transitions. This means:
 * 1. at each (foregoing state m) we first look if the codebit is 0 or 1 for the
 *    infobit 0 (out_info0 & 1)==0 or ==1
 * 2. at each (foregoing state m) we then look if the codebit is 0 or 1 for the
 *    infobit 1 (out_info1 & 1)==0 or ==1
 * 3. If codebit=0 or 1 we calculate the possibility of the transsition L0 or L1 by alpha0 and
 *    alpha1 as well as gammai (i=0,1) and beta
 *
 * Due to foregoing states: alpha_j(m,k) or beta(next state,k) and gamma_j(m,k) are used
 *---------------------------------------------------------------------------------------------*/

      for (j=0;j<n;j++) {

         L0 = 0.0;
         L1 = 0.0;

         for (m = 0; m < num_state; m++) {

             out_info0 = ((int)trellis_out[m])>>j;
             out_info1 = ((int)trellis_out[m+num_state])>>j;

         
         if ((out_info0 & 1) == 0)     /*Infobit=0 Codebit=0*/
                L0 += *(alpha+m*(block_len+1)+k) * *(gamma0+m*block_len+k) * *(beta+((int)trellis_next[m])*block_len+k);

             else                          /*Infobit=0 Codebit=1*/
                L1 += *(alpha+m*(block_len+1)+k) * *(gamma0+m*block_len+k) * *(beta+((int)trellis_next[m])*block_len+k);

         if ((out_info1 & 1) == 0)     /*Infobit=1 Codebit=0*/
                L0 += *(alpha+m*(block_len+1)+k) * *(gamma1+m*block_len+k) * *(beta+((int)trellis_next[m+num_state])*block_len+k);

             else                          /*Infobit=1 Codebit=1*/
                L1 += *(alpha+m*(block_len+1)+k) * *(gamma1+m*block_len+k) * *(beta+((int)trellis_next[m+num_state])*block_len+k);
         }
         L_code[(k*n)+j] = L0 / L1;

     if (L_code[(k*n)+j] > 1.0e200)
             L_code[(k*n)+j] = 1000.0;
         else
         if (L_code[(k*n)+j] < 1.0e-200)
             L_code[(k*n)+j] = -1000.0;
         else
             L_code[(k*n)+j] = log( L_code[(k*n)+j] );
      }   /* for j= */
    }  /* if nlhs==2 */
  }/*End of time instant k*/

  free(alpha);
  free(beta);
  free(gamma0);
  free(gamma1);
  free(metric);
  free(lambda0);
  free(lambda1);
  free(norm);
 }
