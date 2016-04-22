/*
% --------------------------------------------------------------------------------------------
% [signal_out,v_state,v_ipr] = tvc(signal_in,offset,v_state,v_doppler,const_taps,n_tap,t_delay);
% --------------------------------------------------------------------------------------------
%
%
% (c) 2000 Heiko Schmidt, Department of Communications Engineering, University of Bremen
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

  double    *input_real;             /* real input vector */
  double    *input_imag;             /* imag input vector */

  int        sig_len;      	     /* output length */
  int        offset;      	     /* time offset */

  double    *state_real;             /* real echo state */
  double    *state_imag;             /* imag echo state */

  double    *doppler_real;           /* real doppler state */
  double    *doppler_imag;           /* imag doppler state */
  
  int        const_taps;      	     /* number of constant taps (for fast computation) */
  int	     n_tap;		     /* number of taps */

  double    *t_delay;		     /* tap delay */


  double    *output_real;	     /* output vector */
  double    *output_imag;	     /* output vector */

  double    *outstate_real;          /* real new echo state */
  double    *outstate_imag;          /* imag new echo state */

  double    *v_ipr_real;          /* real new impulse response */
  double    *v_ipr_imag;          /* imag new impulse response */

 
  
/*-------------------------------------------------------------------------------------
 * locale declarations:
 *------------------------------------------------------------------------------------- */


  int		lauf,index,n_states,lt,li,block,n_block,sig_pos,n_vipr;
  
  double	tmp;
  
  double    *ipr_real;	             /* current impulse response real part */
  double    *ipr_imag;	             /* current impulse response imag part */
  
  double    *signal_real;	     /* current signal real part */
  double    *signal_imag;	     /* current signal imag part */
  
  double    *newstate_real;	     /* current states real part */
  double    *newstate_imag;	     /* current states imag part */

  int       *zero_mark;	             /* marker to non used taps */

/*--------------------------- Input Declaration --------------------------------*/

  input_real  	= mxGetPr(prhs[0]);
  input_imag 	= mxGetPi(prhs[0]); 


  offset        = (int)mxGetScalar(prhs[1]);
  sig_len       = mxGetN(prhs[0]) * mxGetM(prhs[0]) - offset;
  /* sig_len    = (int)mxGetScalar(prhs[1]); */
  

  state_real  	= mxGetPr(prhs[2]);
  state_imag 	= mxGetPi(prhs[2]);
  
  doppler_real  = mxGetPr(prhs[3]);
  doppler_imag 	= mxGetPi(prhs[3]);
  
  const_taps    = (int)mxGetScalar(prhs[4]);
  n_tap         = (int)mxGetScalar(prhs[5]);

  t_delay       = mxGetPr(prhs[6]);

  
/*-------------------------- Output Declaration ----------------------------*/

  plhs[0]     = mxCreateDoubleMatrix(sig_len,1,mxCOMPLEX);
  output_real = mxGetPr(plhs[0]);
  output_imag = mxGetPi(plhs[0]); 

  n_states	= n_tap;
  plhs[1]       = mxCreateDoubleMatrix(n_states,1,mxCOMPLEX);
  outstate_real = mxGetPr(plhs[1]);
  outstate_imag = mxGetPi(plhs[1]); 

  if (const_taps >= 1) {
      n_block     = ceil((double)sig_len / const_taps);
  }
  else {
      n_block     = 1;
  }


  n_vipr	= (offset+1) * n_block;
  plhs[2]       = mxCreateDoubleMatrix(n_vipr,1,mxCOMPLEX);
  v_ipr_real    = mxGetPr(plhs[2]);
  v_ipr_imag    = mxGetPi(plhs[2]);
  n_vipr        = offset+1;



/*------------------------------ Pre Allocations -----------------------------*/



if ((zero_mark = calloc(offset+1, sizeof(int))) == NULL) {
     fprintf(stderr, "out of memory: zero_mark\n");
     exit(1);
}

if ((ipr_real = calloc(offset+1, sizeof(double))) == NULL) {
     fprintf(stderr, "out of memory: ipr_real\n");
     exit(1);
}

if ((ipr_imag = calloc(offset+1, sizeof(double))) == NULL) {
     fprintf(stderr, "out of memory: ipr_imag\n");
     exit(1);
}

if ((signal_real = calloc(sig_len, sizeof(double))) == NULL) {
     fprintf(stderr, "out of memory: signal_real\n");
     exit(1);
}

if ((signal_imag = calloc(sig_len, sizeof(double))) == NULL) {
     fprintf(stderr, "out of memory: signal_imag\n");
     exit(1);
}

if ((newstate_real = calloc(n_states, sizeof(double))) == NULL) {
     fprintf(stderr, "out of memory: newstate_real\n");
     exit(1);
}

if ((newstate_imag = calloc(n_states, sizeof(double))) == NULL) {
     fprintf(stderr, "out of memory: newstate_imag\n");
     exit(1);
}

  for(lauf=0;lauf<n_states;lauf++){
		newstate_real[lauf] = state_real[lauf];
		if (mxIsComplex(prhs[2])>0) {
			newstate_imag[lauf] = state_imag[lauf]; }
		else {
			newstate_imag[lauf] = 0;  
		}
	}

/* ------------------------------------------------------------------------------- */


 /* fprintf(stderr, "sig_len %d      const_taps %d   n_tap %d \n",sig_len,const_taps, n_tap); */

/* ------------------------- START CONVOLUTION -------------------------- */

if (mxIsComplex(prhs[3]) ==0) {
	const_taps = 0;
	}

if (const_taps >= 1) {

    for(block=0;block<n_block;block++){

 	/*fprintf(stderr, "n_block %d     sig_len %d      const_taps %d \n",n_block,sig_len,const_taps);*/
	for(li=0;li<=offset;li++){
		ipr_real[li] = 0;
		ipr_imag[li] = 0;
		
		zero_mark[li] = 0;
	}
	for(lt=0;lt<n_tap;lt++){
			index = (int)t_delay[lt];
			ipr_real[index] += newstate_real[lt];
			ipr_imag[index] += newstate_imag[lt];
			zero_mark[index] = 1;
								
			tmp               = newstate_real[lt] * doppler_real[lt] - newstate_imag[lt] * doppler_imag[lt];
			newstate_imag[lt] = newstate_real[lt] * doppler_imag[lt] + newstate_imag[lt] * doppler_real[lt];
			newstate_real[lt] = tmp;
	}

	for(li=0;li<=offset;li++){
		v_ipr_real[li+(n_vipr*block)] = ipr_real[li];
		v_ipr_imag[li+(n_vipr*block)] = ipr_imag[li];
	}
 
	for(lauf=0;lauf<const_taps;lauf++) {

		sig_pos = lauf + const_taps*block;
		if (sig_pos < sig_len) {
			signal_real[sig_pos] = 0;
			signal_imag[sig_pos] = 0;
			for(li=0;li<=offset;li++){
				index = sig_pos + offset - li;
				if (zero_mark[li] == 1) {
 				 if (mxIsComplex(prhs[0])>0) {
					signal_real[sig_pos] += ((ipr_real[li] * input_real[index]) - (ipr_imag[li] * input_imag[index]));
					signal_imag[sig_pos] += ((ipr_real[li] * input_imag[index]) + (ipr_imag[li] * input_real[index])); }
				 else {
					signal_real[sig_pos] += (ipr_real[li] * input_real[index]);
					signal_imag[sig_pos] += (ipr_imag[li] * input_real[index]); 
				 }
				}
			} 
			output_real[sig_pos] = signal_real[sig_pos];
			output_imag[sig_pos] = signal_imag[sig_pos];
		}
	}
 }
}	


else {

     for(li=0;li<=offset;li++){
	ipr_real[li] = 0;
	ipr_imag[li] = 0;
	zero_mark[li] = 0;

	}
     for(lt=0;lt<n_tap;lt++){
		index = (int)t_delay[lt];
		ipr_real[index] += newstate_real[lt];
		ipr_imag[index] += newstate_imag[lt]; 
		zero_mark[index] = 1;
     }

     for(li=0;li<=offset;li++){
		v_ipr_real[li] = ipr_real[li];
		v_ipr_imag[li] = ipr_imag[li];
     }

 
     for(lauf=0;lauf<sig_len;lauf++) {

	signal_real[lauf] = 0;
	signal_imag[lauf] = 0;
	
	for(li=0;li<=offset;li++){
		index = lauf + offset - li;
		if (zero_mark[li] == 1) {
		 if (mxIsComplex(prhs[0])>0) {
			signal_real[lauf] += ((ipr_real[li] * input_real[index]) - (ipr_imag[li] * input_imag[index]));
			signal_imag[lauf] += ((ipr_real[li] * input_imag[index]) + (ipr_imag[li] * input_real[index])); }
		 else {
			signal_real[lauf] += (ipr_real[li] * input_real[index]);
			signal_imag[lauf] += (ipr_imag[li] * input_real[index]); 
		 }
		}
			
	} 
	output_real[lauf] = signal_real[lauf];
	output_imag[lauf] = signal_imag[lauf];
     }
}


for(lauf=0;lauf<n_states;lauf++){
	outstate_real[lauf] = newstate_real[lauf];
	outstate_imag[lauf] = newstate_imag[lauf]; 
}
	
/* --------------------- END OF PROGRAM ------------------------ */
free (ipr_real);
free (ipr_imag);

free (signal_real);
free (signal_imag);

free (newstate_real);
free (newstate_imag);

free (zero_mark); 

 }


