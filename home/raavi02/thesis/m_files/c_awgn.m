% Function calculates channel capacity of AWGN channel for defined signal-to-noise-ratio
% and discrete input alphabet
%
% -------------------------------------------------------------------------------------
% syntax: capacity = c_awgn(alphabet,P_x,snr);
%
% -------------------------------------------------------------------------------------
% Input Parameters:
%       alphabet : vector containing signal space at input of AWGN channel
%                  average energy is assumed to be unity  (see signal_space.m)!
%       P_x      : a priori probabilities of input alphabet
%       snr      : vector containing considered signal-to-noise-ratio of AWGN channel
%                  in dB
%
% -------------------------------------------------------------------------------------
% Output Parameter:
%       capacity : vector of same size as snr containing channel capacity
% -------------------------------------------------------------------------------------
% Author: Volker Kuehn, 06.12.01
% -------------------------------------------------------------------------------------
function capacity = c_awgn(alphabet,P_x,snr)


  EsN0   = 10.^(snr/10);
  sigma2 = 1./(2*EsN0);   % variance for each, real and imaginary part of noise

  xmin = -max(abs(real(alphabet)))-10*sqrt(sigma2);
  xmax =  max(abs(real(alphabet)))+10*sqrt(sigma2);

  N_samples = 200;
  tol       = abs(xmax-xmin)/N_samples;
  x = xmin:tol:xmax;
  
  if isreal(alphabet)
%     capacity = quad(@integrand_awgn_1d,xmin,xmax,tol,[],alphabet,P_x,sigma2);
     integrand = integrand_awgn_1d(x,alphabet,P_x,sigma2);
     capacity  = sum(integrand);
  else
     ymin = -max(abs(imag(alphabet)))-10*sqrt(sigma2);
     ymax =  max(abs(imag(alphabet)))+10*sqrt(sigma2);
     y    = ymin:tol:ymax;

%     capacity = dblquad(@integrand_awgn_2d,xmin,xmax,ymin,ymax,tol,[],alphabet,P_x,sigma2);
     integrand = integrand_awgn_2d(x,y,alphabet,P_x,sigma2);
     capacity  = sum(sum(integrand));

  end


