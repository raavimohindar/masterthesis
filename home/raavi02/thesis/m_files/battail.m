%-*- matlab -*-
%
% m-file for soft-output decoding of linear block codes using extended algorithm of Battail (E. Offer)
%
% used functions: make_trellis_bc.m
%
% ------------------------------ input variables ------------------------------
% r:    vector of received data
% L_a:  vector containing a-priori information for each information bit
% H:    parity check matrix of block code
% snr:  signal to noise ratio of the channel (linear, NOT in dB)
% dual: dual==1 indicates decoding with trellis of dual code, otherwise decoding with original code
%       (dual==1: H must be the generator matrix of the original code, that is the parity check matrix for the dual code)
%
% ------------------------------ output variables ------------------------------
% L:   vector containing log likelihood estimates for each information bit
% L_e: vector containing extrinsic information for each information bit
%
function [L,L_e] = battail(r, L_a, H, snr, dual)

[m,n] = size(H);

if (dual == 1)           % k indicates number of information bits in a code word
   k = m;                % H is generator matrix of original code --> m = number of rows != k
else
   k = n - m;            % H is parity check matrix of original code --> m describes number of parity bits
end

r = r(:);
L_a = L_a(:);

% ---------------------------- calculation of log-likelihood values from received data -------------------------------
%
L_e = zeros(k,1);

if ( dual == 1)
  metric      = 4 * snr * r;
  metric(1:k) = metric(1:k) + L_a;
  metric      = tanh( metric ./ 2);
else
  metric        = zeros(n,2);
  metric(:,1)   = exp(-(r-1).^2 * snr);                         % calculating transition probabilities for info = 0
  metric(1:k,1) = metric(1:k,1) ./ (1 + exp(-L_a));             % considering a-priori information only for infobits
  metric(:,2)   = exp(-(r+1).^2 * snr);                         % calculating transition probabilities for info = 1
  metric(1:k,2) = metric(1:k,2) ./ (1 + exp(L_a));              % considering a-priori information only for infobits
end
  
for i=1:k                                   % calculating extrinsic information for each information bit
    trellis = make_trellis_bc(H,i)+1;       % generating trellis discarding i-th column of H
    gamma = metric;
    gamma(i,:) = [];
    alpha = zeros(2^m,n);
    alpha(1,1) = 1.0;                       % trellis starts in all-zero-state
    
    
% --------------------------------------  decoding via dual code ----------------------------------------
  if (dual == 1) 			
    for j=1:n-1
       alpha(:,j+1) = alpha(:,j);                                                        % alpha for info=0
       alpha(trellis(:,j),j+1) = alpha(trellis(:,j),j+1) + alpha(:,j) .* gamma(j);   % alpha for info=1
    end

    L_e(i) = 2 * real( atanh( alpha(1,n) / alpha(int_state(H(:,i)')+1,n) ) );
%    L_e(i) = log( abs( (alpha(int_state(H(:,i)')+1,n)+alpha(1,n))/(alpha(int_state(H(:,i)')+1,n)-alpha(1,n)) ) );


% --------------------------------------  decoding via original code ------------------------------------  
  else                                  
    for j=1:n-1
        alpha(:,j+1) = alpha(:,j) .* gamma(j,1);                                        % alpha for info = 0
        alpha(trellis(:,j),j+1) = alpha(trellis(:,j),j+1) + alpha(:,j) .* gamma(j,2);   % alpha for info = 1
        norm = sum(alpha(:,j+1));
        alpha(:,j+1) =alpha(:,j+1)./ norm;
    end
    
    L_e(i) = - log(alpha(int_state(H(:,i)')+1,n)) + log(alpha(1,n));
  end
end


% -------------------------------------- adding systematic and a priori information ----------------------

L = 4*snr*r(1:k) + L_a + L_e;


