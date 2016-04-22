%-*- matlab -*-
%
% m-file for soft-output decoding of linear block codes using BCJR algorithm (E. Offer)
%
% function [L,L_e] = bcjr_bc(r, L_a, H, snr, dual, info)
%
% used functions: make_trellis_bc.m
%
% ------------------------------ input variables ------------------------------
% r    : vector of received data
% L_a  : vector containing a-priori information for each information bit
% H    : parity check matrix of block code
% snr  : signal to noise ratio of the channel Es/N0 (linear, NOT in dB)
% dual : dual==1 indicates decoding with trellis of dual code, otherwise decoding with original code
%       (dual==1: H must be the generator matrix of the original code, that is the parity check matrix for the dual code)
% info : if info==1, soft-output is calculated only for information bits 1...k
%        if info==0, soft-output is calculated for all coded bits 1...n
%
% ------------------------------ output variables ------------------------------
% L:   vector containing log likelihood estimates for each information bit
% L_e: vector containing extrinsic information for each information bit
%
function [L,L_e] = bcjr_bc(r, L_a, H, snr, dual, info)


[m,n] = size(H);

if (dual == 1)           % k indicates number of information bits in a code word
   k = m;                % H is generator matrix of original code --> m = number of rows != k
else
   k = n - m;            % H is parity check matrix of original code --> m describes number of parity bits
end

if (info == 1)
    n_out = k;           % soft-output only for information bits --> n_out = k;
else
    n_out = n;           % soft-output for all coded bits --> n_out = n;
end

r = r(:);
L_a = L_a(:);
n_La = length(L_a);

trellis = make_trellis_bc(H,0)+1;              % generating trellis with parity check matrix


% ---------------------------- calculation of log-likelihood values from received data -------------------------------
%
L_e = zeros(n_out,1);

if ( dual == 1 )
  metric          = 4 * snr * r;
  metric(1:n_La)  = metric(1:n_La) + L_a;
  metric          = tanh( metric ./ 2);
else
  metric           = zeros(n,2);
  metric(:,1)      = exp(-(r+1).^2 * snr);                             % calculating transition probabilities for info = 0
  metric(1:n_La,1) = metric(1:n_La,1) ./ (1 + exp(L_a));               % considering a-priori information only for infobits
  metric(:,2)      = exp(-(r-1).^2 * snr);                             % calculating transition probabilities for info = 1
  metric(1:n_La,2) = metric(1:n_La,2) .* exp(L_a) ./ (1 + exp(L_a));   % considering a-priori information only for infobits
end

% ------------------ initialize alpha ------------------------
alpha = zeros(2^m,n+1);
alpha(1,1) = 1.0;                           % trellis starts in all-zero-state

% ------------------ initialize beta ------------------------
beta = zeros(2^m,n);
beta(1,n) = 1.0;                           % trellis ends in all-zero-state




% -------------------------------------------------------------------------------------------------------
% --------------------------------------  decoding via dual code ----------------------------------------
%
if (dual == 1) 			

% ------------------ forward recursion -----------------------
%
  for j=1:n
    alpha(:,j+1) = alpha(:,j);                                                     % alpha for info = 0
    alpha(trellis(:,j),j+1) = alpha(trellis(:,j),j+1) + alpha(:,j) .* metric(j);   % alpha for info = 1
  end

% ------------------ backward recursion -----------------------
%
  for j=n:-1:2
    beta(:,j-1) = beta(:,j);                                         % alpha for info = 0
    beta(:,j-1) = beta(:,j-1) + beta(trellis(:,j),j) .* metric(j);   % alpha for info = 1
  end

% ------------------ calculating extrinsic information -----------------------
%
  L_e0 = sum(alpha(:,1:n_out).*beta(:,1:n_out));        % info = 0
  t1 = 2^m*(0:n_out-1);                                 % avoid for ... -loop by generating trellis description
  t2 = ones(2^m,1) * t1;                                % 'trellis_neu' contains linear addresses of each element
  trellis_neu = trellis(:,1:n_out) + t2;                % --> instead of for-loop element-wise matrix-multiplication
  L_e1 = sum(alpha(:,1:n_out).*beta(trellis_neu));      % info = 1
  
  L_e = (real (log( (L_e1 + L_e0) ./ (L_e1 - L_e0) )))';

  
  
  

% -------------------------------------------------------------------------------------------------------
% --------------------------------------  decoding via original code ------------------------------------  
%
else                                  

% ------------------ forward recursion -----------------------
%
  norm = zeros(1,n);
  for j=1:n
    alpha(:,j+1) = alpha(:,j) .* metric(j,1);                                       % alpha for info = 0
    alpha(trellis(:,j),j+1) = alpha(trellis(:,j),j+1) + alpha(:,j) .* metric(j,2);  % alpha for info = 1
    norm(j) = sum(alpha(:,j+1));
    alpha(:,j+1) = alpha(:,j+1)./ norm(j);
  end


% ------------------ backward recursion -----------------------
%

  for j=n:-1:2
    beta(:,j-1) = beta(:,j) .* metric(j,1);                                     % alpha for info = 0
    beta(:,j-1) = beta(:,j-1) + beta(trellis(:,j),j) .* metric(j,2);            % alpha for info = 1
    beta(:,j-1) = beta(:,j-1) / norm(j-1);
  end

  
% ------------------ calculating extrinsic information -----------------------
%
  L_e0 = sum(alpha(:,1:n_out).*beta(:,1:n_out));                % info = 0
  t1 = 2^m*(0:n_out-1);
  t2 = ones(2^m,1) * t1;
  trellis_neu = trellis(:,1:n_out) + t2;
  L_e1 = sum(alpha(:,1:n_out).*beta(trellis_neu));          % info = 1

  L_e = (log(L_e1 ./ L_e0))';
end


% ----------------------------- Adding systematic, a priori and extrinsic information -------------------------
if (n_La > n_out)
   L = 4*snr*r(1:n_out) + L_a(1:n_out) + L_e;
else
   L = 4*snr*r(1:n_out) + [L_a(1:n_La) zeros(n_out-n_La,1)] + L_e;
end