%-*- matlab -*-
%
% m-file for soft-output decoding of single parity check codes (SPC) using the MAP algorithm
% if r is a matrix, decoding is performed for each column
%
% ---------------------   input variables   -------------------------
% LLR:    matrix containing log likelihood values of each received bit
% L_a:    matrix of a priori information for each information bit in r
% approx: appox==1 indicates the usage of approximation in line 33
%
% ---------------------   output variables   -------------------------
% L:   decoding result; matrix of log likelihood values ln[P(r|d_k=1)/P(r|d_k=0)]
% L_e: matrix of extrinsic information (part of L that is independant of information bit)
%
function [L,L_e] = map_spc(LLR,L_a,approx)


[k1,k2] = size(L_a);   % Dimension of information matrix

if (k1 == 1)
   LLR = LLR(:);
   L_a = L_a(:);
end

LLR(1:k1,1:k2) = LLR(1:k1,1:k2) + L_a;   % nach E. Offer: Definition of L(r;c)

L_e = zeros(k1,k2);

for i=1:k1
  LLR_e = LLR;
  LLR_e(i,:)=[];     % log likelihood values
  
  if (approx == 1)
     L_e(i,:) = - prod(sign(LLR_e)) .* min(abs(LLR_e));    
  else
    tmp = 1 - prod(tanh(LLR_e/2));
    ptr = find(tmp == 0);
    tmp(ptr) = 1e-20;
    L_e(i,:) = log(tmp);

    tmp = 1 + prod(tanh(LLR_e/2));
    ptr = find(tmp==0);
    tmp(ptr) = 1e-20;
    L_e(i,:) = L_e(i,:) - log (tmp);
  end
end

L = L_e + LLR(1:k1,1:k2);
