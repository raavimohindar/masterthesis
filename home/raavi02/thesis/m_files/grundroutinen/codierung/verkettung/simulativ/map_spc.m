% ##############################################################################
% ##  map_spc.m : Soft-Output Decodierung eines SPC-Codes (MAP Algorithmus)   ##
% ##############################################################################
%
% function [L,L_e] = map_spc(LLR,L_a,approx)
% ------------------------------------------------------------------------------
% EINGABE:
%     LLR:    Matrix mit den LLRs der Empfangsdaten
%     L_a:    Matrix mit den a-priori Informationen fuer jedes Infobit
%     approx: appox==1 Naeherung in Zeile 33 wird benutzt
%
% AUSGABE
%     L:   Decodierergebnis; Matrix der LLRs ln[P(r|d_k=1)/P(r|d_k=0)]
%     L_e: Matrix der extrinsischen Information
%
function [L,L_e] = map_spc(LLR,L_a,approx)


[k1,k2] = size(L_a);   % Dimension der Informationsmatrix
[n1,n2] = size(LLR);   % Dimension der Empfangsmatrix

if (k1 == 1)
  LLR = LLR(:);
  L_a = L_a(:);
  k1  = k2;
  n1  = n2;
  k2  = 1;
  n2  = 1;
end

LLR(1:k1,1:k2) = LLR(1:k1,1:k2) + L_a;   % nach E. Offer: Definition von L(r;c)

L_e = zeros(n1,n2);

for i=1:n1
  LLR_e     = LLR;
  LLR_e(i,:)=[];     % log likelihood Werte

  if (approx == 1)
    L_e(i,:) = prod(sign(LLR_e)) .* min(abs(LLR_e));
  else
    tmp = 1 + prod(tanh(LLR_e/2));
    ptr = find(tmp == 0);
    tmp(ptr) = 1e-20;
    L_e(i,:) = log(tmp);

    tmp = 1 - prod(tanh(LLR_e/2));
    ptr = find(tmp==0);
    tmp(ptr) = 1e-20;
    L_e(i,:) = L_e(i,:) - log (tmp);
  end
end

L = L_e + LLR;

% ### EOF ######################################################################
