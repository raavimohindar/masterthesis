% ##############################################################################
% ##  turbo_scspc.m : Iterative Decodierung von seriell verketteten           ##
% ##                  SPC Codes (MAP)                                         ##
% ##############################################################################
%
%function u_hat = turbo_scspc(LLR,k1,k2,iterations,Pi,approx)
%
% ------------------------------------------------------------------------------
% EINGABE:
%        LLR:        Spaltenvektor der Empfangsdaten (LLRs)
%        k1:         Anzahl  der Informationsbit des aeusseren Codes
%        k2:         Anzahl  der Informationsbit des inneren Codes
%        iterations: Anzahl der Decodieriterationen
%        Pi:         Spaltenvektor der die Permutation von c1 enthaelt
%                    (Codebit des aeusseren Codes)
%        approx:     approx==1 Approximation (Max-Log-MAP) wird benutzt
%                    (ansonsten exakte SS-MAP)
%
% AUSGABE:
%        u_hat:      Matrix der decodierten Daten
%
% ANMERKUNGEN:
%   - benoetigt Datei: map_spc.m
% ------------------------------------------------------------------------------
function u_hat = turbo_scspc(LLR,k1,k2,iterations,Pi,approx)


LLR   = reshape(LLR(:),k2+1,k1+1);

L_a   = zeros(k2,k1+1);
u_hat = zeros(k1*k2,iterations);

for it=1:iterations
  [L,L_e]     = map_spc(LLR, L_a, approx);         % innere Decodierung

  L_a         = zeros(k1,k2);
  sig(Pi)     = L(1:k2,:);                         % De-Interleaving
  sig         = reshape(sig,k1+1,k2);
  [L,L_e]     = map_spc(sig,L_a, approx);          % aeussere Decodierung  

  L_a         = reshape(L_e(Pi),k2,k1+1);          % Interleaving der a-priori
  %                                                  Information

  tmp         = L(1:k1,1:k2);
  u_hat(:,it) = (1-sign(tmp(:))) / 2;
end

% ### EOF ######################################################################
