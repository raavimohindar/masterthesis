% ##############################################################################
% ##  turbo_pcspc.m :Iterative Decodierung von parallel verketteten           ##
% ##                 SPC Codes (MAP)                                          ##
% ##############################################################################
%
%function u_hat = turbo_pcspc(LLR,k1,k2,iterations,Pi,approx)
%
%
% ------------------------------------------------------------------------------
% EINGABE:
%        LLR:        Matrix der Empfangsdaten (LLRs)
%        k1:         Anzahl  der Informationsbit des ersten Codes
%        k2:         Anzahl  der Informationsbit des zweiten Codes
%        iterations: Anzahl der Decodieriterationen
%        Pi:         Spaltenvektor der die Permutation von u enthaelt
%        approx:     approx==1 Approximation (Max-Log-MAP) wird benutzt
%                    (ansonsten exakte SS-MAP)
%
% AUSGABE:
%        u_hat:      Matrix der decodierten Daten
%
% ANMERKUNGEN:
%   - benoetigt Datei: map_spc.m
% ------------------------------------------------------------------------------
function u_hat = turbo_pcspc(LLR,k1,k2,iterations,Pi,approx)


LLR  = LLR(:);                          % Empfangswerte
L_u  = LLR(1:k1*k2);                    % systematischer Anteil
L_c1 = [LLR(k1*k2+1:k1*k2+k2)]';        % Redundanz vom ersten Code
L_c2 = [LLR(k1*k2+k2+1:length(LLR))]';  % Redundanz vom zweiten Code

L_e   = zeros(k1*k2,1);
u_hat = zeros(k1*k2,iterations);

for it=1:iterations
  L_a(Pi) = L_e;                            % De-interleaving
  L_a     = reshape(L_a,k1,k2);
  sig     = [reshape(L_u,k1,k2);L_c1];
  [L,L_e] = map_spc(sig, L_a, approx);      % Decodierung des ersten Codes
  L_e     = L_e(1:k1,:);

  L_a     = reshape(L_e(Pi),k2,k1);         % Interleaving and re-arranging
  sig     = [reshape(L_u(Pi),k2,k1);L_c2];
  [L,L_e] = map_spc(sig,L_a, approx);       % Decodierung des zweiten Codes
  L_e     = L_e(1:k2,:);

  tmp(Pi) = L(1:k2,1:k1);
  u_hat(:,it) = (1-sign(tmp(:))) / 2;
end

% ### EOF ######################################################################
