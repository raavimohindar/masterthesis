% ##############################################################################
% ##  enc_pcspc.m : Codierung mit parallel verketteten SPC Codes              ##
% ##############################################################################
%
% function c = enc_pcspc(u, k1, k2, Pi)
% ------------------------------------------------------------------------------
% EINGABE:
%   u:  binaerer Spalten Vektor, enthaelt k1*k2 Information Bit (0,1)
%   k1: Anzahl der Infobit von Code 1 (horizontal)
%   k2: Anzahl der Infobit von Code 2 (vertikal)
%   Pi: Interleaver, Spaltenvektor der Laenge k1*k2, enthaelt die Positionen
%       der Bit nach dem Interleaven (z.B.  1 4 7 2 5 8 3 6 9)
%
% AUSGABE:
%   c: Spaltenvektor mit k1*k2+k1+k2 Bit (0,1)
%-------------------------------------------------------------------------
function c = enc_pcspc(u, k1, k2, Pi)

if (length(Pi)~=k1*k2)
  error('enc_pcspc: Interleaverlaenge entspricht nicht den Groessen k1 und k2');
end
u = u(:);

c1 = rem(sum(reshape(u,k1,k2)),2);     % Generiere Parity Bit fuer ersten Code
c2 = rem(sum(reshape(u(Pi),k2,k1)),2); % Generiere Parity Bit fuer zweiten Code

c = [u;c1(:);c2(:)];

% ### EOF ######################################################################
