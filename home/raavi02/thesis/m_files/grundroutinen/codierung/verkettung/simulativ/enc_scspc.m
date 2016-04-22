% ##############################################################################
% ##  enc_scspc.m : Codierung mit seriell verk. SPC Codes                     ##
% ##############################################################################
%
% function c2 = enc_scspc(u, k1, k2, Pi)
% ------------------------------------------------------------------------------
% EINGABE:
%   u:  binaerer Spaltenvektor mit k1*k2 Infobit
%   k1: Anzahl der Infobit von Code 1 (horizontal)
%   k2: Anzahl der Infobit von Code 2 (vertikal)
%   Pi: Interleaver, Spaltenvektor der Laenge (k1+1)*k2
%       enthaelt die Positioinen der Bit nach dem Interleaven
%       (z.B.  1 4 7 2 5 8 3 6 9)
%
% AUSAGEB:
%   c2: binaerer Spaltenvektor mit (k1+1)(k2+1) Bit (0,1)
%-------------------------------------------------------------------------
function c2 = enc_scspc(u, k1, k2, Pi)

if (length(Pi)~=(k1+1)*k2)
  error('enc_scspc: Interleaverlaenge entspricht nicht den Groessen k1 und k2');
end
u = reshape(u(:),k1,k2);


p = rem(sum(u),2);               % Paritybit fuer aeusseren Code

c1 = [u; p];

c1 = reshape(c1(Pi),k2,k1+1);    % Interleaven fuer inneren Code

p = rem(sum(c1),2);              % Paritybit fuer inneren Code

c2 = [c1; p];
c2 = c2(:);

% ### EOF ######################################################################
