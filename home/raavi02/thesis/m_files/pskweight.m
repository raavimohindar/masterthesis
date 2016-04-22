% ##############################################################################
% ##  pskweight.m : Vektor zur Gewichtung fuer die (D)PSK Pb-Berechnung       ##
% ##                ermitteln                                                 ##
% ##############################################################################
%
% function [w] = pskweight(K)
% ------------------------------------------------------------------------------
% EINGABE:
%      K: Anzahl der Bits pro (D)PSK-Symbol
%         (Skalar)
%
% AUSGABE:
%      w: mittlere Anzahl der Bitfehler bei Phasenabweichung
%         (Spaltenvektor)
%    map: Art des Mapping
%         'gray': Gray-Codierung
%         'nat':  natuerliches Mapping
%         (string)
%
% ANMERKUNGEN:
%   - fuer gray-kodierte (D)PSK
%   - benoetigt Dateien: graytable.m, nattable.m
%
% AUTOR: Juergen Rinas,  08.12.1998
% Erweitert von Volker Kuehn, 05.01.2001
% ------------------------------------------------------------------------------

function w=pskweight(K,map)

if (nargin<2)
  map='gray'; % default: Gray-Codierugn
end

M=2^K;

if (upper(map(1))=='N')
  codetable=nattable(K);
else
  codetable=graytable(K);
end

mcoeff=zeros(M,M);

% alle Sendesymbole untersuchen
for i=1:M;
  % Codetabelle Rotieren
  rottab=[codetable(i:M,:);codetable(1:(i-1),:)];
  oksymbol=rottab(1,:); % Voraussetzung: Symbol rottab(1,:) wurde gesendet

  tabelle=xor(ones(M,1)*oksymbol,rottab);

  mcoeff(i,:)=sum(tabelle,2).';
end;

% Ueber alle als Sendesymbol vorausgesetzten Moeglichkeiten mitteln
w=sum(mcoeff,1)/M;
w=w(:); % Spaltenvektor!

% ### EOF ######################################################################
