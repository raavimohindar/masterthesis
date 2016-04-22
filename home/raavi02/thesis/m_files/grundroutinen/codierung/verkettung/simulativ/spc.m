% ##############################################################################
% ##  spc.m : Berechnung des Single Parity Check Bits p des Eingangsvektors u ##
% ##############################################################################
%
%  function c = spc(u, parity)
% ------------------------------------------------------------------------------
% EINGABE:
%        u     :  Enthaelt nur Elemente 0/1
%                 Falls Matrix, dann wird fuer jede Zeile und Spalte Parity
%                 Bit angehaengt
%                 (product code)
%        parity:  1, p wird mit gerader Paritaet generiert
%                 0, p wird mit ungerader Paritaet generiert
% ------------------------------------------------------------------------------

function c = spc(u, parity)

[row,col] = size(u);

if (row==1)   % Zeilenvektor in Spaltenformat bringen
  u = u(:);
end

p = rem(sum(u),2);   % Parity Bit fuer jede Spalte

if (parity == 0)     % anhaengen der Parity Bit spaltenweise
  c = [u; p];
else
  c = [u; ones(size(p))-p];
end

if (col > 1 & row > 1)   % Parity Bit fuer jede Zeile
  p = rem(sum(c'),2)';

  if (parity == 0)
    c = [c p];
  else
    c = [c ones(size(p))-p];
  end

end

% ### EOF ######################################################################
