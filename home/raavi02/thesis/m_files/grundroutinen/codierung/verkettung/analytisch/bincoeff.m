% ##############################################################################
% ##  bincoeff.m : Berechnet die Binomialkoeffizienten                        ##
% ##############################################################################
%
% function b = bincoeff (n, k)
% ------------------------------------------------------------------------------
% EINGABE:
%      n (Skalar)
%      k (Skalar oder Vektor)
%
% AUSGABE:
%   b: Binomialkoeffizient   (n)
%      (Spaltenvektor)       (k)
% BEISPIEL:
%    » bincoeff(4,2);
%    ans =
%       6
%
% ANMERKUNGEN:
%   - groessere Binomialkoeffizienten werden mit dem Pascal`schen
%     Dreieck berechnet => Die Angabe von mehreren Werten fuer k
%                          als Vektor und einmaliger Aufruf von
%                          bincoeff spart Rechenzeit.
%
% AUTOR: Juergen Rinas,  31.05.1999
% ------------------------------------------------------------------------------

function b = bincoeff (n, k)

if (nargin ~= 2)
  error ('bincoeff (n, k)');
end;

if (n>60)
  % Palcal´sches Dreieck
  b=[1;1];
  for nloop=2:n; 
    b=conv(b,[1 1].'); 
  end;
  b=b(k+1);
else
  b=gamma(n+1)./(gamma(k+1).*gamma(n-k+1));
  b=b(:);
end;

% ### EOF ######################################################################
