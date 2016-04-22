% ##############################################################################
% ##  rselp.m: Berechnung des Fehlerstellenpolynoms bei Reed-Solomon Codes    ##
% ##############################################################################
%
% function [sigma, error] = rselp(syndrome, k, m);
%
% Eingabe: Syndrom Vektor, den die Funktion 'rssyndrom' liefert
%          k : Laenge der Nachricht
%          m : 2^m-1 ist die Dimension des GF(2^m)
% Ausgabe  sigma: Fehlerstellenpolynom mit Grad t;
%
% benoetigt gftuple.m, errlocp.m (communications toolbox)
%
%      extracted from
%      Wes Wang 8/11/94, 10/11/95.
%      Copyright (c) 1996-98 by The MathWorks, Inc.
%      $Revision: 1.5 $

function [sigma, error] = rselp(syndrome, k, m);

if nargin < 3
  error('Not enough input parameters.')
end;

n  = 2^m-1;
tp = gftuple([-1:n-1]',m,2);    % Liste aller m-Tuple des GF(2^m)


t2 = n - k;
t = floor(t2 / 2);

% For non-binary case, using Berlekamp's iterative method to do the computation.

[sigma, error] = errlocp(syndrome, t, tp, n, 0, 1);

% ### EOF ######################################################################
