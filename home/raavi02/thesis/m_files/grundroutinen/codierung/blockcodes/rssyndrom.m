% ##############################################################################
% ##  rssyndrom.m: Berechnung des Syndroms bei Reed-Solomon Codes             ##
% ##############################################################################
%
% function syndrome = rssyndrom(y, k, m);
%
% Eingabe: y : Empfangsvektor
%              y muss eines der folgenden Formate besitzen:
%              (1) Spaltenvektor mit 2^m-1 Elementen im exponentiellen
%                  Format
%              (3) Spaltenvektor mit m*(2^m-1) binaeren Elementen
%
%          k : Laenge der Nachricht
%          m : 2^m-1 ist die Dimension des GF(2^m)

%      extracted from
%      Wes Wang 8/11/94, 10/11/95.
%      Copyright (c) 1996-98 by The MathWorks, Inc.
%      $Revision: 1.5 $

function syndrome = rssyndrom(y, k, m);

if nargin < 3
  error('Not enough input parameters.')
end;

n  = 2^m-1;
tp = gftuple([-1:n-1]',m,2);    % Liste aller m-Tuple des GF(2^m)

t2 = n - k;
t = floor(t2 / 2);

% Konvertierung von y in exponentielles Format
l = length(y);
if (l==n*m)
  y = vec2mat(y,m);
  y = bi2de(y)-1;
end
if (length(y)~=n)
  error('Empfangssignal hat ungueltiges Format!')
end

% (1) find syndrome
for sy_i = 1 : t2
  % The ith element of syndrome equals the result of
  % dividing code(X) by alpha^i+X
  [tmp, syndrome(sy_i)] = gfdeconv(y.', [sy_i, 0], tp);
end;

% ### EOF ######################################################################
