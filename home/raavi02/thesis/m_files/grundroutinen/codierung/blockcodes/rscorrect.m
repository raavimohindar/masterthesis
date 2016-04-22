% ##############################################################################
% ##  rscorrect.m: Addiert den geschaetzten Fehler zum Empfangswort bei       ##
% ##               RS Codes                                                   ##
% ##############################################################################
%
% function codehat = rscorrect(y,loc_err,k,m);
%
% Eingabe: y : Empfangsvektor
%              y muss eines der folgenden Formate besitzen:
%              (1) eine Matrix mit 2^m-1 Reihen und M Zeilen mit 
%                  binaeren Elementen
%              (2) einen Spaltenvektor mit 2^m-1 Elementen im exponentiellen
%                  Format
%              (3) einen Spaltenvektor mit m*(2^m-1) binaeren Elementen%
%          loc_err : Fehlerstellenamplituden nach rserramp.m
%          k : Laenge der Nachricht
%          m : 2^m-1 ist die Dimension des GF(2^m)
%
%      extracted from
%      Wes Wang 8/11/94, 10/11/95.
%      Copyright (c) 1996-98 by The MathWorks, Inc.
%      $Revision: 1.5 $

function code_hat = rscorrect(y,loc_err, m);

if nargin < 3
  error('Not enough input parameters.')
end;

n      = 2^m-1;
tp     = gftuple([-1:n-1]',m,2);    % Liste aller m-Tuple des GF(2^m)
tp_num = tp * 2.^[0:m-1]';
tp_inv(tp_num+1) = 0:n;



% Konvertierung von y in exponentielles Format
[r,c] = size(y);
if (c==m)
  y = bi2de(y)-1;
elseif ((r==n*m) & (c==1))
  y = vec2mat(y,m);
  y = bi2de(y)-1;
end
if (length(y)~=n)
  error('Empfangssignal hat ungueltiges Format!')
end

code_hat = gfplus(loc_err, y', tp_num, tp_inv);
code_hat(find(code_hat==-inf)) = -1;
code_hat = code_hat(:)+1;

code_hat = de2bi(code_hat,m);
code_hat = code_hat';
code_hat = code_hat(:);

% ### EOF ######################################################################
