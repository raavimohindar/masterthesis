% ##############################################################################
% ##  rschien.m: Findet die Fehlerstellen eines Reed-Solomon Codewortes       ##
% ##############################################################################
%
% function [pos_err, err] = rschien(SIGMA, ERROR, K, M)
%
% Eingabe sigma, error : Fehlerstellenpolynom nach Funktion rselp.m
%         k : Laenge der Nachricht
%         m : 2^m-1 ist die Dimension des GF(2^m)
%
%       extracted from
%       Wes Wang 8/11/94, 10/11/95.
%       Copyright (c) 1996-98 by The MathWorks, Inc.
%       $Revision: 1.5 $
function [pos_err, err] = rschien(sigma, error, k, m);

if nargin < 4
  error('Not enough input parameters.')
end;

n  = 2^m-1;
tp = gftuple([-1:n-1]',m,2);    % Liste aller m-Tuple des GF(2^m)

t2 = n - k;
t  = floor(t2 / 2);

tp_num = tp * 2.^[0:m-1]';
tp_inv(tp_num+1) = 0:n;

err = error;

% Loesungen des Polynoms
loc_err = zeros(1, n) - Inf;
num_err = length(sigma) - 1;
if num_err > t
  % fuer den Fall dass es mehr als korrigierbare Fehler gibt
  err = 1;
end;
if (~err) & (num_err > 0)
  cnt_err = 0;
  pos_err = [];
  er_i = 0;
  while (cnt_err < num_err) & (er_i < n * m)
    test_flag = sigma(1);
    for er_j = 1 : num_err
      if sigma(er_j + 1) >= 0
        tmp = er_i * er_j;
        if (tmp < 0) | (sigma(er_j+1) < 0)
          tmp = -1;
        else
          tmp = rem(tmp + sigma(er_j + 1), n);
        end;

        test_flag = gfplus(test_flag, tmp, tp_num, tp_inv);
      end;
    end;
    if test_flag < 0
      cnt_err = cnt_err + 1;
      pos_err = [pos_err, rem(n-er_i, n)];
    end;
    er_i = er_i + 1;
  end;
  pos_err = rem(n+pos_err, n);
  pos_err = pos_err + 1; % wegen x^0=1 eins nach rechts ruecken
  err = num_err;
else
  if err
    err = -1;
  end;
end;

% ### EOF ######################################################################
