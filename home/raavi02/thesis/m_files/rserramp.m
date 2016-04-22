% ##############################################################################
% ##  rserramp.m: Berechnung der Fehlerstellenamplituden des Syndroms bei     ##
% ##              RS-Codes                                                    ##
% ##############################################################################
%
%  function loc_err = rserramp(syndrome, sigma, pos_err, err, k, m);
%
%  Eingabe: syndrome : Syndrom Vektor nach rssyndrom.m
%           sigma    : Fehlerstellenvektor
%           pos_err  : Fehlerpositionen nach rserrloc.m
%           err      : Flag nach rserrloc.m
%           k        : Laenge der Nachricht
%           m        : Dimension ist 2^m-1
%
%  Ausgabe: loc_err  : Amplituden der Fehlerstellen

%       extracted from
%       Wes Wang 8/11/94, 10/11/95.
%       Copyright (c) 1996-98 by The MathWorks, Inc.
%       $Revision: 1.5 $
function loc_err = rserramp(syndrome, sigma, pos_err, err, k, m);

if nargin < 6
  error('Not enough input parameters.')
end;

n  = 2^m-1;
tp = gftuple([-1:n-1]',m,2);    % Liste aller m-Tuple des GF(2^m)

t2 = n - k;
t = floor(t2 / 2);

tp_num = tp * 2.^[0:m-1]';
tp_inv(tp_num+1) = 0:n;

loc_err = zeros(1, n) - Inf;
num_err = length(sigma) - 1;

% (4) find the amplitude of the error
if (err > 0)
  % Construct Z(X) in (6.34)
  Z = 1;
  for am_i = 1 : num_err
    Z(am_i+1) = gfplus(syndrome(am_i), sigma(am_i + 1), tp_num, tp_inv);
    if am_i > 1
      for am_j = 1 : am_i - 1
        Z(am_i + 1) = gfplus(Z(am_i + 1),...
                             gfmul(sigma(am_j + 1),...
                                   syndrome(am_i - am_j), tp), tp_num, tp_inv);
      end;
    end;
  end;

  % use pos_err here.
  pos_err_1 = pos_err - 1;    % considering position starting from zero.
  er_loc = [];
  for am_i = 1 : length(pos_err_1)
    num = 0;
    den = 0;
    pos_err_inv = rem(n - pos_err_1(am_i), n);
    for am_j = 1 : num_err
      % num = gfadd(num, gfmul(Z(am_j + 1), pos_err_inv * am_j, tp), tp);
      tmp = pos_err_inv * am_j;
      if (tmp < 0) | (Z(am_j + 1) < 0)
        tmp = -1;
      else
        tmp = rem(tmp + Z(am_j + 1), n);
      end;
      num = gfplus(num, tmp, tp_num, tp_inv);
      if am_i ~= am_j
        % den = gfmul(den, gfadd(0, gfmul(pos_err_1(am_j),...
        %             pos_err_inv, tp), tp), tp);
        if (den < 0)
          den = -1;
        else
          if (pos_err_1(am_j) < 0) | (pos_err_inv < 0)
            tmp = -1;
          else
            tmp = rem(pos_err_1(am_j) + pos_err_inv, n);
          end;
          tmp = gfplus(0, tmp, tp_num, tp_inv);
          if (tmp < 0)
            den = -1;
          else
            den = rem(den + tmp, n);
          end;
        end;
      end;
      %
    end;
    %        er_loc(am_i) = gfmul(num, n-den, tp);
    tmp = n - den;
    if (tmp < 0) | (num < 0)
      er_loc(am_i) = -1;
    else
      er_loc(am_i) = rem(tmp  + num, n);
    end;
    %
  end;
  loc_err(pos_err) = er_loc;
end;

% ### EOF ######################################################################
