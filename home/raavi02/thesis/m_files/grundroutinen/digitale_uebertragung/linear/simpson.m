% ##############################################################################
% ##  simpson.m : Integration nach der Simpson-Methode                        ##
% ##              (quadr. Approximation)                                      ##
% ##############################################################################
%
% function [y_area] = simpson (y_mat, [x_diff])
% ----------------------------------------------------------------------
% Eingabe:
% y_mat   : Vektor der Laenge N oder (N x M) Matrix
%           In beiden Faellen sollte N ungerade sein, ansonsten
%           wird der letzte Funktionswert nicht beruecksichtigt
% [x_diff]: Optionales Argument, das die Abtastrate zwischen
%           zwei Spalten von y_mat festlegt
%           Der Standardwert ist x_diff = 1.0
% ----------------------------------------------------------------------
% Ausgabe:
% ----------------------------------------------------------------------
% y_area:  Wenn y_mat ein Vektor ist, ist y_area das Ergebnis der
%          Integration nach der Simpson methode
%          Wenn y_mat eine (N x M) Matrix ist, stehen im M-reihigen Vektor
%          y_area die Ergebnisse der Integration ueber die Spalten
% ----------------------------------------------------------------------
% Beispiele:
% ----------------------------------------------------------------------
% \int_0^1     1  dt != 1.0:   simpson(ones(101,1),0.01) = 1.0
% \int_0^pi sin(x)dx != 2.0:   simpson(sin(pi*(0:0.01:1)),pi/100) = 2.0
% ----------------------------------------------------------------------
% AUTHOR:   Dieter Boss,   18-jun-96
%           Dirk Nikolai,  08-nov-97
% ----------------------------------------------------------------------

function [y_area] = simpson (y_mat, x_diff);

if nargin < 2, x_diff = 1.0; end;

y_size = size(y_mat);
if min(y_size)==1             % if y_mat is a vector
  y_mat = y_mat(:);           % force column vector
  Ly = max(y_size);
else                          % if y_mat is a matrix
  Ly = y_size(1);
end;

if rem(Ly,2) == 0,            % if even length
  disp(['WARNING (simpson.m):',...
        '  Discarding last sample(s) due to even length input.']);
  y_mat(Ly,:) = [];           % force odd length
  Ly = Ly - 1;                % force odd value of Ly
end;

y0     = y_mat(1,:);
y2N    = y_mat(Ly,:);
y_even = y_mat(3:2:Ly-1,:);   % length (Ly-3)/2, where Ly-3 is even
y_odd  = y_mat(2:2:Ly,:);     % length (Ly-1)/2, where Ly-1 is even
y_area = y0 + y2N + 2*sum(y_even) + 4*sum(y_odd);
y_area = (x_diff/3)*y_area;

% ### EOF ######################################################################
