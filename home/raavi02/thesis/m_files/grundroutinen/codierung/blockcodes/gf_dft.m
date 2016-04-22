% ##############################################################################
% ##  gf_dft.m: diskrete Fouriertransformation auf Galoisfeldern (GF2^m)      ##
% ##############################################################################
%
% function:     A = gf_dft(a,m)
%
% Eingabe:  a  : Spaltenvektor der Laenge n=2^m-1 mit Elementen aus G(2^m) in
%                Exponentialform oder (n,m)-Matrix mit binaeren Elementen (0,1)
%           m  : Dimension des Galoisfeldes
%
% Ausgabe:  A  : Spaltenvektor der Laenge n=2^m-1 mit Elementen aus GF(2^m) in
%                Exponentialform
%-------------------------------------------------------------------------------

function A = gd_dft(a,m)
n     = 2^m-1;
tp = gftuple([-1:n-1]',m,2);    % Liste aller m-Tuple des GF(2^m)

% Konvertiere Eingangsvektor a in Exponentialform
[r,c] = size(a);
if (c>1)
  [a,index] = gftuple(a,m,2);     % a ist in Polynomform
  a = index;
elseif (r==n*m)
  a = vec2mat(a,m);
  a = bi2de(a) - 1;
end
if (length(a)~=n)
  error('Eingangsvektor hat ungueltiges Format!');
end


% Transformationsmatrix
T = zeros(n);
u = 1:n-1;
T(2:n,2:n) = rem(u'*u,n);

% Ausgangsvektor
A = -ones(n,1);
for i=1:n
  for j=1:n
    A(i) = gfadd(A(i),gfmul(T(i,j),a(j),tp),tp);
  end
end

% ### EOF ######################################################################
