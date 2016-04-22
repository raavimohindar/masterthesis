% ##############################################################################
% ##  kkf.m : Schaetzung der Kreuzkorrelationsfunktion                        ##
% ##############################################################################
%
% Aufruf:    [tau,rxy]=kkf(M,x,y);
%
% Eingabe:   x,y = Signalvektoren duerfen eine beliebige gleiche Dimension
%            annehmen und koennen als Spalten- oder Zeilenvektor uebergeben
%            werden.
%            M = max{tau} - 1
%
% Ausgabe:   rxy = [rxy(-(M-1)), ..., rxy(M-1)] (Spaltenvektor)
%
% Anmerkung: Nicht erwartungstreue KKF-Schaetzung durch Vektormultiplikation

function [tau,rxy]=kkf(M,x,y)

x=x(:);
y=y(:);
N=length(x);
r0=x'*y;
r_plus=zeros(M-1,1);
for ii=1:M-1
  r_plus(ii)=[zeros(1,ii) x'] * [y; zeros(ii,1)];
end
r_minus=zeros(M-1,1);
for ii=1:M-1
  r_minus(ii)=[zeros(1,ii) y'] * [x; zeros(ii,1)];
end

rxy=[conj(r_minus(M-1:-1:1)); r0; r_plus]/N;
tau=[-(M-1):M-1];

% ### EOF ######################################################################
