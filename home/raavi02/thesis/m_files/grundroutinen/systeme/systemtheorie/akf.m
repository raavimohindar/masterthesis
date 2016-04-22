% ##############################################################################
% ##  akf.m : Schaetzung der Autokorrelationsfunktion                         ##
% ##############################################################################
%
% Aufruf:    [tau,rxx]=akf(M,x)
%
% Eingabe:   M = max{tau} - 1
%            x = Signalvektor darf eine beliebige Dimension annehmen und kann
%            als Spalten- oder Zeilenvektor uebergeben werden
%
% Ausgabe:   rxx(tau) = [rxx(-(M-1)), ..., rxx(M-1)] (Spaltenvektor)
%
% Anmerkung: Nicht erwartungstreue AKF-Schaetzung durch Vektormultiplikation

function [tau,rxx]=akf(M,x)

x=x(:);
N=length(x);
r0=x'*x;
r=zeros(M-1,1);
for ii=1:M-1
  r(ii)=[zeros(1,ii) x'] * [x; zeros(ii,1)];
end
rxx=[conj(r(M-1:-1:1)); r0; r]/N;
tau=[-(M-1):M-1];

if size(x,2)~=1
  tau = tau.';
  rxx = rxx.';
end;

% ### EOF ######################################################################
