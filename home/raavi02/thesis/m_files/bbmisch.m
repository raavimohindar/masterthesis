% ##############################################################################
% ##  bbmisch.m : {B}asis{b}and{misch}ung (Quadraturmischer)                  ##
% ##############################################################################
%
% Aufruf:      [t,s,xplus] = bbmisch(t,x,f0);
%
% Eingabe:     t  = zeitliche Stuetzstellen in ms
%              x  = reelles Bandpass-Signal
%              f0 = Traegerfrequenz in kHz
%
% Ausgabe:     t     = zeitliche Stuetzstellen in ms
%              s     = komplexe Einhuellende
%              xplus = analytisches Signal zu x
%
% K.D. Kammeyer, D. Nikolai   15-jan-97

function [t,s,xplus] = bbmisch(t,x,f0)

fA=1/(t(2)-t(1));

% Entfernung des Spektrums im negativen Frequenzbereich
X = fft(x(:));
X(round(length(X)/2)+1:length(X)) = zeros(fix(length(X)/2),1);

% analytisches Zeitsignal
xplus = 2*ifft(X);

% komplexe Einhuellende
Omega0 = 2*pi*f0/fA;
s = xplus.*exp(-j*Omega0.*t(:)*fA);

if size(x,2)~=1
  xplus = xplus.';
  s     = s.';
end;

% ### EOF ######################################################################
