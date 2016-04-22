% ##############################################################################
% ##  bpmisch.m : {B}and{p}ass{misch}ung                                      ##
% ##############################################################################
%
% Aufruf:      [t,x] = bpmisch(t,s,f0);
%
% Eingabe:     t  = zeitliche Stuetzstellen in ms
%              s  = komplexe Einhuellende
%              f0 = Traegerfrequenz in kHz
%
% Ausgabe:     t  = zeitliche Stuetzstellen in ms
%              x  = reelles Bandpass-Signal
%
% D. Boss 15-jan-97

function [t,x] = bpmisch(t,s,f0)

fA=1/(t(2)-t(1));

Omega0 = 2*pi*f0/fA;
x = real(s(:).*exp(j*Omega0.*t(:)*fA));

if size(s,2)~=1
  x = x.';
end;

% ### EOF ######################################################################
