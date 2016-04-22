% ##############################################################################
% ## fm_mod.m : {FM}-{Mod}ulator                                              ##
% ##############################################################################
%
% Aufruf:      [t,xfm] = fm_mod(t,v,dF,f0,phi0);
%              benoetigt integr.m
%
% Eingabe:     t    = Zeitvektor
%              v    = NF-Signalvektor
%              dF   = Frequenzhub in kHz
%              f0   = Traegerfrequenz in kHz
%              phi0 = Anfangsphase
%
% Ausgabe:     t    = Zeitvektor
%              xfm  = FM-moduliertes Bandpass-Signal
%
% Anmerkung:   Bei f0 = 0 wird die komplexe Einhuellende erzeugt
%
% K.D. Kammeyer    17-jan-97

function  [t,xfm] = fm_mod(t,v,dF,f0,phi0)

[t,vint]=integr(t,v);
if f0==0
  xfm=exp(j*(dF*2*pi.*vint+phi0));         % komplexe Einhuellende
else
  xfm=cos(2*pi*f0.*t+dF*2*pi.*vint+phi0);  % Bandpasssignal
end;

% ### EOF ######################################################################
