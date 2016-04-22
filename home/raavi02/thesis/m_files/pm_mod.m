% ##############################################################################
% ## pm_mod.m : {PM}-{Mod}ulator                                              ##
% ##############################################################################
%
% Aufruf:      [t,xpm] = pm_mod(t,v,dphi,f0,phi0);
%
% Eingabe:     t    = Zeitvektor
%              v    = NF-Signalvektor
%              dphi = Phasenhub
%              f0   = Traegerfrequenz in kHz
%              phi0 = Anfangsphase
%
% Ausgabe:     t    = Zeitvektor
%              xpm  = PM-moduliertes Bandpass-Signal
%
% Anmerkung:   Bei f0 = 0 wird die komplexe Einhuellende erzeugt
%
% K.D. Kammeyer    17-jan-97

function [t,xpm] = pm_mod(t,v,dphi,f0,phi0);

if f0==0
  xpm=exp(j*(dphi.*v+phi0));            % komplexe Einhuellende
else
  xpm=cos(2*pi*f0.*t+dphi.*v+phi0);     % Bandpasssignal
end;

% ### EOF ######################################################################
