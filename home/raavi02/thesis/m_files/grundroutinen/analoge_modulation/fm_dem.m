% ##############################################################################
% ## fm_dem.m : {FM}-{Dem}odulator                                            ##
% ##############################################################################
%
% Aufruf:      [t,vtilde] = fm_dem(t,xFM,dF,f0,phi0);
%              benoetigt bbmisch.m
%
% Eingabe:     t    = Zeitvektor
%              xfm  = FM-moduliertes Bandpass-Signal
%              dF   = Frequenzhub in kHz (zur korrekten Skalierung)
%              f0   = Traegerfrequenz in kHz
%              phi0 = Anfangsphase
%
% Ausgabe:     t      = Zeitvektor
%              vtilde = demoduliertes Signal
%
% K.D. Kammeyer    17-jan-97

function [t,vtilde] = fm_dem(t,xFM,dF,f0,phi0);

[t,s,xplus] = bbmisch(t,xFM,f0);    % Basisbandmischung
s=s*exp(-j*phi0);                   % Anfangsphase

Ta=t(2)-t(1);
AM=(real(s)).^2+(imag(s)).^2;
s=s./sqrt(AM);
s0=s;
s1=[s(length(s)) s];
s1(length(s1))=[];
vtilde=imag(s0).*real(s1) - real(s0).*imag(s1);
vtilde=asin(vtilde)/(2*pi*dF*Ta);

% ### EOF ######################################################################
