% ##############################################################################
% ##  multisin.m : Multisinussignal mit zufaelligen Phasen                    ##
% ##############################################################################
%
% Aufruf:      [t,x] = multisin(t,fi,ai);
%
% Eingabe:     t  = Zeitvektor
%              fi = Vektor der diskreten Frequenzen
%              ai = Vektor der Frequenzgewichte
%
% Ausgabe:     t = Zeitvektor
%              x = Multisinussignal
%
% Normieung:   Das Signal wird auf den |Maximalwert| eins normiert
%
% Beispiel:    fa = 100;
%              t  = [-2.5:1/fa:2.5-1/fa];
%              fi = [0.4:0.2:3.0];
%              agrenz = 0.25;
%              ai = [1:(agrenz-1)/(length(fi)-1):agrenz];
%              [t,v_multisin] = multisin(t,fi,ai);

function [t,x] = multisin(t,fi,ai)

if size(ai) == [1,1]
  ai = ai * ones(size(fi));
end

phi = rand(size(fi))*2*pi;
x   = zeros(size(t));

for cnt=1:length(fi)
  x = x + ai(cnt)*sin(2*pi*fi(cnt).*t + phi(cnt));
end

% Normierung
[maxwert, maxindex] = max(abs(x));
x = x/maxwert;

% ### EOF ######################################################################
