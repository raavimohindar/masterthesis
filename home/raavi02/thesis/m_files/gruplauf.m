% ##############################################################################
% ##  gruplauf.m : Berechnung der Gruppenlaufzeit eines FIR-Filters           ##
% ##############################################################################
%
% Aufruf:    [f,taug] = gruplauf(t,h,Nt);
%
% Eingabe:   t  = Zeitvektor
%            h  = Abtastwerte des Zeitsignals
%            Nt = Zeitintervall, ueber das integriert wird
%            - Nt=1 : Lange des eingegebenen Zeitsignals
%                     zweckmaessig, falls periodisch fortgesetzt
%            - Nt=2,3,4 ...  Zeitsignal auf Nt-fache Laenge mit Nullen 
%                            aufgefuellt
%                     je groesser Nt, desto feiner wird Frequenzaufloesung
%
% Ausgabe:   f    = Frequenzvektor [-f_A/2 ... +f_A/2]
%            taug = Gruppenlaufzeit

function [f,taug]=gruplauf(t,h,Nt)

N=Nt*length(t);
DeltaT=t(2)-t(1);
t0=t(1);

n=length(h);
for i=1:n
  hs(i)=h(i)*(i-1);
end

h(N)=0;                % Auffuellen mit Nullen
hs(N)=0;               % Auffuellen mit Nullen

H=fft(h);              % Uebertragungsfuntion des FIR-Filters
H2=imag(H);
H1=real(H);

Hs=fft(hs);            % abgeleitete Uebertragungsfunktion des FIR-Filters
Hs2=imag(Hs);
Hs1=real(Hs);

taug=(H1.*Hs1+H2.*Hs2)./(H1.^2+H2.^2);  % Herleitung KK98
taug=DeltaT*taug+t0;

if 2*round(N/2)==N     % Frequenzvektor
  f=[-N/2:N/2-1]/(N*DeltaT);
else
  f=[-N/2+0.5:N/2-0.5]/(N*DeltaT);
end

% ### EOF ######################################################################
