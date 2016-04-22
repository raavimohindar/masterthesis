% ##############################################################################
% ##  cosroll.m : Erzeugung eines Signals mit Cosinus-Roll-Off                ##
% ##              Charakteristik                                              ##
% ##############################################################################
%
% Aufruf:    [t,g]=cosroll(r,w,L);
%
% Eingabe:   r = Roll-Off-Faktor
%                bei r=0 wird eine bei +/- TL/2 abgebrochene si-Funktion
%                eingesetzt
%            w = Anzahl der Abtastintervalle pro Symbol ( T/TA ); w = gerade
%            L = Anzahl der Symboltakte T; Bedg.: L >= 1/r; L = gerade
%
% Ausgabe:   t = Zeitvektor auf T normiert: t --> t/T = -L/2 ... +L/2
%            g = cos-roll-off-Impuls
%
%                                                   Kammeyer '90

function [t,g]=cosroll (r,w,L);

LL=w*L;
t=[-LL/2:LL/2]/w;

if r == 0                                   % si-Funktion
  warning off;                              % Fehlermeldung unterdruecken
  g=sin(pi*t)./(pi*t);
  warning on;
  g(find(isnan(g)))=1;                      % da 0/0 bei sin-Term -> Grenzwert 1
  g(LL/2+1)=1;

else                                        % beliebige Roll-Off-Faktor
  warning off;                              % Fehlermeldung unterdruecken
  g=sin(pi*t)./(pi*t).*cos(r*pi*t)./(1-(2*r*t).^2);
  warning on;
  g(LL/2+1)=1;                              % da 0/0 bei sin-Term -> Grenzwert 1
  k0=w/(2*r);                               % 0/0 bei cos-Term
  if k0==fix(k0)                            % -> Grenzwert r/2*sin(pi/(2*r))
    g(LL/2+1+w/(2*r))=r/2*sin(pi/(2*r));
    g(LL/2+1-w/(2*r))=r/2*sin(pi/(2*r));
  end
end

% ### EOF ######################################################################
