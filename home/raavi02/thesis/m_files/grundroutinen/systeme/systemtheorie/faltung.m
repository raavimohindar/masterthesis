% ##############################################################################
% ##  faltung.m : Faltung                                                     ##
% ##############################################################################
%
% Aufruf:      [ty,y] = faltung(tx1,x1,tx2,x2);
%
% Eingabe:     tx1,tx2 = [t_min:delta_t:t_max] = zeitliche Stuetzstellen
%
%              Stuetstellen-Abstaende delta_t muessen fuer
%              tx1,tx2 und ty identisch sein.
%
%              x1,x2 = [x(t_min), ... ,x(t_max)] = Stuetzwerte der Zeitsignale
%
% Ausgabe:     y(t)=x1(t)*x2(t)
%
% Anmerkung:   Wegen der Uebergabe des Zeitvektors t koennen
%           kausale und nichtkausale Signale gefaltet werden.
%
% Normierung:  Die Zeiteinheit ist in ms einzugeben.


function [ty,y]=faltung(tx1,x1,tx2,x2);

delta_t=tx1(2)-tx1(1);

if ((tx2(2)-tx2(1))-delta_t) > 1e-10
  error(['ERROR (faltung.m): tx1 und tx2 muessen gleichen ',...
         'Stuetzstellenabstand haben.']);
end;
t_min=tx1(1)+tx2(1);

% Faltung
y=conv(x1,x2)*delta_t;

% Zeitvektor
ty=[1:length(y)]*delta_t;
ty=ty+t_min-delta_t;

% ### EOF ######################################################################
