% ##############################################################################
% ## integr.m : Integrierer                                                   ##
% ##############################################################################
%
% Aufruf:      [t,xint] = integr(t,x);
%
% Eingabe:     t = Zeitvektor
%              x = zu integrierender Vektor
%
% Ausgabe:     t    = Zeitvektor
%              xint = integrierter Vektor

function [t,xint] = integr(t,x)

delta_t = t(2)-t(1);
x_tmp = x(:).';        % Zeilenvektor erzwingen
Lx = length(x_tmp);

xx = reshape([[0 x_tmp(2:Lx)]; x_tmp], 1, 2*Lx);
xint = delta_t * cumsum(xx)/2;
xint = xint(1:2:length(xint));

if size(x,2)==1
  xint = xint.';
end

% ### EOF ######################################################################
