% ##############################################################################
% ##  trunc.m : Darstellung eines waehlbaren Zeitausschnitts                  ##
% ##############################################################################
%
% Aufruf:    [t_trunc,x_trunc]= trunc(t,x,t_u,t_o);
%
% Eingabe:   t = [t_min:delta_t:t_max]      = urspruenglicher Zeitvektor
%            x = [x(t_min), ... , x(t_max)] = urspruenglicher Signalvektor
%            t_u = untere Grenze des Zeitausschnitts
%            t_o = obere Grenze des Zeitausschnitts
%
% Ausgabe:   t_trunc = [t_u:delta_t:t_o]      = neuer Zeitvektor
%            x_trunc = [x(t_u), ... , x(t_o)] = neuer Signalvektor

function [t_trunc,x_trunc]=trunc(t,x,t_u,t_o);

L=length(x);
delta_t=t(2)-t(1);

deltamin=-min(t)+t_u;
deltamin=deltamin/delta_t;
deltamin=round(deltamin);

deltamax=max(t)-t_o;
deltamax=deltamax/delta_t;
deltamax=round(deltamax);

t_trunc=t;
t_trunc(L-deltamax:L)=[];
t_trunc(1:deltamin)=[];

x_trunc=x;
x_trunc(L-deltamax:L)=[];
x_trunc(1:deltamin)=[];

% ### EOF ######################################################################
