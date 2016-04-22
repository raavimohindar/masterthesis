% ##############################################################################
% ##  i_f_trafo.m : Inverse Fouriertransformation                             ##
% ##############################################################################
%
% Aufruf:     [t,x] = i_f_trafo(f,X)
%
% Eingabe:    f = [f_min:delta_f:f_max]       = Frequenz-Stuetzstellen
%             X = [X(f_min, ... , X(f_max)]  = Stuetzwerte des Spektrums
%
% Ausgabe:    t = [t_min:delta_t:t_max]      = zeitliche Stuetzstellen
%             x = [x(t_min), ... , x(t_max)] = Stuetzwerte des Zeitsignals
%
% Anmerkung:  Es koennen kausale und nichtkausale Signale erzeugt werden.
%
% Normierung: Wird die Frequenz in kHz eingegeben,
%             so ergibt sich die Zeiteinheit in ms.

function [t,x]=i_f_trafo(f,X);

N=length(X);
fA=2*max(f);
X=ifftshift(X);
x=fA*ifft(X,N);

% Rundungsfehler abschneiden
r=real(x)*10^10;
r=fix(r)/(10^10);
im=imag(x)*10^10;
im=fix(im)/(10^10);
x=r+j*im;
x=fftshift(x);

% Zeitvektor
t=[-N/2:N/2-1]/fA;

% ### EOF ######################################################################
