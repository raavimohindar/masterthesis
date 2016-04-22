% ##############################################################################
% ##  f_trafo.m : Fouriertransformation                                       ##
% ##############################################################################
%
% Aufruf:     [f,X] = f_trafo(t,x,Nt)
%
% Eingabe:    t = [t_min:delta_t:t_max]      = zeitliche Stuetzstellen
%             x = [x(t_min), ... , x(t_max)] = Abtastwerte des Zeitsignals
%             Nt = Zeitintervall, ueber das integriert wird
%             Nt=1 : Lange des eingegebenen Zeitsignals
%                    zweckmaessig, falls x(t) periodisch fortgesetzt
%             Nt=2,3,4 ... Zeitsignal auf Nt-fache Laenge mit Nullen aufgefuellt
%                          je groesser Nt, desto feiner wird Frequenzaufloesung
%
% Ausgabe:    f = [f_min:delta_f:f_max]       = Frequenz-Stuetzstellen
%             X = [X(f_min, ... , X(f_max)]  = Stuetzwerte des Spektrums
%
% Anmerkung:  Wegen der Uebergabe des Zeitvektors t koennen
%             kausale und nichtkausale Signale transformiert werden.
%
% Normierung: Wird die Zeiteinheit in ms eingegeben,
%             so ergibt sich die Frequenz in kHz.
%
%
% D.Boss, K.D.Kammeyer, D.Nikolai   24-jan-97

function [f,X] = f_trafo(t,x,Nt)

N = Nt*length(t);
X = fft(x(:),N);
DeltaT = t(2)-t(1);
t0 = t(1);

% Phasenkorrektur
phi = [0:length(X)-1]*2*pi*t0/(N*DeltaT);
X = DeltaT*X.*exp(-j*phi(:));

% Rundungsfehler abschneiden
R = real(X)*10^10;
R = fix(R)/(10^10);
I = imag(X)*10^10;
I = fix(I)/(10^10);
X = R+j*I;
X = fftshift(X);

if size(x,2)~=1
  X = X.';
end;

% Frequenzvektor
if 2*round(N/2)==N
  f = [-N/2:N/2-1]/(N*DeltaT);
else
  f = [-N/2+0.5:N/2-0.5]/(N*DeltaT);
end

% ### EOF ######################################################################
