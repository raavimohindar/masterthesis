% ##############################################################################
% ##  traegerreg.m : Entscheidungsrueckgekoppelte Traegerphasenregelung fuer  ##
% ##                 QPSK                                                     ##
% ##############################################################################
%
% Aufruf:    [y_TP,ddach] = traegereg(y_channel,a,ref);
%
% Eingabe:   y_channel = Zu korrigierendes Signal
%                        gestoertes QPSK
%                        --> Erzeugung z.B. mit "phas_stoer"
%            a         = Koeffizient(en) der Regelung
%                        a=a0       --> erster Ordnung
%                        a=[a1 a2]  --> zweiter Ordnung
%            ref       = QPSK-Referenzsignal (in pi/4-Lage)
%                        Laenge = Lref: kleiner oder gleich Laenge y_TP
%                        bis Lref=Traininngs-Praeambel
%                        danach entschiedene Daten zur Phasenschaetzung
%
% Ausgabe:   y_TP      = Korrigiertes Signal
%            ddach     = entschiedene Symbole
%                        die ersten Lref Werte == ref
%                        danach = entschiedene Werte
%
%                                                   Kammeyer

function [y_TP,ddach] = traegerreg(y_channel,a,ref);

y_TP=y_channel;
detect = 0;
ddach = zeros(length(y_TP),1);
lref = length(ref);
phi = 0;         % Initialisierung mit den beiden Anfangsphasen
phi_old = 0;
DO = 0;

if (length(a) == 1) % Trgaegerphasenregelung 1. Ordnung
  % Fuer die Referenzdaten
  for i = 1:lref
    y_TP(i)= y_TP(i) .* exp(-j*phi); % Korrektur des Signals
    % Messung der Phasen Differenz und Multiplikation mit
    % Filterkoeffizient
    ddach(i) = (sign(real(y_TP(i))) + j*sign(imag(y_TP(i))))/sqrt(2);
    DO= a*imag(y_TP(i) * conj(ref(i)));
    % Integrator
    phi = DO + phi;
  end;
  % Fuer alle anderen Daten
  for i = lref+1:length(y_TP)
    y_TP(i)= y_TP(i) .* exp(-j*phi); % Korrektur des Signals
    % Entscheidung des Symbols + konj. Komplex (QPSK)
    detect = (sign(real(y_TP(i))) - j*sign(imag(y_TP(i))))/sqrt(2);
    ddach(i) = conj(detect);
    % Messung der Phasen Differenz und Multiplikation mit Filterkoeffizient
    DO= a*imag(y_TP(i) * detect); % diff_phase(k-1)
    % Integrator
    phi = DO + phi;
  end;
  % Ausgabe der letzten Phase
  corr = phi;

elseif (length (a) == 2) % Trgaegerphasenregelung 2. Ordnung
  for i = 1:lref
    y_TP(i)= y_TP(i) .* exp(-j*phi); % Korrektur des Signals
    % Entscheidung des Symbols + konj. Komplex (QPSK)
    detect = conj(ref(i));
    ddach(i) = (sign(real(y_TP(i))) + j*sign(imag(y_TP(i))))/sqrt(2);
    DO_old = DO; % diff_phase(k-2)
    DO = imag(y_TP(i) * detect); % diff_phase(k-1)
    phi_neu = 2*phi - phi_old + a(1)*DO + a(2)*DO_old; % Phase berechnen
    phi_old = phi; % Phasen umschreiben
    phi = phi_neu;
  end;

  for i = lref+1:length(y_TP)  % Schleife ueber den Signalblock
    y_TP(i)= y_TP(i) .* exp(-j*phi); % Korrektur des Signals
    % Entscheidung des Symbols + konj. Komplex (QPSK)
    detect = (sign(real(y_TP(i))) - j*sign(imag(y_TP(i))))/sqrt(2);
    ddach(i) = conj(detect);
    DO_old = DO ; % diff_phase(k-2)
    DO = imag(y_TP(i) * detect) ; % diff_phase(k-1)
    phi_neu = 2*phi - phi_old + a(1)*DO + a(2)*DO_old; % Phase berechnen
    phi_old = phi; % Phasen umschreiben
    phi = phi_neu;
  end;
else                                    % Abfangen von Fehlern
  disp('Fehler bei der Parametereingabe') ;
end;

% ### EOF ######################################################################
