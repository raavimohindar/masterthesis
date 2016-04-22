% ##############################################################################
% ##  lms.m : Adaptiver Entz. mit iterativer Einstellung nach dem             ##
% ##          LMS-Algorithmus                                                 ##
% ##############################################################################
%
% Aufruf:    [e_H,epsilon,y] = lms(x_channel,d_ref,n,L_Train,i_0,delta,e_H_in);
%
% Eingabe:   x_channel : Empfangsvektor
%            d_ref     : Vektor der Referenzsignale bestehend aus
%                        L_Train Trainingssymbolen und L_Data
%                        Informationssymbolen Gesamtlaenge L=L_Train+L_Data
%            n         : Anz. der Entzerrerkoeffizienten
%            L_Train   : Anz. der Trainingssymbole
%            i_0       : Position der Eins, Wertebreich 0,...,m+n
%            delta     : Schrittweise des LMS
%            e_H_in    : Eingangswert der Entzerrerkoeffizienten (optional)
%
% Ausgabe:   e_H       : Matrix der iterativ bestimmten Entzerrerkoeffizienten
%            epsilon   : Differenzsignal y(i)-d_ref(i-i_0) der Laenge L
%            y         : Ausgangssignal des Entzerrers der Laenge L
%
% Erlaeuterung:
%
% Mit der Entzerrung wird eine moeglichst gute Approximation von d_ref(i-i_0)
% mit y(i) angestrebt. Die Differenz zwischen gesendetem Symbol und
% approximaierten Symbol ergibt sich demanch zu epsilon(i) = y(i)-d_ref(i-i_0).
%
% Insgesamt liegen L = L_Train+L_Data Referenzsignale d_ref vor, von denen die
% ersten L_Train Symbole zur iterativen Einstellung der Entzerrerkoeffizienten
% verwendet werden:
% d_ref = [d_ref(1) d_ref(2) ... d_ref(i_0+1) d_ref(i_0+1) ... d_ref(L)]^T
%
% Demnach wird das erste Referenzsignale d_ref(1) zum Zeitpunkt i=i_0+1 mit
% dem Ausgangssignal y(i_0+1) des Entzerrers verglichen. Ebenso wird das letzte
% Referenzsignale d_ref(L) zum Zeitpunkt i=L+i_0 mit dem y(L+i_0) verglichen,
% wie die folgende Tabelle verdeutlicht.
%
% i         d_ref          y             epsilon
% ------------------------------------------------------------------------------
% 1         d_ref(1)       0             0
% 2         d_ref(2)       0             0
% :         :              :             :
% i_0       d_ref(i_0)     0             0
% i_0+1     d_ref(i_0+1)   y(i_0+1)      epsilon(i_0+1) = y(i_0+1) - d_ref(1)
% i_0+2     d_ref(i_0+2)   y(i_0+2)      epsilon(i_0+2) = y(i_0+2) - d_ref(2)
% :         :              :             :
% L         d_ref(L)       y(L)          epsilon(L)   = y(L)   - d_ref(L-i_0)
% L+1       0              y(L+1)        epsilon(L+1) = y(L+1) - d_ref(L-i_0+1)
% :         :              :             :
% L+i_0-1   0              y(L+i_0-1)    epsilon(L+i_0-1) = y(L-1) - d_ref(L-1)
% L+i_0     0              y(L+i_0)      epsilon(L+i_0)   = y(L)   - d_ref(L)
%
% Folglich werden nur die Ausgangssignale y(i_0+1) bis y(L+i_0) als
% Approximation der Daten d_ref(1) bis d_ref(L) und ebenso epsilon(i_0+1) bis
% epsilon(L+i_0) als Differenz der Approximation zu den gesendeten Daten von
% der Funktion lms.m ausgeben. Damit liegen der Testfunktion die Referenzdaten
%  d_ref und die approximierten Daten y "synchronisiert" vor.
% ==============================================================================

function [e_H,epsilon,y] = lms(x_channel,d_ref,n,L_Train,i_0,delta,e_H_in);

% Falls keine Entzerrerkoeffizienten uebergeben wurden, werden zur
% Initialisierung alle Koeffizienten zu Null gesetzt und nur einer zu Eins.
% Im ersten Zeitpunkt der Approximation i=i_0+1 sind die ersten i_0+1 Werte des
% Vektors x_i ungleich Null und die letzten n-i_0-1 Werten gleich Null. Damit
% y(i) = e_H*x_i einen Werte ungleich Null annimmt (verhindert Fehlermeldung bei
% der Decodierung), sollte die Eins des Entzerrers zwischen e_(0) und e_(i_0)
% liegen, wie folgende Tabelle fuer den Zeitpunkt i=i_0+1 verdeuticht.
%
% x_i               e_H
% ------------------------------------------------------------------------------
% x(i_0+1)          e_(0)
% x(i_0)            e_(1)
% x(i_0-1)          e_(2)
% :                 :
% x(1)              e_(i_0)       -> spaetestens hier sollte die Eins liegen
% x(0)              e_(i_0+1)     -> ab hier enthaelt x_i nur Werte gleich Null
% :                 :
% x(i_0+1-n)        e_(n)
%
% ==> Gilt i_0 >= [n+1/2]-1, so wird e_([n+1/2])=1 gesetzt, ansonsten  e_(i_0)=1

if nargin<7;
  e_H_in = zeros(1,n+1);
  if i_0 >= ceil((n+1)/2)-1
    e_H_in(ceil((n+1)/2)) = 1;
  else
    e_H_in(i_0+1) = 1;
  end
end;

if nargin<6; delta = 0.1;      end;      % Schrittweite des LMS
if nargin<5; i_0   = fix(n/2); end;      % Position der Eins der
%                                          Gesamtimpulsantwort

% Deklarieren der Variablen
e_H      = zeros(L_Train,n+1);          % Deklarieren der Entzerrerkoeffizienten
e_H(1,:) = e_H_in;                      % Entzerrerkoeffizienten initialisieren

epsilon  = zeros(length(d_ref)+i_0,1);   % Spaltenvektor der Differenzsignale
y        = zeros(length(d_ref)+i_0,1);   % Spaltenvektor der Filterausgabewerte
x_i      = zeros(n+1,1);                 % Spaltenvektor der betrachteten
%                                          Empfangswerte

x_channel = [x_channel; zeros(i_0,1)];   % Anhaengen von Nullen an den
%                                          Empfangsvektor
%                                          [x(1) x(2) ... x(L) 0 ... 0]^T
% Adaptive Entzerrung
for i=1:length(d_ref)+i_0
  x_i = [x_channel(i); x_i(1:length(x_i)-1)]; % Auswahl der Empfangswerte
  %                                             ("Fenster")
  % x_i = [x(i) x(i-1) ... x(i-n)]^T

  if i<=i_0                             % Ueberspringen der ersten i_0 Werte
    e_H(i+1,:) = e_H(i,:);              % Kopieren der Entzerrerkoeffizienten

  elseif i<=L_Train                     % Schleife durch Trainingssequenz
    y(i)       = e_H(i,:)*x_i;          % Ausgangswert des Enzerrers
    epsilon(i) = y(i)-d_ref(i-i_0);     % Bestimmen des Diff.-Signals
    e_H(i+1,:) = e_H(i,:)-delta*epsilon(i)*x_i'; % e_H(i+1) aus e_H(i) bestimmen

  else                                  % Schleife durch Informationsdaten
    y(i)       = e_H(L_Train,:)*x_i;    % Ausgangswert des Enzerrers
    epsilon(i) = y(i)-d_ref(i-i_0);     % Bestimmen des Diff.-Signals
  end
end

y(1:i_0)       = [];            % Entfernen der ersten i_0 Werte (vgl. Tabelle)
epsilon(1:i_0) = [];            % Entfernen der ersten i_0 Werte (vgl. Tabelle)

% ### EOF ######################################################################
