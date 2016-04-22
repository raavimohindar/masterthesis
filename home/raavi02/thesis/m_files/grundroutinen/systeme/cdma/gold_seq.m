% ##############################################################################
% ## gold_seq.m : Berechnung von {Gold}-{Seq}quenzen                          ##
% ##############################################################################
%
% Aufruf y = gold_seq(m,delay,fb_poly1,fb_poly2,...
%                       start1,start2[,n_chip[,oversmpl[,g_p]]])
%
% Eingabe:
%
%   m       :  Laenge des Schieberegisters 
%   delay   :  Verzoegerung des zweiten Muttercodes (innerhalb [0 ... 2^m-2])
%   fb_poly1:  Feedback Polynom von Code 1, stellt Feedback Pol. in GF(2) dar 
%   fb_poly2:  Feedback Polynom von Code 2, stellt Feedback Pol. in GF(2) dar 
%   start1  :  Start Vektor (Laenge m) von Mutter Code 1
%   start2  :  Start Vektor (Laenge m) von Mutter Code 2
%   n_chip  :  Anzahl der auzugebenen Chips (optional)
%   oversmpl:  Ueberabtastungsfaktor (optional)
%   g_p     :  Impulsform der Chips (optional)
%
% Ausgabe:
%
%   y       :  Matrix mit Goldsequenz in jeder Spalte
%              erste und zweite Spalte enthaelt Mutter-m-Sequenzen
%              Anzahl der Spalten entspricht der Laenge der Verzoegerung + 2
% Bemerkungen:
%
% Weitere Details siehe  Kammeyer "Nachrichtenuebertragung", p. 628 ff
%                        (weisses Buch)
%
% Author: Armin Dekorsy  19.06.96, Volker Kuehn, 06.07.01

function y = gold_seq(m,delay,fb_poly1,fb_poly2,start1,...
                      start2,n_chip,oversmpl,g_p)

m = length(fb_poly1)-1;

fb_poly1 = fb_poly1(:); 
fb_poly2 = fb_poly2(:);

start1 = start1(:);
start2 = start2(:);

if (nargin<7)
  n_chip = 2^m-1;     % Default: eine Periode
end

if (nargin<8)
  oversmpl = 1;
  g_p      = 1;
end

if (nargin==8)
  g_p = ones(oversmpl,1);   % Default: Rechteckimpulse
end

y = zeros(n_chip*oversmpl,length(delay));

y1 = sign(m_seq(fb_poly1,start1,n_chip));
y2 = sign(m_seq(fb_poly2,start2,n_chip));

for tau = 1:length(delay)
  tmp = y1 .* [y2(delay(tau)+1:n_chip);y2(1:delay(tau))];
  tmp = [tmp'; zeros(oversmpl-1,n_chip)];               % Ueberabtastung
  tmp = tmp(:)./sqrt(n_chip);                           % Normierung

  y(:,tau) = conv(tmp(1:oversmpl*(n_chip-1)+1),g_p);
end;

% Anfuegen der Muttersequenz zur Ausgangssequenz
y1 = [y1'; zeros(oversmpl-1,n_chip)];                    % Ueberabtastung
y1 = y1(:)./sqrt(n_chip);                                % Normierung
y1 = conv(y1(1:oversmpl*(n_chip-1)+1),g_p);

y2 = [y2'; zeros(oversmpl-1,n_chip)];                    % Ueberabtastung
y2 = y2(:)./sqrt(n_chip);                                % Normierung
y2 = conv(y2(1:oversmpl*(n_chip-1)+1),g_p);

y = [y1 y2 y];

% ### EOF ######################################################################
