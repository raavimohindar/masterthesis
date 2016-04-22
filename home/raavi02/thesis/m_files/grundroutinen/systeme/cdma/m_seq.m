% ##############################################################################
% ## m_seq.m : Generierung von {m}-{Seq}uenzen                                ##
% ##############################################################################
%
% Aufruf  y = m_seq(fb_poly,start[,n_chip[,oversmpl[,g_p]]])
%
% Eingabe :
%
%   fb_poly :  Vektor mit Feedback Polynom in GF(2)
%   start   :  Startvektor der Laenge m (Register Initialisierung) 
%   n_chip  :  Anzahl der Ausgabechips (optional)
%   oversmpl:  Ueberabtastungsfaktor (optional)
%   g_p     :  Impulsform der Chips (optional, Vorsicht mit Impulsenergie!)
%
% Ausgabe:
%
%   y       :  Spaltenvektor (n-Elemente) mit m-Sequenzen {0,1} 
%
% Bemerkungen:
% Weitere Details siehe Kammeyer "Nachrichtenuebertragung", p. 625 ff
%
% Author: Armin Dekorsy  19.06.96, Volker Kuehn, 02.07.01

function [y] = m_seq(fb_poly,start,n_chip,oversmpl,g_p)

fb_poly(1) = [];       % loeschen der '1' fuer Feedback-Schleife
m          = length(fb_poly);

if (length(start)~=m)
  disp('m_seq.m: Length of fb_poly and start do not match!');
  break
end

if (nargin<3)
  n_chip = 2^m-1;     % Default: eine Periode
end

if (nargin<4)
  oversmpl = 1;
  g_p      = 1;
end

if (nargin==4)
  g_p = ones(oversmpl,1)/sqrt(oversmpl);   % Default: Rechteckimpulse
end

reg     = start(:);
fb_poly = find(fb_poly(:));

y       = zeros(n_chip,1);                 % Bereitstellung des Vektors y      
y(1)    = start(m);                        % erstes Element

for lauf = 2:n_chip                        % Chips entdecken
  fb   = rem(sum(reg(fb_poly)),2);
  reg  = [fb; reg(1:m-1)];
  y(lauf) = reg(m);         
end;

y = 1-2*y(:)';                             % BPSK Modulation
y = [y; zeros(oversmpl-1,n_chip)];         % Ueberabtastung
y = y(:)./sqrt(n_chip);                    % Normierung

y = conv(y(1:oversmpl*(n_chip-1)+1),g_p);

% ### EOF ######################################################################
