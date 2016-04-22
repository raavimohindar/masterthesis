% ##############################################################################
% ##  akf_cpm.m : AKF der komplexen Einhuellenden von CPM-Signalen            ##
% ##############################################################################
%
% akf_cpm bestimmt die Autokorrelation R(t/T) der komplexen Einhuellenden
% von CPM Signalen als Vektor.
%
% Aufruf:    [r, struct_cpm_dat] = akf_cpm(L,impulsform,eta,[M],[w],[mmax],[P])
%
% Beispiele: [r, struct_cpm_dat] = akf_cpm(1,'RC',0.5,4,64);
%            [r, struct_cpm_dat] = akf_cpm([5 0.3],'GAUSS',0.5);
%
% Eingabe:   L          = Impulslaenge L*T     (L = [mu f3dBT] bei 'GAUSS')
%                          mu    = Gaussimpuls im +/- "mu*sigma" Intervall.
%                          f3dBT = 3dB-Frequenz des Gaussfilters * T
%            impulsform = 'REC','RC' oder 'GAUSS'
%            eta        = Modulationsindex
%            [M]        = Stufigkeit des gesendeten Datensignals; default: M=2
%            [w]        = Abtastwerte der AKF im Symbolintervall T;
%                         default: w=32
%            [mmax]     = Gew. Laenge der AKF in Einheiten von T;
%                         default: mmax=L
%            [P]        = Vektor der a-priori Wahrsch.; default: [1/M 1/M ...]
%                         z.B. zweistufig:  [P(1) p(-1)]
%
% Ausgabe:   r(0)...r(mmax)
%            struct_cpm_dat = Uebergabeparameter fuer spec_cpm
%
% [.] sind optionale Angaben und werden mit default Einstellungen belegt, wenn
% sie nicht explizit gesetzt werden. Argumente "von hinten" weglassen !!
%
% Ausfuehrliche Erlaeuterungen finden sich im Quelltext des M-Files.
% Weitere Funktionen, die aufgerufen werden: qpuls
%
%                        Benthin 1994
%                        Kammeyer 2000

function [R, struct_cpm_dat] = akf_cpm(L,impulsform,eta,M,w,mmax,P);

% ..........................................................................
%
%         SPEKTRALEIGENSCHAFTEN VON CONTINUOUS PHASE SYSTEMEN
%               Symbolische Darst. des CPM Systems
%
%    d(k)   -->     g(t)    -->   q(t) = integral[g(t)]
%   Daten       Impulsformer       "FM-Modulator"
%
%                        j * pi * eta * sum[d(i)*q(t-i*T)]
%              s(t)  =  e
%                      komplexe Einhuellende
%
% AKFCPM berechnet die mittlere Autokorrelation R(t/T) der komplexen Ein-
% huellenden s(t) eines CPM Signals (vgl. [Sund 86], p.150 ff.). Die AKF
% wird als Vektor ausgegeben. Die Funktion SPEKCPM benoetigt diese AKF
% als Argument.
%
% Erlaeuterungen zu den Argumenten:
%
% T:
% Das Zeitintervall T kennzeichnet den Symboltakt, mit dem die Daten d(k)
% gesendet werden. L ist eine ganze Zahl und L*T ist die Laenge der Impuls-
% antwort des Impulsformers g(t).
%
% impulsform:
% Es sind REC, RC und GAUSS Impulse fuer den Impulsformer g(t) moeglich.
% REC = Rechteckimpuls, d.h. harte Umtastung der Sendefrequenz ,
% RC = Raised Cosine Impuls ,  GAUSS = Gausssimpuls
% Fuer Gaussimpulse ist es nicht sinnvoll die Impulslaenge L unabhaengig
% von der 3dB Bandbreite "f3dBT" des Impulsformers festzulegen. Der Gauss-
% impuls entsteht aus der Faltung eines Rechteckes mit einer Gaussverteilung.
% Enstprechend ist die Argumentation ueber die Standardabweichung, die hier
% die zeitliche Ausdehnung beschreibt anschaulich. Die Faltung soll so
% durchgefuehrt werden, dass ein +/- mu * sigma Intervall der Gaussverteilung
% voll erfasst wird, (das Rechteck soll zu diesen Zeitpunkten nicht in die
% Gaussverteilung hineinragen). Die "Standardabweichung" sigma soll dabei
% eine Zeitangabe in t/T sein.
% Aus diesen Gruenden wird bei der Wahl 'GAUSS' nicht eine Laenge L des
% Impulsformers g(t) angegeben sondern ein Vektor [mu f3dBT], der die
% beschriebenen Groessen enthaelt. Das Programm AKFCPM bestimmt automatisch
% die naechste ganze Zahl L entsprechend der urspruenglichen Def. von L
% als Impulslaenge.
%
% eta:
% Der Modulationsindex "eta" ist gemaess eta = 2 * df * T definiert, wenn
% df den Frequenzhub bezeichnet.
%
% M:
% M bezeichnet die Stufigkeit der gesendeten Datenfolge. Ein Datum d(k)
% entstammt der Menge { -(M-1), .. , 3-, -1, 1, 3, .. (M+1) }, wenn
% M die Stufigkeit ist. M ist stets als gerade Zahl anzusetzen.
%
% w:
% Die Autokorrelation R(t) wird nur an diskreten Stellen t bestimmt,
% und zwar im Abtasttakt To. Es gilt T = w * To. w ist eine ganze
% Zahl, die durch zwei teilbar sein muss. Dies haengt damit zusammen,
% dass zur Integration spaeter die Simpsonformel verwendet wird und
% hierbei eine gerade Stuetzstellenanzahl benoetigt wird.
%
% mmax: Die Autokorrelation wird ausgewertet in einem Zeitintervall
% [0, (mmax+1)*T]. Fuer die Berechnung des PSD genuegt, die Festlegung
% mmax = L, die auch defaultmaessig verwendet wird, wenn nicht eine
% explizite Angabe fuer mmax erfolgt.
%
% P:
% Vektor der a-priori Wahrscheinlichkeiten mit denen ein bestimmtes
% Datum gesendet wird. Zur Illustration ein Beispiel (M = 8) :
%      P = [ P(1) P(2) P(3) P(4) P(5) P(6) P(7) P(8) ];
%              1   -1    2   -2    3   -3    4   -4
% In der letzten Zeile sind die zugehoerigen Daten angegeben. P(4)
% beschreibt in diesem Beispiel die a-priori Wahrscheinlichkeit, dass
% das gesendete Datum eine -2 ist. Die default Einstellung fuer P ist
% die Annahme von gleichverteilten Daten.
%
%
% Hinweise zur Benutzung:
% Die optionalen Argumente beim Funktionsaufruf koennen nur "von hinten"
% weggelassen werden, weil die Funktion eine Zurordnung gemaess der Reihen-
% folge der uebergebenen Argumente zu den lokalen Argumenten innerhalb der
% Funktion vornimmt.
%
%
% Literatur:
% [Sund 86]  J.B. Anderson, T. Aulin, C.E. Sundberg.: Digital Phase
%            Modulation, Plenum Press, New York 1986.
%
% Angabe der default Werte fuer die optionale Argumente:

w_def = 32;
M_def = 2;

% ........................................................................

if strcmp(impulsform,'GAUSS') == 1,      % Behandlung des Arguments L ab-
  mu = L(1);                             % haengig  davon ,ob ein GAUSS-
  f3dBT = L(2);                          % Impuls vorliegt oder nicht.
  L = ceil(1 + mu * sqrt(log(2))/(pi*f3dBT) );
else
  f3dBT = 0;
end;

if nargin == 3,                  % Behandlung der optionalen Argumente
  M = M_def;
  w = w_def;
  mmax = L;
  P = 1/M * ones(1,M);
elseif nargin == 4,
  w = w_def;
  mmax = L;
  P = 1/M * ones(1,M);
elseif nargin == 5,
  mmax = L;
  P = 1/M * ones(1,M);
elseif nargin == 6,
  P = 1/M * ones(1,M);
end;

struct_cpm_dat = struct('L', L, 'impulsform', impulsform, 'eta', eta,...
                        'M', M, 'w', w, 'P', P, 'f3dBT', f3dBT);

%save cpm_dat L impulsform eta M w P f3dBT;
% Diese Variablen werden in cpm_dat gespeichert, um
% sie in einer anderen Matlab Funktion einfach wieder
% zur Verfuegung stellen zu koennen, ohne sie als
% Argumente uebergeben zu muessen.

qq = qpuls(L,impulsform,w,0,f3dBT);    % Erzeugung des Grundimpulses q(t)
pause(2);                              % im  Zeitintervall [0,L*T]

q = [zeros(1,((mmax+1)*w)),qq,ones(1,(1+mmax)*w),1];
% Achtung: die letzte "1" von mir!!!
% Verlaengerte und verschobene Darstellung von q(t)

sh = (mmax+1)*w+1;     % Verschiebungskonstante
% Man findet den interessierenden Wert q(t) nun in
% q(t+sh), weil MATLAB keine negativen Zeiten bzw.
% Vektorindizes zulaesst.

R = zeros(0,((mmax+1)*w));     % Mittlere Autokorrelationsfunktion
abw = zeros(0,((mmax+1)*w));   % Schaetzungen der Abweichungen von der
% berechneten Autokorrelation zur wahren
% Autokorrelation zu jedem Zeitpunkt tau.
M_1 = M-1;
kk=1:2:M_1;
Ppos = P(1:2:M_1);       % Aufteilung des Vektors der a-priori Wahrschein-
Pneg = P(2:2:M);         % lichkeiten fuer pos. bzw. neg. Sendesymbole.
jpieta = j*pi*eta;


for m=0:mmax,                         % Berechnung der Autokorrelation
  for k_rel=0:w,
    yy = ones(1,w+1);
    for l=1-L:m+1,
      lm1 = -l*w +sh;
      lm2 = lm1 + k_rel + m*w;

      dq = q(lm2:lm2+w) - q(lm1:lm1+w);
      xx = zeros(1,w+1);

      X = exp( jpieta * kk'*dq );        % Summe ueber k, bezueglich der
      xx = Ppos*X + Pneg*conj(X);        % Wahrscheinlichkeiten P

      yy = yy .* xx;
    end; % of l
    % Integral ueber ein Zeitintervall T zur Mit-
    % telung der Autokorrelation zur Unterdrueckung
    % der Tatsache, dass die AKF zyklostationaer ist.
    Sn = (yy(1) + 4*sum(yy(2:2:w)) + 2*sum(yy(3:2:(w-1)))  ...
          + yy(w+1)) / (3*w);
    R(m*w+k_rel+1) = Sn;
    zz = [yy(1),yy(3:2:w+1)];

    % Berechnung des Integrals mit doppelter Schritt-
    % weite, zur Fehlerabschaetzung
    %    wh = w/2;
    %    Sm = (zz(1) + 4*sum(zz(2:2:wh)) + 2*sum(zz(3:2:(wh-1)))  ...
    %         + zz(wh+1)) / (3*wh);
    %    abw(m*w+k_rel+1) = abs(Sn-Sm)/17;

  end; % of k_rel
end; % of m

% Darstellung der Fehlerabschaetzung der Integralaus-
% wertung zur Bestimmung der mittleren Autokorrelation
%plot(abw);
%title('Fehlerabschaetzung der mittl. Autokorrelation');
%pause(5);

LR=length(R)-1;
tau=[-LR:LR]/w;
R_plot=[conj(R((LR+1):-1:2)) R];
plot(tau,real(R_plot),tau,imag(R_plot),'--');
%   konjugiert gerade Ergaenzung der AKF
%   zur  Plot-Ausgabe
%   Programmausgabe:  R(0)....R(taumax)
grid;
title('- Realteil   -- Imaginaerteil');
ylabel('r_{SS}(\tau)  \rightarrow');
xlabel('\tau/T  \rightarrow');
pause(2);

% ### EOF ######################################################################
