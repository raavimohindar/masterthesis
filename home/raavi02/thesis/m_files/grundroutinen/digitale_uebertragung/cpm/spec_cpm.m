% ##############################################################################
% ##  spec_cpm.m : berechtnet Leistungsdichtespektrum S(f*Tb)                 ##
% ##               von CPM Signalen                                           ##
% ##############################################################################
%
% Aufruf:    Spc = spec_cpm(R, struct_cpm_dat, fTb);
%
% Beispiele: S = spec_cpm(R, struct_cpm_dat, [-2.5 0.01 2.5]);
%
% Eingabe:   R     = Vektor der Autokorrelation, erzeugt mit akf_cpm
%            struct_cpm_dat = Daten aus akf_cpm
%            [fTb] = [start step stop] Festlegung des Frequenzausschnitts
%                    default: [-2 0.01 2]
%
% Ausgabe:   t = Zeitvektor
%            xH = hilberttransformiertes Signal
%
% Anmerkung: Die Funktion akf_cpm muss vor spec_cpm augerfen werden
%
% Weitere Funktionen, die aufgerufen werden: fourint
%
% spec_cpm berechnet das Leistungsdichtespektrum (PSD, power spectrum density
% fuer CPM Signale (continous phase modulation (FM)). Eine analytische
% Loesung existiert nur bedingt. Der hier implementierte Algorithmus
% stuetzt sich auf die Ableitung des PSD in [Sund 86], p.147 ff. .
% Die Berechnung des PSD wird dabei auf die Fouriertransformation einer
% zeitbegrenzten Autokorrelation (AKF) zurueckgefuehrt.
% SPEKCPM fuehrt im wesentlichen die Fouriertransformation der AKF durch.
% Diese Integration erfolgt mit einer modifizierten Simpsonformel. Dabei
% wird nur der Integrand quadratisch interpoliert und die Exponential-
% schwingung bezueglich der Integration exakt behandelt. Dies hat zur
% Folge, dass die Stuetzstellenanzahl bei wachsender Frequenz nicht
% mit wachsen muss, um etwa den Integranden angemmessen zu interpolieren.
% Der numerische Aufwand wird damit erheblich begrenzt.
%
% Hinweise zur Benutzung:
% Die optionalen Argumente beim Funktionsaufruf koennen nur "von hinten"
% weggelassen werden, weil die Funktion eine Zurordnung gemaess der Reihen-
% folge der uebergebenen Argumente zu den lokalen Argumenten innerhalb der
% Funktion vornimmt.
%
% Literatur:
% [Sund 86]  J.B. Anderson, T. Aulin, C.E. Sundberg.: Digital Phase
%            Modulation, Plenum Press, New York 1986.
%
%                     Benthin 1994
%                     Kammeyer 2000

function Spc = spec_cpm(R, struct_cpm_dat, fTb);

%fTb_def = [ 0 0.01 2.5 ];
fTb_def = [ -2.0 0.01 2.0 ]; % Kammeyer geaendert: Symmetrischer Frequenzbereich

if nargin == 2,           % Behandlung der optionalen Argumente
  fTb = fTb_def;
  %elseif nargin == 2,
end;


%load cpm_dat;            % cpm_dat wird durch akf_cpm erzeugt
% Variablen: L impulsform eta M w P f3dBT
L = struct_cpm_dat.L;
impulsform = struct_cpm_dat.impulsform;
eta = struct_cpm_dat.eta;
M = struct_cpm_dat.M;
w = struct_cpm_dat.w;
P = struct_cpm_dat.P;
f3dBT = struct_cpm_dat.f3dBT;

To = 1/w;
% Maximalanzahl der Elemente einer Matrix
comtyp = 6000;       % werden in comtyp abgelegt.

M_1 = M-1;
kk=1:2:M_1;
Ppos = P(1:2:M_1);
Pneg = P(2:2:M);

x = exp( j*pi*eta * kk' );
Calpha = Ppos*x + Pneg*conj(x);
jpi2ln = -j*2*pi*log(M)/log(2);
Spc = [];
fTb = fTb(1):fTb(2):fTb(3);

mumax = fix(comtyp/((L+1)*w));
lf = length(fTb);
II = fix(lf/mumax);
if rem(lf,mumax) == 0,
  II = II-1;
end
fB =[];
for I=0:II,          % Aufteilung des Frequenzbereiches in Teilabschnitte,
  % um mit nicht zu grossen Matrizen arbeiten zu koennen.
  lm = I*mumax +1;
  rm = (I+1)*mumax;
  fT = fTb(lm:min([lf rm]));
  fB = [fB,fT];
  % Hier werden mit dem Frequenzausschnitt fT die noetigen
  % Matrixoperationen  durchgefuehrt.

  S0L = fourint(0,L*w,R,fT,M,w,To);
  SL1 = fourint(L*w,((L+1)*w),R,fT,M,w,To);

  one = ones(1,length(fT));
  kal = one./(one - Calpha*exp(jpi2ln*fT));


  Spc =  [Spc abs( 2*real(S0L + kal.*SL1))];
end;
%plot(fTb-fB);              % Zu Testzwecken einblendbar. Es wird geprueft
%pause(2);                  % ob insgesamt genau fTb erfasst wird.
plot(fTb,10*log10(Spc));
grid;

etatext = sprintf('%.2f',eta);
if f3dBT ~= 0,
  f3dBTtext = sprintf('%.3f',f3dBT);
  f3dBTtext = [', f3dB*T = ' f3dBTtext];
else
  f3dBTtext = '';
end;

text = [int2str(L),' ',impulsform,f3dBTtext,',  eta= ',etatext,  ...
        ',  M= ',int2str(M),'                f * Tb  \rightarrow'];
title('Leistungsdichetespektrum    CPM-Signale')
xlabel(text);
ylabel('S_{SS}(j\omega)/T   in dB  \rightarrow');
pause(2);

% ### EOF ######################################################################
