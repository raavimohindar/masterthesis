% ##############################################################################
% ##  awgn_sim.m : Simulation einer AWGN-Uebertragung im komplexen            ##
% ##               TP-Bereich mit vorgegebenen Sende- und Empfangsfiltern     ##
% ##############################################################################
%
% Aufruf:      [y,d_ref] = awgn_sim(g_S,g_E,w,d,Es_N0,kappa,A);
%
% Eingabe:     g_S   = Sendeimpuls (fA=W/T); z.B. wurzcos
%              g_E   = Empfangsfilter (fA=W/T); z.B. wurzcos
%              w     = Abtastwerte pro Symbolintervall (T/TA);  w:gerade
%              d     = eingegebener Datenvektor (0/1 oder 1/-1)
%              Es_N0 = Es/N0 in dB
%              kappa = Versatz bei der Symbolabtatsung (i.a. kappa =0)
%              A     = 1 -> Ein- und Ausschwinger abschneiden
%                           dann sind die ersten und die letzten L/2 Werte
%                           des eingegebenen Datenvektors irrelevant: z.B. L=4
%                           d     = [d(*) d(*) d(1) d(2) ... d(N) d(*) d(*)]
%                           d_ref =           [d(1) d(2) ... d(N)]
%                           x(1:w:length(x))= [d_dach(1) d_dach(2)... d_dach(N)]
%
% Ausgabe:     y     = Empfangssignal im Symboltakt
%              d_ref = Datenvektor zum Vergleich (eingegebener Vektor mit
%                      (bei A=1) abgeschnittenen Vor-u. Nachlaeufern)
%
% Anmerkung:   Der Vergleich von d_ref und sign(x(1:w:length(x)) fuehrt direkt
%              zur Messung der BER


function [y,d_ref]=awgn_sim(g_S,g_E,w,d,Es_N0,kappa,A)
N=length(d);
%d2_m=0.5*(max(d)+min(d));       % quad. MW. bei zweistfg. Signal
d2_m=sum(abs(d).^2)/length(d);

NN=(N-1)*w+1;
d_ref=d;
dd=zeros(1,NN);
dd(1:w:NN)=d(1:N);
x=conv(dd,g_S);

% Rauschueberlagerung
if Es_N0 >=100
  sig2_r=0;
else
  Es_N0=10^(Es_N0/10);     % Entlogarithmieren
  sig2_r= d2_m/2*sum(g_S.^2)/Es_N0;
end
sig_r=sqrt(sig2_r);
r=sig_r*randn(1,length(x))+j*sig_r*randn(1,length(x));

x=x+r;

y=conv(x,g_E);        % Empfangsfilter

% Einschwingvorgang abschneiden

if A == 1
  Lg=length(g_S)+length(g_E)-1;
  y(1:Lg-1)=[];
  y(length(y)-Lg+w+1:length(y))=[];
  L=(Lg-1)/w;
  d_ref(1:round(L/2))=[];
  d_ref(length(d_ref)-fix(L/2)+1:length(d_ref))=[];
end

% Symbolabtastung

y=y(1+kappa:w:length(y));  % kappa = 0 fuer L=gersde, sonst  kappa = w/2

% ### EOF ######################################################################
