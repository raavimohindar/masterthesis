% ##############################################################################
% ##  ezsim.m :  Simulation eines Entzerrers mit Uebergabe der Kanal- u.      ##
% ##              Entzerrer-Impulsantwort (wahlweise w=1 und w=2)             ##
% ##############################################################################
%
% Aufruf:    [y,d_ref]=ez_sim(d,h,w,e,ph,EsN0,di0);
%
% Eingabe:   d          = Datenvektor (auf Leistung eins normiert!)
%            h          = Kanalimpulsantwort im w-fachen Takt
%            w          = 1 --> Symboltaktentzerrer
%                         2 --> T/2-Entzerrer
%            e          = Entzerrer-Impulsantwort (w=1 oder 2)
%            ph         = Abtastphase am Entzerrer-Ausgang: ph = 0/1
%            EsN0 in dB = Angabe fuer Rauschueberlagerung
%            di0        = Verschiebung der "1" aus der Mittellage
%
% Ausgabe:   y          = Enterrtes Signal im Symboltakt
%                        Einschwingvorgaenge werden abgeschnitten
%            d_ref      = Datenvektor mit abgschnittenen Ein- und Ausschwinger
%
% Anmerkung: Bei w=2 ist am Matchedfilterausgang das Rauschen i.a. nicht weiss;
%            dieser Effekt wird nicht beruecksichtigt

function [y,d_ref]=ez_sim(d,h,w,e,ph,EsN0,di0);

N=length(d);
dd=zeros(1,w*N);
dd(1:w:w*N)=d(1:N);
% Bestimmen des Kanalausgangssignals x
h=h(:);
h=h/sqrt(h'*h);         % Normierung auf die Energie 1
x=conv(dd,h);
x(1:length(h)-1)=[];
x(length(x)-length(h)+2:length(x))=[];
sig2_N=10^(-EsN0/10);   % Leistung des komplexen Kanalrauschens
n=sqrt(sig2_N/2)*(randn(1,length(x))+j*randn(1,length(x)));
% Ueberlagerung von Rauschen
x=x+n;

% Bestimmen des Entzerrerausganssignals y
y=conv(x,e);
y(1:length(e)-1)=[];
y(length(y)-length(e)+2:length(y))=[];
y1=y(1+ph:w:length(y));
y=y1;

% Synchronisieren der Referenzdaten mit den Entzerreraussganssignalen
d_ref=d;
Speicher = length(h)+length(e)-1;               % Speicheranzahl

if w==1                                         % w=1
  i0 = fix(Speicher/(2*w)) + di0 +1;            % Umrechnung von di0 -> io
  d_ref(1:Speicher-i0)=[];                      % Loeschen der Referenzdaten
else                                            % w=2
  i0 = round(Speicher/4)+di0+1;                 % Umrechnung von di0 -> io
  d_ref(1:round(Speicher/2)-i0+ph)=[];          % Loeschen der Referenzdaten
end

d_ref(length(y)+1:length(d_ref))=[];            % Loeschen der Aussschwinger

% ### EOF ######################################################################
