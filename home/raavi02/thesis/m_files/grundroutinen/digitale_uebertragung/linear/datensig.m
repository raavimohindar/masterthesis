% ##############################################################################
% ##  datensig.m : Erzeugung eines Datensignals mit Impulsformung             ##
% ##############################################################################
%
% Aufruf:    [x,d_ref] = datensig(g,w,d,A);
%
% Eingabe:   g = Sendeimpuls (fA=W/T); z.B. cosroll oder wurzcos
%            w = Abtastwerte pro Symbolintervall (T/TA)
%            d = eingegebener Datenvektor (auch komplex zulaessig)
%            A = 1 -> Ein- und Ausschwinger abschneiden
%                     dann sind die ersten und die letzten L/2 Werte
%                     des eingegebenen Datenvektors irrelevant: z.B. L=4
%                     d     = [d(*) d(*) d(1) d(2) ... d(N) d(*) d(*)]
%                     d_ref =           [d(1) d(2) ... d(N)]
%                     x(1:w:length(x))= [d_dach(1) d_dach(2) .... d_dach(N)]
%
% Ausgabe:   x     = hochabgetastetes Datensignal
%            d_ref = Datenvektor zum Vergleich (eingegebener Vektor mit
%                    (bei A=1) abgeschnittenen Vor-u. Nachlaeufern)
%
% Anmerkung: Der Vergleich von d_ref und sign(x(1:w:length(x)) fuehrt direkt
%            zur Messung der BER

function [x,d_ref]=datensig(g,w,d,A)

N=length(d);
NN=(N-1)*w+1;
dd=zeros(1,NN);
dd(1:w:NN)=d(1:N);

x=conv(dd,g);
d_ref=d;

% Einschwingvorgang abschneiden
if A == 1
  x(1:length(g)-1)=[];
  x(length(x)-length(g)+w+1:length(x))=[];
  L=(length(g)-1)/w;
  d_ref(1:round(L/2))=[];
  d_ref(length(d_ref)-fix(L/2)+1:length(d_ref))=[];
end;

% ### EOF ######################################################################
