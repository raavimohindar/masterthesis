% ##############################################################################
% ##  modquel.m : Datenquelle mit verschied. Modulationsformen                ##
% ##              (zufaellige Daten)                                          ##
% ##############################################################################
%
% Aufruf:    d=modquel(N,ART);
%
% Eingabe:   d   = Symbolfolge
%
% Ausgabe:   N   = Anzahl der Symbole
%            ART = Modulationsform
%              1 = antipodal, 2-stufig
%              2 = ASK, 4-stufig
%              3 = QPSK (lambda =0)
%              4 = QPSK (lambda = pi/4)
%              5 = 8PSK
%              6 = 16PSK
%              7 = 16QAM
%              8 = 16ASK/PSK (V29)

function d=modquel(N,ART);

if ART == 1
  d=sign(randn(1,N));

elseif ART == 2
  dvz=sign(randn(1,N));
  d0=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d=dvz.*(1+2*d0);
  %d=d/sqrt(5);     % Normierung der mittleren Leistung auf eins

elseif ART == 3
  d0=sign(randn(1,N));
  d1=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d1=0.5*(d1+ones(1,N));
  d=d0+2*d1;
  d=exp(j*2*pi .*d/4);

elseif ART == 4
  d0=sign(randn(1,N));
  d1=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d1=0.5*(d1+ones(1,N));
  d=d0+2*d1;
  d=exp(j*pi/4)*exp(j*2*pi .*d/4);

elseif ART == 5
  d0=sign(randn(1,N));
  d1=sign(randn(1,N));
  d2=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d1=0.5*(d1+ones(1,N));
  d2=0.5*(d2+ones(1,N));
  d=d0+2*d1 +4*d2;
  d=exp(j*2*pi .*d/8);

elseif ART == 6
  d0=sign(randn(1,N));
  d1=sign(randn(1,N));
  d2=sign(randn(1,N));
  d3=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d1=0.5*(d1+ones(1,N));
  d2=0.5*(d2+ones(1,N));
  d3=0.5*(d3+ones(1,N));
  d=d0+2*d1+4*d2+8*d3;
  d=exp(j*2*pi .*d/16);

elseif ART == 7
  z(1)= 1-j;
  z(2)= 1+j;
  z(3)=-1+j;
  z(4)=-1-j;
  z(5)= 1-3*j;
  z(6)= 3-3*j;
  z(7)= 3- j;
  z(8)= 3+ j;
  z(9)= 3+3*j;
  z(10)= 1+3*j;
  z(11)=-1+3*j;
  z(12)=-3+3*j;
  z(13)=-3+j;
  z(14)=-3-j;
  z(15)=-3-3*j;
  z(16)=-1-3*j;
  d0=sign(randn(1,N));
  d1=sign(randn(1,N));
  d2=sign(randn(1,N));
  d3=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d1=0.5*(d1+ones(1,N));
  d2=0.5*(d2+ones(1,N));
  d3=0.5*(d3+ones(1,N));
  d=d0+2*d1+4*d2+8*d3;
  d=d+ones(1,N);
  d=z(d);
  %d=d/sqrt(40/4);   % Normierung der Leistung auf eins

elseif ART == 8
  z(1)= 1-j;
  z(2)= 1+j;
  z(3)=-1+j;
  z(4)=-1-j;
  z(5)= 0-3*j;
  z(6)= 3-3*j;
  z(7)= 3;
  z(8)= 3+3*j;
  z(9)= 0+3*j;
  z(10)=-3+3*j;
  z(11)=-3;
  z(12)=-3-3*j;
  z(13)=-5*j;
  z(14)= 5;
  z(15)=5*j;
  z(16)= -5;
  d0=sign(randn(1,N));
  d1=sign(randn(1,N));
  d2=sign(randn(1,N));
  d3=sign(randn(1,N));
  d0=0.5*(d0+ones(1,N));
  d1=0.5*(d1+ones(1,N));
  d2=0.5*(d2+ones(1,N));
  d3=0.5*(d3+ones(1,N));
  d=d0+2*d1+4*d2+8*d3;
  d=d+ones(1,N);
  d=z(d);
  %d=d/sqrt(54/4);   % Normierung der Leistung auf eins

end;

% ### EOF ######################################################################
