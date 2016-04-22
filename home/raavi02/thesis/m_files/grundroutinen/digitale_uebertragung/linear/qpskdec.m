% ###########################################################################
% ##  qpskdec.m : QPSK-Entscheidung und Bitzuordnung                       ##
% ###########################################################################
%
% Aufruf:    [sdach,bdach]=qpskdec(sqpsk,typ);
%
% Eingabe:   sqpsk = entschiedenes QPSK-Signal bei kohaerenter Demodulation
%                    (DECPSK)
%            sqpsk = gestoertes QPSK-Signal bei inkohaerenter Demodulation
%            typ = 1 --> lambda=pi/4, Gray-Zuordnung
%            typ = 2 --> lambda=0,    sequentielle Bitzuordnung
%            typ = 3 --> lambda=0,    Gray-Zuordnung
%            typ = 4 --> lambda=pi/4, sequentielle Bitzuordnung
%
% Ausgabe:   sdach = entschiedenes QPSK-Signal
%            bdach = zugeordnete Bits

function [sdach,bdach]=qpskdec(sqpsk,typ);

L=2*length(sqpsk);  %  Laenge der Bitfolge

if typ == 1
  sdach = sign(real(sqpsk))+j*sign(imag(sqpsk));
  bdach(1:2:L)= 0.5*(-real(sdach) + 1);
  bdach(2:2:L)= 0.5*(-imag(sdach) + 1);
  sdach=sdach/abs(sdach(1));

elseif typ==2
  sdach=sqpsk*exp(j*pi/4);
  s=sign(real(sdach))+j*sign(imag(sdach));
  sdach=s*exp(-j*pi/4);
  sdach=sdach/abs(sdach(1));
  phi=angle(s)*2/pi;
  phi=0.5*(sign(phi)+1).*phi - 0.5*(sign(phi)-1).*(4+phi);
  phi=phi-0.5;
  bdach(1:2:L)=fix(phi/2);
  bdach(2:2:L)=phi-2*bdach(1:2:L);

elseif typ==3
  sdach = sqpsk*exp(j*pi/4);
  sdach = sign(real(sdach))+j*sign(imag(sdach));
  bdach(1:2:L)= 0.5*(-real(sdach) + 1);
  bdach(2:2:L)= 0.5*(-imag(sdach) + 1);
  sdach=sdach*exp(-j*pi/4);
  sdach=sdach/abs(sdach(1));

elseif typ==4
  sdach=sqpsk;
  sdach=sign(real(sdach))+j*sign(imag(sdach));
  sdach=sdach/abs(sdach(1));
  phi=angle(sdach)*2/pi;
  phi=0.5*(sign(phi)+1).*phi - 0.5*(sign(phi)-1).*(4+phi);
  phi=phi-0.5;
  bdach(1:2:L)=fix(phi/2);
  bdach(2:2:L)=phi-2*bdach(1:2:L);

end;

% ### EOF ######################################################################
