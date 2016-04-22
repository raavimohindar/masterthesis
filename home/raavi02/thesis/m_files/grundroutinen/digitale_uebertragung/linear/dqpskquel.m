% ##############################################################################
% ##  dqpskquel.m : Erzeugung eines DQPSK-Signals im Symboltakt               ##
% ##############################################################################
%
% Aufruf:    d_dqpsk=dqpskquel(b,typ);
%
% Eingabe:   b   = binaerer Quellbitvektor gerader Laenge ('0', '1')
%            typ = 1 --> lambda=pi/4, Gray-Zuordnung
%            typ = 2 --> lambda=0,    sequentielle Bitzuordnung
%            typ = 3 --> lambda=0,    Gray-Zuordnung
%            typ = 4 --> lambda=pi/4, sequentielle Bitzuordnung
%
% Ausgabe:   d_dqpsk = komplexer QPSK-Signalvektor
%
% Anmerkung: Die Datei qpskquel.m wird benoetigt
%
%                                                   Kammeyer '90

function d_dqpsk=dqpskquel(b,typ);

d_qpsk=qpskquel(b,typ);

dphi=angle(d_qpsk)/pi;
phi=[0 cumsum(dphi)];
d_dqpsk=exp(j*pi.*phi);

% ### EOF ######################################################################
