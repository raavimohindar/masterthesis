% ##############################################################################
% ##  qpskquel.m : Erzeugung eines QPSK-Signals im Symboltakt                 ##
% ##############################################################################
%
% Aufruf:    d_qpsk=qpskquel(b,typ);
%
% Eingabe:   b   = binaerer Quellbitvektor gerader Laenge ('0', '1')
%            typ = 1 --> lambda=pi/4, Gray-Zuordnung
%            typ = 2 --> lambda=0,    sequentielle Bitzuordnung
%            typ = 3 --> lambda=0,    Gray-Zuordnung
%            typ = 4 --> lambda=pi/4, sequentielle Bitzuordnung
%
% Ausgabe:   d_qpsk = komplexer QPSK-Signalvektor
%
%                                                   Kammeyer '90

function d_qpsk=qpskquel(b,typ);

% Fehlermeldung bei ungeradem b
if mod(length(b),2)
  disp('Fehler: Quellbitvektor b muss gerade sein, (2 Bit --> 1 Symbol)');
else
  % lambda=pi/4, Gray-Zuordnung
  if typ==1
    d_qpsk=-2*(b(1:2:length(b))-0.5 +j*(b(2:2:length(b))-0.5));

    % lambda=0, numerische Bitzuordnung
  elseif typ==2
    d_qpsk(1:length(b)/2)=2*b(1:2:length(b)) + b(2:2:length(b));
    d_qpsk=exp(j*pi/2.*d_qpsk);

    % lambda=0, Gray-Zuordnung
  elseif typ==3
    d_qpsk=-2*(b(1:2:length(b))-0.5 +j*(b(2:2:length(b))-0.5));
    d_qpsk=d_qpsk.*exp(-j*pi/4);

    % lambda=pi/4, numerische Bitzuordnung
  elseif typ==4
    d_qpsk(1:length(b)/2)=2*b(1:2:length(b)) + b(2:2:length(b));
    d_qpsk=exp(j*pi/2.*d_qpsk);
    d_qpsk=d_qpsk.*exp(j*pi/4);
  end
  d_qpsk=d_qpsk/abs(d_qpsk(1));
end;

% ### EOF ######################################################################
