% ##############################################################################
% ##  dpskdem.m : DPSK-Demodulation                                           ##
% ##############################################################################
%
% kohaerent:    (DECPSK),
%               --> Phasenregelung u. Uebergabe des entschiedenen Phasensterns
% inkohaerent:  (Differentieller Demodulation)
%               --> Uebergabe des ungeregelten verrauschten Pahsenensterns
%
% Aufruf:       sdem=dpskdem(sdpsk);
%
% Eingabe:      sdpsk = entschiedener Phasenstern bei DECPSK
%               sdpsk = gestoerte DPSK-Symbolfolge bei inkohaerenter
%                       Demodulation
%
% Ausgabe:      sdem = demoduliertes PSK-Signal
%               zur Entscheidung bzw. Bitzuordnung -> qpskdec.m

function sdem=dpskdem(sdpsk);

s1=sdpsk;
s1(1)=[];
s2=sdpsk;
s2(length(sdpsk))=[];
sdem=s1.*conj(s2);

% ### EOF ######################################################################
