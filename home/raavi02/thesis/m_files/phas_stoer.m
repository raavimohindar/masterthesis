% ##############################################################################
% ##  phas_stoer.m : Drehung eines Vektors mit der Geschwindigkeit do*iT      ##
% ##############################################################################
%
% Aufruf:    function y_channel = Phas_stoer(d_quel,dfT,theta,dphi,fjT,EbNo);
%
% Eingabe:   d_quel    = QPSK-Quelldaten
%            dfT       = Frequenzoffset Delta_f*T
%            theta     = Phasenoffset (RAD)
%            dphi      = Phasenhub (RAD) des Jitters
%            fjT       = Jitterfrequenz f1*T
%            EbNo      = E_b/N_0 des komplexen Rauschens
%
% Ausgabe:   y_channel = Kanalausgang im Symboltakt
%
%                                                   Kammeyer

function y_channel = phas_stoer(d_quel,dfT,theta,dphi,fjT,EbNo);

M=4;        % Stufigkeit QPSK
Snutz = 1;  % Leistung Datum
l_s = length(d_quel);

% Phasen- und Frequenzoffset einfuegen
k = 0:l_s -1;
dreh = exp(j.*(2*pi*dfT.* k + theta + dphi.*cos(2*pi*fjT.*k) ));
y_channel = d_quel(:) .* dreh(:);

% Rauschen hinzufuegen
EsN0=10^(EbNo/10)*log2(M);
sigma=Snutz/EsN0;
sigma1=sigma/2;
skal=sqrt(sigma1);
r=skal*(randn(l_s,1)+j*randn(l_s,1));
y_channel = y_channel + r;

% ### EOF ######################################################################
