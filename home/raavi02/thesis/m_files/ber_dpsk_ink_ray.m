% ##############################################################################
% ##  ber_dpsk_ink_ray.m : BER fuer DPSK beim Rayleigh Kanal,                 ##
% ##                       inkoh. Detektion                                   ##
% ##############################################################################
%
% function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
%  = ber_dpsk_ink_ray (EbN0db,M,L,map,maxbiterror,maxbitcount)
% ------------------------------------------------------------------------------
% EINGABE:
%      EbN0db: Signal/Noise pro Infobit in db
%              (Skalar)
%           M: Stufigkeit der Modulation (default=2)
%              (Skalar)
%           L: Anzahl der Pfade (default=1)
%              (Skalar)
%         map: Art des Mapping  (dafault='gray')
%              'gray' = Gray-Codierung
%              'nat' = natuerliches Mapping
%              (String)
% maxbiterror: Abbruchkriterium: maximale Anzahl der Bitfehler (default=1e4)
%              (Skalar)
% maxbitcount: Abbruchkriterium: maximale Anzahl der Bits
%              (default=100*maxbiterror)
%              (Skalar)
%
%
% AUSGABE:
%   manber: Bitfehlerrate (gemittelt ueber alle Bits)
%           (Skalar)
%      ber: Bitfehlerrate (bit-spezifisch)
%           (Zeilenvektor)
%      ser: Symbolfehlerrate
%           (Skalar)
% biterror: Anzahl der Bitfehler (bit-spezifisch)
%           (Zeilenvektor)
% bitcount: Anzahl der Bits (Summe aller uebertragenen Bits)
%           (Skalar)
% symerror: Anzahl der Symbolfehler
%           (Slalar)
% symcount: Anzahl der Symbole
%           (Skalar)
%
%
% ANMERKUNGEN:
%   - Simulation!
%   - benoetigt Datei: graytable.m
%                      nattable.m
%
%
% AUTOR:
% * Juergen Rinas,  08.12.1998
% * aufgeteilt in Simulations- und Speicherroutine,
%   Symbolfehlerrate hinzugefuegt - Juergen Rinas, 21.06.01
% ------------------------------------------------------------------------------

function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
  = ber_dpsk_ink_ray (EbN0db,M,L,map,maxbiterror,maxbitcount)

% Defaultwerte
%%%%%%%%%%%%%%
if (nargin<2) M=2; end;
if (nargin<3) L=1; end;
if (nargin<4) map='gray'; end;
if (nargin<5) maxbiterror=1e4; end;
if (nargin<6) maxbitcount=100*maxbiterror; end;


disp([mfilename,' simuliere...']);
disp([' EbN0db=',num2str(EbN0db,'%8.3f')]);
disp([' M=',int2str(M)]);
disp([' L=',int2str(L)]);
disp([' map=',map]);
disp([' maxbiterror=',int2str(maxbiterror)]);
disp([' maxbitcount=',int2str(maxbitcount)]);


% Plausibilitaetstests
%%%%%%%%%%%%%%%%%%%%%
K=log2(M);
if (K~=fix(K))
  disp(['# Achtung(',mfilename,') M ist keine Potenz von 2.']);
end;


% Initialisierung
%%%%%%%%%%%%%%%%%
symanz=10000; % Anzahl der Symbole pro Simulationsschleifendurchlauf

% Mapping fuer PSK
if (upper(map(1))=='N')
  codetable=nattable(K);
else
  codetable=graytable(K);
end

EbN0=10^(EbN0db/10);
EsN0=K*EbN0;
gammas=EsN0/L;


symcount=0;
bitcount=0;
symerror=0;
biterror=zeros(1,K);

alphabet=exp(j*(pi*2/M*[0:M-1]'+pi/M));

% Simulationsschleife
%%%%%%%%%%%%%%%%%%%%%
while ((biterror<maxbiterror) & (bitcount<maxbitcount))
  % Sender
  symno1TX=fix(M*rand(symanz,1));
  symno2TX=fix(M*rand(symanz,1));
  s1TX=exp(j*pi*2/M*symno1TX);
  s2TX=exp(j*pi*2/M*symno2TX);
  sdemTX=s1TX.*conj(s2TX);
  arcTX=angle(sdemTX.*exp(j*pi/M)); % Signalraum um pi/m drehen, damit
  %                                   Entscheidungsgrenzen nicht auf die Re-
  %                                   oder Im-Achse des Signalraumes faellt
  %                                   und einfache Entscheidung moeglich ist
  arcTX(arcTX<0)=arcTX(arcTX<0)+2*pi;
  symTX=fix(arcTX/2/pi*M);
  BitTX=codetable(symTX'+1,:);


  % Uebertragungskanal
  % symanz * L Kanalkoeffizienten auswuerfeln
  a=sqrt(gammas/2)*(randn(symanz,L)+j*randn(symanz,L));
  % Symbole s1 und s2 ueber GLEICHEN Rayleigh-Kanal uebertragen
  % == LANGSAMES Fading
  Rauschen=sqrt(1/2)*(randn(symanz,L)+j*randn(symanz,L));
  s1RX=(s1TX*ones(1,L)).*a+Rauschen;
  Rauschen=sqrt(1/2)*(randn(symanz,L)+j*randn(symanz,L));
  s2RX=(s2TX*ones(1,L)).*a+Rauschen;


  % Empfaenger
  sRX=sum(s1RX.*conj(s2RX),2);
  arcRX=angle(sRX*exp(j*pi/M));
  arcRX(arcRX<0)=arcRX(arcRX<0)+2*pi;
  symRX=fix(arcRX/2/pi*M);
  BitRX=codetable(symRX'+1,:);


  % Bitvergleich
  biterror=biterror+sum(BitTX~=BitRX,1);
  % Symbolvergleich
  symerror=symerror+sum(symTX~=symRX);
  symcount=symcount+symanz;
  bitcount=bitcount+symanz*K;
end;


ber=biterror/symcount;
meanber=mean(ber);
ser=symerror/symcount;

% ### EOF ######################################################################
