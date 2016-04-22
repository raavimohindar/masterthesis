% ##############################################################################
% ##  ber_psk_ray.m : BER und SER fuer PSK beim Rayleigh Kanal,               ##
% ##                  koh. Detektion                                          ##
% ##############################################################################
%
% function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
%  = ber_psk_ray (EbN0dbin,M,L,map,maxbiterror,maxbitcount)
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
%   Symbolfehlerrate hinzugefuegt                        Juergen Rinas, 21.06.01
% ------------------------------------------------------------------------------

function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
  = ber_psk_ray (EbN0db,M,L,map,maxbiterror,maxbitcount)

% Defaultwerte
%%%%%%%%%%%%%%
if (nargin<2) M=2; end;
if (nargin<3) L=1; end;
if (nargin<4) map='gray'; end;
if (nargin<5) maxbiterror=1e4; end;
if (nargin<6) maxbitcount=100*maxbiterror; end;


% Plausibilitaetstests
%%%%%%%%%%%%%%%%%%%%%
K=log2(M);
if (K~=fix(K))
  disp(['# Achtung(',mfilename,') M ist keine Potenz von 2.']);
end;


disp([mfilename,' simuliere...']);
disp([' EbN0db=',num2str(EbN0db,'%8.3f')]);
disp([' M=',int2str(M)]);
disp([' L=',int2str(L)]);
disp([' map=',map]);
disp([' maxbiterror=',int2str(maxbiterror)]);
disp([' maxbitcount=',int2str(maxbitcount)]);


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
  symnoTX=fix(M*rand(symanz,1));
  BitTX=codetable(symnoTX'+1,:);
  sTX=alphabet(symnoTX+1);


  % Uebertragungskanal
  % symanz x L Kanalkoeffizienten auswuerfeln
  a=sqrt(gammas/2)*(randn(symanz,L)+j*randn(symanz,L));
  % Symbole sTX ueber L-Rayleigh-Kanal uebertragen
  Rauschen=sqrt(1/2)*(randn(symanz,L)+j*randn(symanz,L));
  sRX=(sTX*ones(1,L)).*a+ Rauschen;

  % Empfaenger
  sRX=sum(sRX.*conj(a),2); % Maximum Ratio Combining
  arcRX=angle(sRX);
  arcRX(arcRX<0)=arcRX(arcRX<0)+2*pi;
  symnoRX=fix(arcRX/2/pi*M);
  BitRX=codetable(symnoRX'+1,:);


  % Bitvergleich
  biterror=biterror+sum(BitTX~=BitRX,1);
  % Symbolvergleich
  symerror=symerror+sum(symnoTX~=symnoRX);
  symcount=symcount+symanz;
  bitcount=bitcount+symanz*K;
end;


ber=biterror/symcount;
meanber=mean(ber);
ser=symerror/symcount;

% ### EOF ######################################################################