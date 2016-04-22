% ##############################################################################
% ##  ber_dpsk_ink_awgn.m : BER fuer DPSK beim AWGN Kanal,                    ##
% ##                        inkohaerente Detektion                            ##
% ##############################################################################
%
% function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
%  = ber_dpsk_ink_awgn(EbN0db,M,map,maxbiterror,maxbitcount)
% ------------------------------------------------------------------------------
% EINGABE:
%      EbN0db: Signal/Noise pro Infobit in db
%              (Skalar)
%           M: Stufigkeit der Modulation (default=2)
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
% QUELLE:
%   [Kam96,S.378f.]
%
%
% AUTOR:
% * Juergen Rinas,  08.12.1998
% * aufgeteilt in Simulations- und Speicherroutine,
%   Symbolfehlerrate hinzugefuegt - Juergen Rinas, 21.06.01
% ------------------------------------------------------------------------------

function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
  = ber_dpsk_ink_awgn(EbN0db,M,map,maxbiterror,maxbitcount)


% Defaultwerte
%%%%%%%%%%%%%%
if (nargin<2) M=2; end;
if (nargin<3) map='gray'; end;
if (nargin<4) maxbiterror=1e4; end;
if (nargin<5) maxbitcount=100*maxbiterror; end;


disp([mfilename,' simuliere...']);
disp([' EbN0db=',num2str(EbN0db,'%8.3f')]);
disp([' M=',int2str(M)]);
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

% Mapping fuer DPSK
if (upper(map(1))=='N')
  codetable=nattable(K);
else
  codetable=graytable(K);
end

EbN0=10^(EbN0db/10);
EsN0=K*EbN0;

symcount=0;
bitcount=0;
symerror=0;
biterror=zeros(1,K);

alphabet=exp(j*(pi*2/M*[0:M-1]'+pi/M));


% Simulationsschleife
%%%%%%%%%%%%%%%%%%%%%
while ((biterror<maxbiterror) & (bitcount<maxbitcount))
  % Realisierung eines Demodulators auf Sender- und Empfangsseite, da
  % es zu aufwendig ist einen Bitvector in einen
  % Symboldifferenzvektor umzurechnen


  % Sender
  symnoTX=fix(M*rand(symanz+1,1));
  sTX=alphabet(symnoTX+1);
  sdemodTX=sTX(2:symanz+1,1).*conj(sTX(1:symanz,1));
  arcTX=angle(sdemodTX.*exp(j*pi/M)); % Signalraum um pi/m drehen, damit
  %                                     Entscheidergrenzen nicht auf die
  %                                     Re- oder Im-Achse des Signalraumes
  %                                     fallen
  arcTX(arcTX<0)=arcTX(arcTX<0)+2*pi;
  symTX=fix(arcTX/2/pi*M);
  BitTX=codetable(symTX'+1,:);


  % Uebertragungskanal
  Rauschen=sqrt(1/2)*(randn(symanz+1,1)+j*randn(symanz+1,1));
  sRX=sqrt(EsN0)*sTX+Rauschen;


  % Empfaenger
  sdemodRX=sRX(2:symanz+1,1).*conj(sRX(1:symanz,1));
  arcRX=angle(sdemodRX.*exp(j*pi/M));
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