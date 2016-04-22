% ##############################################################################
% ##  ber_dpsk_ink_awgn_save.m : Speichern und Laden von                      ##
% ##                             Simulationsergebnissen der Routine           ##
% ##                             ber_dpsk_ink_awgn.m                          ##
% ##############################################################################
%
% function [meanber, ber, ser, biterror, bitcount, symerror, symcount] = ...
%  ber_dpsk_ink_awgn_save(EbN0db,M,map,maxbiterror,maxbitcount)
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
%      ber: Bitfehlerrate (pro bit)
%           (Zeilenvektor)
%      ser: Symbolfehlerrate
%           (Skalar)
% biterror: Anzahl der Bitfehler (pro bit)
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
%   Simulationsergebnisse werden im Unterverzeichnis ./result/
%   gespeichert!
%
%
% AUTOR:
% * Juergen Rinas,  08.12.1998
% * Erweitert um bit-spezifische Fehlerraten - Volker Kuehn, 04.01.01
% * aufgeteilt in Simulations- und Speicherroutine - Juergen Rinas, 21.06.01
% ------------------------------------------------------------------------------

function [meanber, ber, ser, biterror, bitcount, symerror, symcount] = ...
  ber_dpsk_ink_awgn_save (EbN0db,M,map,maxbiterror)

global tempdir;

% Defaultwerte
%%%%%%%%%%%%%%
if (nargin<2) M=2; end;
if (nargin<3) map='gray'; end;
if (nargin<4) maxbiterror=1e4; end;
if (nargin<5) maxbitcount=100*maxbiterror; end;


% Plausibilitaetstests
%%%%%%%%%%%%%%%%%%%%%
K=log2(M);
if (K~=fix(K))
  disp(['# Achtung(',mfilename,') M ist keine Potenz von 2.']);
end;


% bereits berechnete Werte laden, falls Datei vorhanden!
% ansonsten Simulation starten
matfilename=[tempdir,'/results/',mfilename,...
             '-EbN0=',num2str(EbN0db),...
             '-M=',int2str(M),...
             '-map=',upper(map(1)),...
             '-mbe=',int2str(maxbiterror),...
             '-mbc=',int2str(maxbitcount),...
             '.mat'];
if exist(matfilename)~=0
  load(matfilename);
else
  [meanber, ber, ser, biterror, bitcount, symerror, symcount] = ...
    ber_dpsk_ink_awgn(EbN0db,M,map,maxbiterror,maxbitcount);
  save(matfilename,'meanber','ber','ser','biterror','bitcount',...
       'symerror','symcount');
end;

% ### EOF ######################################################################
