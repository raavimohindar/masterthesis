% ##############################################################################
% ##  ber_dpsk_ink_ray_save.m : Speichern und Laden von                       ##
% ##                       Simulationsergebnissen der Routine                 ##
% ##                       ber_dpsk_ink_ray.m                                 ##
% ##############################################################################
%
% function [meanber, ber, ser, biterror, bitcount, symerror, symcount] ...
%  = ber_dpsk_ink_ray_save (EbN0dbin,M,L,map,maxbiterror,maxbitcount)
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
%   Simulationsergebnisse werden im Unterverzeichnis ./result/
%   gespeichert!
%
%
% AUTOR:
% * Juergen Rinas,  08.12.1998
% * aufgeteilt in Simulations- und Speicherroutine,
%   Symbolfehlerrate hinzugefuegt                        Juergen Rinas, 21.06.01
% ------------------------------------------------------------------------------

function [meanber, ber, ser, biterror, bitcount, symerror, symcount] = ...
  ber_dpsk_ink_ray_save (EbN0db,M,L,map,maxbiterror)

global tempdir;

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


% bereits berechnete Werte laden, falls Datei vorhanden!
% ansonsten Simulation starten
matfilename=[tempdir,'/results/',mfilename,...
             '-EbN0=',num2str(EbN0db),...
             '-M=',int2str(M),...
             '-L=',int2str(L),...
             '-map=',upper(map(1)),...
             '-mbe=',int2str(maxbiterror),...
             '-mbc=',int2str(maxbitcount),...
             '.mat'];
if exist(matfilename)~=0
  load(matfilename);
else
  [meanber, ber, ser, biterror, bitcount, symerror, symcount] = ...
    ber_dpsk_ink_ray(EbN0db,M,L,map,maxbiterror,maxbitcount);
  save(matfilename,'meanber','ber','ser','biterror','bitcount',...
       'symerror','symcount');
end;

% ### EOF ######################################################################
