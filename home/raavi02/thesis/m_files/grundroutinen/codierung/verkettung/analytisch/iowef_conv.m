% ##############################################################################
% ##  iowef_conv.m : Distanzfunktion eines Faltungscode berechnen             ##
% ##############################################################################
%
% function [dfuncwh]=iowef_conv(Cl,n,k,Gin,Rin,dmax,jmax,punctin)
%
% EINGABE:
%   Cl: Constraint-length
%       (Skalar)
%
%    n: Anzahl der Codebits
%       (Skalar)
%
%    k: Anzahl der Infobits
%       (Skalar)
%
%  Gin: Spaltenvektor mit den Generatorpolynomen
%         dezimal codiert (1+D+D^3 -> 11) als Zahl
%         oktal codiert   (1+D+D^3 -> 15) als string
%       (Spaltenvektor)
%
%  Rin: Spaltenvektor mit den Polynomen der Rueckkopplung
%         dezimal codiert (1+D+D^3 -> 11) als Zahl
%         oktal codiert   (1+D+D^3 -> 15) als string
%       (Spaltenvektor)
%       default: keine Rueckkopplung
%
% dmax: bis zur Distanz dmax berechnen
%       (Skalar)
%       default: 30
%
% jmax: maximale Anzahl von Trellissegmenten, um Abbruch der Berechnung zu
%       erzwingen
%       (Skalar)
%       default: unendlich
%
% punctin: Punktierungsmuster des Code als Spaltenvektor
%          Je Zeile ein Punktierungsmuster fuer ein Codewort
%          dezimal codiert, als Zahl
%          (2 -> [1 0] = 2^1+0*2^0 -> Bit c_0 punktiert, c_1 uebertragen) 
%          oktal codiert, als String, Zuordnung siehe dezimal 
%          default: keine Punktierung
%
% AUSGABE:
% dfuncwh: Distanzfunktion
%          w=Infowortgewicht (Matrix-Index 1 entspricht w=0 !!!)
%          h=Codewortgewicht (Matrix-Index 1 entspricgt h=0 !!!)
%         (Matrix)
%
% ANMERKUNGEN:
%   - Interface zur mex-Funktion distfunc_conv_mex
%
%
%   Es existieren zwei Bedingungen, die Alternativ zum Abbruch fuehren koennen:
%    1.) Die Anzahl der zu beruecksichtigenden Trellissegmente kann
%        auf jmax begrenzt werden.
%    2.) Falls bei der Berechnung kein zusaetzliches Trellissegment die
%        Distanzen <= dmax veraendern wird.
%
% AUTOR: Juergen Rinas,  31.05.1999
%                        10.11.2000
%        erweitert Volker Kuehn, 13.09.2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dfuncwh]=iowef_conv(Cl,n,k,Gin,Rin,dmax,jmax,punctin)

% Default-Werte vergeben
if (nargin<5) Rin=zeros(size(Gin,1),1); end;
if (nargin<6) dmax=30; end;
if (nargin<7) jmax=inf; end;
if (nargin<8) punctin=(2.^n)-1; end;

if (isstr(Gin)) 
  G=base2dec(Gin,8);           % oktale Eingabe von G und R ermoeglichen
else 
  G=Gin;
  ptr=find(Gin~=1);
  G(ptr)=bi2de(fliplr(de2bi(Gin(ptr),Cl))); 
end;
if (isstr(Rin)) 
  R=base2dec(Rin,8);           % oktale Eingabe von G und R ermoeglichen
else 
  R=Rin;
  ptr=find(Rin~=1);
  R(ptr)=bi2de(fliplr(de2bi(Rin(ptr),Cl))); 
end;
if (isstr(punctin)) 
  punct=base2dec(punctin,8);   % oktale Eingabe von G und R ermoeglichen
else 
  punct=punctin; 
end;

[dfuncwh]=iowef_conv_mex(Cl,n,k,G,R,dmax,jmax,punct);

% ### EOF ######################################################################
