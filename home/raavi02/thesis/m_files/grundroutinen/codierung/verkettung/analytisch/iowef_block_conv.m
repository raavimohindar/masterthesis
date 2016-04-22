% ##############################################################################
% ##  iowef_block_conv.m : Distanzfunktion eines terminierten Faltungscode    ##
% ##                       berechnen, der als Blockcode interpretiert wird    ##
% ##############################################################################
%
% function [iowef]=iowef_block_conv(Cl,n,k,Gin,Rin,htrunc,wtrunc,wmax,punctin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% htrunc: Begrenzen auf Codegewicht htrunc
%         (Skalar)
%         default: 30
%
% wtrunc: Begrenzen auf Infogewicht wtrunc
%         (Skalar)
%         default: 30
%
%   wmax: Anzahl der Infobits des terminierten Faltungscode
%         (Skalar)
%         default: unbegrenzt
%
% punctin: Punktierungsmuster des Code als Spaltenvektor
%          Je Zeile ein Punktierungsmuster fuer ein Codewort
%          dezimal codiert, als Zahl
%          (2 -> [1 0] = 2^1+0*2^0 -> Bit c_0 punktiert, c_1 uebertragen) 
%          oktal codiert, als String, Zuordnung siehe dezimal 
%          default: keine Punktierung
%
% AUSGABE:
%   iowef: IOWEF
%          w=Infowortgewicht (Matrix-Index 1 entspricht w=0 !!!)
%          h=Codewortgewicht (Matrix-Index 1 entspricht h=0 !!!)
%         (Matrix)
%
% ANMERKUNGEN:
%   - Interface zur mex-Funktion iowef_block_conv_mex
%
% AUTOR: Juergen Rinas,  31.05.1999
%                        03.11.2000
%        erweitert Volker Kuehn, 13.09.2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iowef]=iowef_block_conv(Cl,n,k,Gin,Rin,htrunc,wtrunc,wmax,punctin)

% Default-Werte vergeben
if (nargin<5) Rin=zeros(size(Gin,1),1); end;
if (nargin<6) htrunc=30; end;
if (nargin<7) wtrunc=30; end;
if (nargin<8) wmax=inf; end;
if (nargin<9) punctin=(2.^n)-1; end;

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

[iowef]=iowef_block_conv_mex(Cl,n,k,G,R,htrunc,wtrunc,wmax,punct);

% ### EOF ######################################################################
