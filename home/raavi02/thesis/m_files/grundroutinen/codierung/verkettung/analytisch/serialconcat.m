% ##############################################################################
% ##  serialconcat.m : Serielle Verkettung von Codes Uniform-Interleaver      ##
% ##                   (Laenge N)                                             ##
% ##############################################################################
%
% function ACwhs=serialconcat(ACwho,ACwhi,N)
% ------------------------------------------------------------------------------
% EINGABE:
%  Awho,ACwhi: Matrix der Input Output Weight Enumerating Function
%              w=Infowortgewicht (Matrix-Index 1 entspricht w=0 !!!)
%              h=Codewortgewicht (Matrix-Index 1 entspricgt h=0 !!!)
%                o=outer Code, i=inner Code
%              (Matrix)
%  N:          Interleaverlaenge
%              (Skalar)
%
% AUSGABE:
%   ACwhs: Matrix der resultierenden IOWEF
%          (Matrix)
%
% ANMERKUNGEN:
%   - benoetigt Datei bin_coef
%   - Matrizen werden in ihrer Groesse angepasst, um Verknuepfung zu
%     ermoeglichen
%
% QUELLE:
%   [BDMP96,S.4,(4)]
%
% AUTOR: Juergen Rinas,  31.05.1999
% ------------------------------------------------------------------------------

function [ACwhs]=serialconcat(ACwho,ACwhi,N)

hmaxo=size(ACwho,2);
wmaxi=size(ACwhi,1);

if (hmaxo>wmaxi)
  disp([mfilename,': ausseres Codegewicht > inneres Infogewicht ',...
        '(hmaxo>wmaxi)']);
  ACwho(:,wmaxi+1:hmaxo)=[];
end;
if (wmaxi>hmaxo)
  disp([mfilename,': inneres Infogewicht > ausseres Codegewicht ',...
        '(wmaxi>hmaxo)']);
  ACwhi(hmaxo+1:wmaxi,:)=[];
end;

hmaxo=size(ACwho,2);
wmaxi=size(ACwhi,1);

b=bin_coef(N*ones(hmaxo,1),0:hmaxo-1);

ACwhs=((ACwho ./repmat(b',size(ACwho,1),1)) * ACwhi);

% ### EOF ######################################################################
