% ##############################################################################
% ##  iowef_spc.m : Berechnung der IOWEF eines Parity-Check-Codes             ##
% ##############################################################################
%
% function [ACwh]=iowef_spc(n)
% ------------------------------------------------------------------------------
% EINGABE:
%        n: Laenge der Codeworte
%           (Skalar)
%
% AUSGABE:
% ACwh: Input Output Weight Enumerating Function
%       w=Infowortgewicht (Matrix-Index 1 entspricht w=0 !!!)
%       h=Codewortgewicht (Matrix-Index 1 entspricgt h=0 !!!)
%       (Matrix)
%
% ANMERKUNGEN:
%   - Infowortlaenge k=n-1
%   - benoetigt Datei bin_coef.m
%
% AUTOR: Juergen Rinas,  31.05.1999
% ------------------------------------------------------------------------------

function [ACwh]=iowef_spc(n)

k=n-1;

ACwh=zeros(k+1,n+1);
b=bin_coef(k*ones(k+1,1),0:k);
for w=0:k
  ACwh(w+1,2*ceil(w/2)+1)=b(w+1);
end;

% ### EOF ######################################################################
