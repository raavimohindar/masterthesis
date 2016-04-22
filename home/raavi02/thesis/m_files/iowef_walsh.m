% ##############################################################################
% ##  iowef_walsh.m : Berechnung der IOWEF eines Walsh-Codes                  ##
% ##############################################################################
%
% function [ACwh]=iowef_walsh(n)
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
%   - benoetigt Datei bin_coef.m
%   - Infowortlaenge k=log2(n)
%   - Beschreibung der Walsh-Codierung: [Kam96,S.639ff.]
%
% AUTOR: Juergen Rinas,  31.05.1999
% ------------------------------------------------------------------------------

function [ACwh]=iowef_walsh(n)

k=log2(n);

ACwh=zeros(k+1,n/2+1);
ACwh(:,n/2+1)=bin_coef(k*ones(k+1,1),0:k);
ACwh(1,n/2+1)=0;
ACwh(1,1)=1;

% ### EOF ######################################################################
