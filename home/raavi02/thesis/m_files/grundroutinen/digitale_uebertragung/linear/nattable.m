% ##############################################################################
% ##  nattable.m : Tabelle natuerliches Mapping erstellen                     ##
% ##############################################################################
%
% function [nt] = nattable(K)
% ------------------------------------------------------------------------------
% EINGABE:
%      K: Stellen in der NAT-Code-Tabelle
%         (Skalar)
%
% AUSGABE:
%   gt: NAT-Code-Tabelle
%       ( (2^K,K)-Matrix  )
%
% BEMERKUNGEN:
%   benoetigt Routine de2bi der communications library von Mathworks
%
% BEISPIEL:
%    » nattable(2)
%    ans =
%       0     0
%       0     1
%       1     0
%       1     1
%
% AUTOR: Volker Kuehn,  04.01.2001
% ------------------------------------------------------------------------------

function [nt] = nattable(K)

nt=0:2^K-1;
nt=de2bi(nt(:));

% ### EOF ######################################################################
