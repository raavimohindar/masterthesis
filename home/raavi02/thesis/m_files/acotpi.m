% ##############################################################################
% ##  acotpi.m: Berechnung des Arkuskotangens mit Wertebereich [0.. pi]       ##
% ##            (Hauptwert)                                                   ##
% ##############################################################################
%
% function [y] = acotpi(z)
% ------------------------------------------------------------------------------
% EINGABE:
%      z: Argument
%         (Skalar, Vektor oder Matrix)
%
% AUSGABE:
%      y: Arkuskotangens
%         (Skalar, Vektor oder Matrix)
%
% ANMERKUNGEN:
%   - Diese Funktion liefert den HAUPTWERT des Arkuscotangens
%   - Die acot-Funktion von Matlab hat einen Wertebereich von
%    [-pi/2.. pi/2] und macht eine Division durch Null bei acot(0).
%
% QUELLE:
%   [BSMM95,S.73,(2.148)]
%
% AUTOR: Juergen Rinas,  08.12.1998
% ------------------------------------------------------------------------------

function y = acotpi(z)

y = pi/2-atan(z);

% ### EOF ######################################################################
