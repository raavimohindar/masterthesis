% ##############################################################################
% ##  block_interleave.m : Block-Interleaver                                  ##
% ##############################################################################
%
%function IL = block_interleave(Zeilen,Spalten)
% ------------------------------------------------------------------------------
% EINGABE:
%        Zeilen: Anzahl der Zeilen der Matrix
%        Spalten: Anzahl der Spalten der Matrix
%
% AUSGABE:
%        IL: Vektor mit der Permutation
% ------------------------------------------------------------------------------
% Author: Volker Kuehn  17.03.2000

function IL = block_interleave( rows, cols )

N_total = rows*cols;
list    = reshape(1:N_total,rows,cols);
list    = list';
IL      = list(:);

% ### EOF ######################################################################
