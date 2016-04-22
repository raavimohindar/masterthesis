% ##############################################################################
% ##  gen_punc.m : generiert Punktierungsvektoren fuer Turbo-Codes            ##
% ##############################################################################
%
% function punc_vector = gen_punc(pattern,N,num_codes)
% ------------------------------------------------------------------------------
% EINGABE:
%       pattern:   Matrix mit dem Punktierungsschema fuer eine Periode
%                  (in jeder Spalte zuerst die system. Bit, dann Parity Bit der
%                   anderen Codes, Anzahl der Spalten gleich der
%                   Punktierungsperiode)
%       N:         Anzahl der Codeworte pro Rahmen (Interleaver Groesse)
%
% AUSGABE:
%       punc_vector: Vektor mit den Positionen der zu uebertragenen Bit im
%                    unpunktierten Datenstrom
%-------------------------------------------------------------------------------
function punc_vector = gen_punc(pattern,N)

[n,period] = size(pattern);

punc_vector = repmat(pattern(:),ceil(N/period),1);

% Abschneiden auf Blocklaenge
punc_vector = punc_vector(1:N*n);

punc_vector = find(punc_vector);

% ### EOF ######################################################################
