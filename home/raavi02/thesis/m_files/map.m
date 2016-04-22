% ##############################################################################
% ## map.m : SS-MAP-Decodierung fuer RSC- und NSC-Codes Codes                 ##
% ##         liefert LLR's fuer Infobit und fuer Codebit (optional)           ##
% ##############################################################################
%
% function [L_info,L_code]=map(signal, code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EINGABE:
%  signal              : Struct mit Angaben zum Eingangssignal
%    signal.sig        : empfangenes Signal (Vektor mit LLR's)
%    signal.L_a        : a-priori-Information fuer jedes Informationsbit
%    signal.last_state : letzter Zustand des Trellisdiagramms
%                        (relevant, falls code.term==1)
%
%  code                : Struct mit Angaben zum Code
%    code.trellis_out  : beschreibt Codeworte der einzelnen Zustandsuebergaenge
%    code.trellis_next : beschreibt Uebergaenge zwischen den einzelnen
%                        Zustaenden
%    code.num_state    : Anzahl der Zustaende im Trellis
%    code.block_len    : Anzahl der Codeworte pro Block
%    code.word_len     : Anzahl der Bits je Codewort
%    code.term         : terminiertes Trellisdiagramm: term == 1, sonst term~=1
%
% Ausgabe:
%  L_info              : LLR's fuer decodierte Informationsbits
%  L_code              : LLR's fuer Codebits (optional)
%
% ANMERKUNGEN:
%
%
% AUTOR: Volker Kuehn, 29.08.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ### EOF ######################################################################
