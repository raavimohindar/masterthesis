% ##############################################################################
% ## Grundroutinen zu linearen Blockcodes                                     ##
% ##############################################################################
%
% gf_dft.m    : diskrete Fouriertransformation auf Galoisfeldern (GF2^m)
%
% rschien.m   : Chien-suche, findet die Fehlerstellen eines RS-Codewortes
%               (Suche der Nullstellen im Fehlerstellenpolynom -> rselp.m)
%
% rscorrect.m : Subtrahiert den geschaetzten Fehler vom Empfangswort
%               bei RS-Codes
%
% rselp.m     : Berechnung der Fehlerstellenpolynoms bei Reed-Solomon Codes
%               (elp=error location polynomial)
%
% rserramp.m  : Berechnung der Fehlerstellenamplituden des Syndroms bei RS-Codes
%
% rssyndrom.m : Berechnung des Syndroms bei Reed-Solomon Codes
%
% ### EOF ######################################################################
