% ##############################################################################
% ## Grundroutinen zur Entzerrung                                             ##
% ##############################################################################
%
% dfe.m                   : MSE-Loesung fuer FIR-DF-Entzerrer
%
% ez_sim.m                : Simulation eines Entzerrers mit Uebergabe der Kanal-
%                            und Entzerrer-Impulsantwort (wahlweise w=1 und w=2)
%
% fir_dfe_entwurf.m       : MSE-Loesung fuer FIR-DF-Entzerrer
%
% frac_ez_entwurf.m       : Berechnung eines T/2 Entzerrers aus  Kanalimpuls-
%                            antwort mit Nebenbedingung minimaler
%                            Koeffizientenenergie
%
% lms.m                   : Adaptiver Entzerrer mit iterativer Einstellung nach
%                             dem LMS-Algorithmus
%
% t_ez_entwurf.m          : Zero-Forcing oder MMSE-Entwurf der Impulsantwort
%                            eines T-Entzerrers
%
% viterbi_demo.m          : Demonstrationsprogramm fuer 'viterbi_entzerrer.m'
%
% viterbi_entzerrer.m     : Viterbi-Entzerrung
%
% viterbi_fehlervec.m     : Fehleranalyse des Viterbi-Entzerrers fuer QPSK
%
% viterbi_sim.m           : Simulationsprogramm fuer Viterbi-Entzerrer
%                            Erzeugung von Fehlervektoren und Symbolfehlerrate
%
% ### EOF ######################################################################
