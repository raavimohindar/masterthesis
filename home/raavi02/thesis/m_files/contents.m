% ##############################################################################
% ## Grundroutinen zu Faltungscodes                                           ##
% ##############################################################################
%
% block_interleave.m : erzeugt Vektor fuer Block-Interleaving
%
% demultiplex.m      : Demultiplexer fuer turbo-codierte Sequenz mit RSC-Codes
%                       der Rate 1/n
%
% enc_pccc.m         : Turbo-Encoder fuer n parallel verkettete RSC-Codes
%
% enc_pcspc.m        : Codierung mit zwei parallel verketteten SPC-Codes
%
% enc_scspc.m        : Codierung mit zwei seriell verkettetem SPC-Codes
%
% gen_punc.m         : generiert Punktierungsvektoren fuer Turbo-Codes
%
% log_map.m          : BCJR-Algorithmus fuer Faltungscodes im logarithmischen
%                      Bereich (SS-Log-MAP)
%
% map.m              : BCJR-Algorithmus fuer Faltungscodes (SS-MAP)
%
% map_spc.m          : Soft-Output-Decodierung eines SPC-Codes
%
% max_log_map.m      : Naeherungsloesung des  BCJR-Algorithmus fuer Faltungscodes
%                      im log. Bereich (SS-Max-Log-MAP)
%
% spc.m              : Berechnung des Paritaetsbits fuer bin. Eingangsvektor
%
% turbo_pccc.m       : Turbo-Decodierung von parallel verketteten RSC-Codes,
%                      (RSC-Codes der Rate 1/n, diverse APP-Algorithmen,
%                       erster Code terminiert und nicht interleaved)
%
% turbo_scspc.m      : Turbo-Decodierung von seriell verketteten SPC-Codes
%
% turbo_sccc.m       : Turbo-Decodierung zweier seriell verketteten
%                      Faltungscodes, (aeusserer Code terminiert,
%                      Eingangssequenz nicht interleaved)
%
% turbo_pcspc.m      : Turbo-Decodierung von parallel verketteten SPC-Codes
%
% ### EOF ######################################################################
