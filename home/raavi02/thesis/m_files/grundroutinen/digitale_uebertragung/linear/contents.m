% ##############################################################################
% ## Grundroutinen zur linearen digitalen Modulation                          ##
% ##############################################################################
%
% acotpi.m                : Berechnung des Arcuskotangens mit 
%                            Wertebereich [0.. pi](Hauptwert)
%
% auge.m                  : erzeugt ein Augendiagramm
%
% awgn_sim.m              : Simulation einer AWGN-Uebertragung im komplexen TP-
%                           Bereich mit vorgegebenen Sende- und Empfangsfiltern
%
% ber_dpsk_ink_awgn.m     : BER fuer DPSK beim AWGN-Kanal, inkoh. Detektion
%
% ber_dpsk_ink_awgn_save.m: Speichern und Laden von Simulationsergebnissen
%                            der Routine ber_dpsk_ink_awgn.m
%
% ber_dpsk_ink_ray.m      : BER fuer DPSK beim Rayleigh-Kanal, inkoh. Detektion
%
% ber_dpsk_ink_ray_save.m : Speichern und Laden von Simulationsergebnissen
%                            der Routine ber_dpsk_ink_ray.m
%
% ber_dpsk_koh_awgn.m     : BER fuer DPSK beim AWGN Kanal, kohaerente Detektion
%
% ber_dpsk_koh_awgn_save.m: Speichern und Laden von Simulationsergebnissen
%                            der Routine ber_dpsk_koh_awgn.m
%
% ber_psk_awgn.m          : BER und SER bei kohaerenter PSK und AWGN-Kanal
%
% ber_psk_awgn_save.m     : Speichern und Laden von Simulationsergebnissen
%                            der Routine ber_psk_awgn.m
%
% ber_psk_ray.m           : BER und SER fuer kohaerente PSK beim Rayleigh-Kanal
%
% ber_psk_ray_save.m      : Speichern und Laden von Simulationsergebnissen
%                            der Routine ber_psk_ray.m
%
% bin_coef.m              : Berechnung des Binomialkoeffizienten n ueber k
%
% cosroll.m               : Erzeugung von Signal mit Cos-Roll-Off-Charakteristik
%
% datensig.m              : Erzeugung eines Datensignals mit Impulsformung
%
% dpskdem.m               : DPSK-Demodulation
%
% dqpskquel.m             : Erzeugung eines DQPSK-Signals im Symboltakt
%
% get_impulse             : Erzeugung von beliebig abgetasteten Impulsantworten
%
% graytable.m             : Tabelle mit Gray-Code erstellen
%
% modquel.m               : Datenquelle mit verschiedenen Modulationsformen
%                           (zufaellige Daten)
%
% nattable.m              : liefert Tabelle fuer natuerliches Mapping
%
% pb_dpsk_ink_awgn        : Bitfehlerwahrscheinlichkeit fuer inkohaerente DPSK
%                            beim AWGN-Kanal
%
% pb_dpsk_ink_ray         : Bitfehlerwahrscheinlichkeit fuer inkohaerente DPSK
%                            beim Rayleigh-Kanal
%
% pb_dpsk_koh_awgn        : Bitfehlerwahrscheinlichkeit fuer kohaerente DPSK
%                            beim AWGN-Kanal
%
% pb_psk_awgn             : Bitfehlerwahrscheinlichkeit fuer kohaerente PSK
%                            beim AWGN-Kanal
%
% pb_psk_ray.m            : Bitfehlerwahrscheinlichkeit bei kohaerente PSK
%                            beim Rayleigh-Kanal
%
% phas_stoer.m            : Drehung eines Vektors mit der Geschwindigkeit do*iT
%
% pmu_psk_awgn.m          : Bestimmung der paarweisen Fehlerwahrscheinlichkeiten
%                             beim AWGN-Kanal
%
% pr_decod.m              : Partial-Response-Decodierung
%
% pr_vorcod.m             : Partial-Response-Vorcodierung
%
% pskweight.m             : Ermittlung Vektor zur Gewichtung fuer die Berechnung
%                            der Bitfehlerwahrscheinlichkeit bei (D)PSK
%
% qpskdec.m               : QPSK-Entscheidung und Bitzuordnung
%
% qpskquel.m              : Erzeugung eines QPSK-Signals im Symboltakt
%
% simpson.m               : Integration nach der Simpson-Methode
%                           (quadratische Approximation)
%
% traegerreg.m            : Entscheidungsrueckgekoppelte Traegerphasenregelung
%                            fuer QPSK
%
% wurzcos.m               : Erzeugung eines Signals mit Wurzel-Cos-Roll-Off-
%                            Charakteristik
%
% ### EOF ######################################################################
