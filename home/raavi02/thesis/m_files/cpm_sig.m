% ##############################################################################
% ##  cpm_sig.m : komplexe Einhuellende eines CPM-Signals im Abtasttakt To    ##
% ##############################################################################
%
%  cpm_sig bestimmt die komplexe Einhuellende s eines CPM-Signals im
%  Abtasttakt To.  Die charakteristischen Parameter der CPM-Art werden
%  vorher durch cpm_ini festgelegt und mittels struct_cpm_ini uebergeben.
%  Der Datenvektor d im Symboltakt wird
%  mit ausgegeben.
%
% Aufruf:    s = cpm_sig(d);
%
% Beispiel:  struct_cpm_ini = cpm_ini(L,'GAUSS',eta,M,w,f(3dB)*T)
%            s = cpm_sig(d, struct_cpm_ini);
%
% Eingabe:   d = Vektor der zu uebertragenen Daten, Laenge m;
%                kann z.B. mit 'cpmdaten' erzeugt werden
%                m soll  ganzahliges Vielfaches von L=length(q)/w sein.
%            struct_cpm_ini = CPM-Parameter aus cpm_ini
%
% Ausgabe:   s = komplexe Einhuellende
%
% Benoetigte m-Files: cpm_ini, qpuls, cpm_mod
%
%                                                      Benthin 9/91

function s = cpm_sig(d, struct_cpm_ini);

versatz = 0;
isis = 0;
q = qpuls(struct_cpm_ini.L,struct_cpm_ini.impulsform,...
          struct_cpm_ini.w,versatz,struct_cpm_ini.f3dBT);
q = [0 q]; % Korrektur Ka: Phasenanfang null
q(length(q))=[];
s = cpm_mod(q,d,struct_cpm_ini.eta,struct_cpm_ini.M,struct_cpm_ini.w,isis);

% ### EOF ######################################################################
