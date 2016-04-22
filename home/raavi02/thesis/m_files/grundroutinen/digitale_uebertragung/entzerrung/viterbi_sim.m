% ##############################################################################
% ##  viterbi_sim.m : Simulationsprogramm fuer Viterbi-Entzerrer              ##
% ##                  Erzeugung von Fehlervektoren und Symbolfehlerrate       ##
% ##############################################################################
%
% Aufruf:    [Ps,de] = viterbi_sim(f,EbN0,N);
%
% Eingabe:   h    = Kanalimpulsantwort
%            EbN0 = Eb/N0-Verhaeltnis in dB
%            N    = Anzahl von 8192-Datenbloecken
%
% Ausgabe:   Ps   = Symbolfehlerwahrscheinlichkeit
%            de   = (d - d_dach)/d_min
%
%                                      Maerz 2001, Heiko Schmidt, K.D. Kammeyer


function [Ps,de]=viterbi_sim(h,EbN0,N);

% ------------ Auswuerfeln der Symbole -----------
%
% wichtig: Einschwingen des Kanals mit Symbol Nr. 0
%          Terminierung der Sequenz mit Symbol Nr. 0

h=h(:);
h=h/sqrt(h'*h);
all_sym = [(-1 -j) (1 -j) (1 +j) (-1 +j)].'/sqrt(2);  % Symbolalphabet
n_data  = 8192-length(h)+1;

gamma   = 1/sqrt(2)*10^(-EbN0/20)/sqrt(2);    %  incl.  Eb  --> Es-Umrechnung
Ps=0;
disp('  ')
disp('  ')

disp('bitte etwas Geduld ...');
disp('  ')

for i=1:N
  %disp(sprintf('Schleife %d ',i))       funktioniert nicht da, wo's soll!!!!
  Schleife=i
  sym_num = [zeros(length(h)-1,1);...
             floor(length(all_sym) * rand(n_data,1));...
             zeros(length(h)-1,1)];
  s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
  %
  r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
  %
  % ---- Rauschen
  r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));
  %
  %
  s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
  r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden

  % ------------ Viterbi-Aufruf -------------
  term_flag   = 1;
  r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag);

  Ps_i = mean(sign(abs(s_sym-r_sym)));

  Ps=Ps+Ps_i;
end
Ps=Ps/N;
% normierter Fehler:
de=1/sqrt(2)*(r_sym-s_sym);

% ### EOF ######################################################################
