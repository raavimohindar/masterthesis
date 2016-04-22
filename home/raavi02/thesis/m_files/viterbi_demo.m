% ##############################################################################
% ##  viterbi_demo.m : Demonstrationsprogramm fuer 'viterbi_entzerrer.m'      ##
% ##############################################################################
%
% Aufruf:    [t,g]=cosroll(r,w,L);
%
% Bedienung: Durch Wahl der Menuepunkte 1 ... 8 koennen folgende
%            Experimente durchgefuehrt werden:
%
%            1) BPSK,   Zufallskanal 2. Ordnung, Eingabe: Eb/N0, iv,
%               sukzessive graphische Ausgabe
%            2) BPSK,   Zufallskanal 3. Ordnung, Eingabe: Eb/N0, iv
%               sukzessive graphische Ausgabe
%            3) QPSK,   Zufallskanal 2. Ordnung, Eingabe: Eb/N0, iv
%               sukzessive graphische Ausgabe
%            4) QPSK,   Zufallskanal 3. Ordnung, Eingabe: Eb/N0, iv,
%               mit schneller grafischer Ausgabe
%            -------------------------------------------------------------------
%            5) BPSK,   waehlbarer Kanal 2. Ordnung; Eingabe: f, Eb/N0, iv,
%               Trellis-Terminierung ja/nein; sukzessive graphische Ausgabe
%            6) QPSK,   waehlbarer Kanal 2. Ordnung; Eingabe: f, Eb/N0, iv,
%               Trellis-Terminierung ja/nein; sukzessive graphische Ausgabe
%            -------------------------------------------------------------------
%            7) 16-QAM, Zufallakanal 1. Ordnung, Eingabe: Eb/N0, iv,
%               sukzessive graphische Ausgabe
%            8) 16-QAM, Zufallskanal 2. Ordnung, Eingabe: Eb/N0, iv,
%               mit schneller grafischer Ausgabe
%
% Anmerkung: benoetigt viterbi_entzerrer.m
%
%                                                   Heiko Schmidt, K.D. Kammeyer

auswahl = 10;
while auswahl>0
  disp(' ');
  disp(' ');
  disp(' ');
  disp('-------------------------------------------------------------------');
  disp(' ');
  disp(['1) BPSK,   Zufallskanal 2. Ordnung, Eingabe: Eb/N0, iv,   ',...
        'sukzessive graphische Ausgabe']);
  disp(['2) BPSK,   Zufallskanal 3. Ordnung, Eingabe: Eb/N0, iv    ',...
        'sukzessive graphische Ausgabe']);
  disp(['3) QPSK,   Zufallskanal 2. Ordnung, Eingabe: Eb/N0, iv    ',...
        'sukzessive graphische Ausgabe']);
  disp(['4) QPSK,   Zufallskanal 3. Ordnung, Eingabe: Eb/N0, iv,   ',...
        'mit schneller grafischer Ausgabe']);
  disp('-------------------------------------------------------------------');
  disp(['5) BPSK,   waehlbarer Kanal 2. Ordnung; Eingabe: f, Eb/N0, iv, ',...
        'Trellis-Terminierung ja/nein; sukzessive graphische Ausgabe']);
  disp(['6) QPSK,   waehlbarer Kanal 2. Ordnung; Eingabe: f, Eb/N0, iv, ',...
        'Trellis-Terminierung ja/nein; sukzessive graphische Ausgabe']);
  disp('-------------------------------------------------------------------');
  disp(['7) 16-QAM, Zufallakanal 1. Ordnung, Eingabe: Eb/N0, iv,   ',...
        'sukzessive graphische Ausgabe']);
  disp(['8) 16-QAM, Zufallskanal 2. Ordnung, Eingabe: Eb/N0, iv,   ',...
        'mit schneller grafischer Ausgabe']);
  disp(' ');
  disp('0) Ende ');
  disp(' ');
  disp('-------------------------------------------------------------------');
  disp(' ');
  %
  %
  auswahl_txt = input('Bitte Beispiel waehlen --> ','s');
  auswahl     = str2num(auswahl_txt);
  if length(auswahl) == 0
    auswahl = 9999;
  end;

  %
  switch auswahl

  case 1
    %
    % ------------ Beispiel 1 -----------
    %
    close all;
    disp('------------------ Beispiel 1 ------------------');
    disp(' ');
    disp('- BPSK');
    disp('- zufaelliger Kanal 2. Ordnung');
    disp('- 8 gesendete Symbole (inkl. Terminierung)');
    disp(' ');
    %
    %
    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');


    %
    h_tmp   = randn(3,1) + j*randn(3,1);    % Kanal auswuerfeln
    h       = h_tmp / sqrt(sum(abs(h_tmp).^2)); % Kanal normieren
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    all_sym = [-1 1].';                     % Symbolalphabet
    subplot(223);
    plot(all_sym+(j*1e-10),'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    axis('square')
    grid;
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 12-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal

    % ---- Rauschen
    gamma   = 1/sqrt(2) * 10^(-snr/20);
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));
    %
    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 5;
    grafik_flag = 1;
    term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %
    %
  case 2
    %
    % ------------ Beispiel 2 -----------
    %
    close all;
    disp('------------------ Beispiel 2 ------------------');
    disp(' ');
    disp('- BPSK');
    disp('- zufaelliger Kanal 3. Ordnung');
    disp('- 16 gesendete Symbole (inkl. Terminierung)');
    disp('- grafische Darstellung ohne PAUSE');
    disp(' ');
    %
    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');
    %
    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');

    %
    %
    h_tmp   = randn(4,1) + j*randn(4,1);    % Kanal auswuerfeln
    h       = h_tmp / sqrt(sum(abs(h_tmp).^2)); % Kanal normieren
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    all_sym = [-1 1].';                     % Symbolalphabet
    subplot(223);
    plot(all_sym+(j*1e-10),'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 16-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
    %
    % ---- Rauschen
    gamma   = 1/sqrt(2) * 10^(-snr/20);
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));

    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 10;
    grafik_flag = 1;   % grafische Darstellung ohne PAUSE
    term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %
    %
    %
  case 3
    %
    % ------------ Beispiel 3 -----------
    %
    close all;
    disp('------------------ Beispiel 3 ------------------');
    disp(' ');
    disp('- QPSK');
    disp('- zufaelliger Kanal 2. Ordnung');
    disp('- 8 gesendete Symbole (inkl. Terminierung)');
    disp(' ');
    %
    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');

    %
    %
    h_tmp   = randn(3,1) + j*randn(3,1);    % Kanal auswuerfeln
    h       = h_tmp / sqrt(sum(abs(h_tmp).^2)); % Kanal normieren
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    all_sym = [(-1 -j) (1 -j) (1 +j) (-1 +j)].'/sqrt(2);  % Symbolalphabet
    subplot(223);
    plot(all_sym+(j*1e-10),'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 8-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
    %
    % ---- Rauschen
    gamma   = 1/sqrt(2)*10^(-snr/20)/sqrt(2); % incl Eb --> Es-Korrektur
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));

    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 10;
    grafik_flag = 1;
    term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %
    %
    %
    %
  case 4
    %
    % ------------ Beispiel 4 -----------
    %
    close all;
    disp('------------------ Beispiel 4 ------------------');
    disp(' ');
    disp('- QPSK');
    disp('- zufaelliger Kanal 3. Ordnung');
    disp('- 16 gesendete Symbole (inkl. Terminierung)');
    disp('- grafische Darstellung ohne PAUSE');
    disp(' ');
    %
    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');

    %
    %
    h_tmp   = randn(4,1) + j*randn(4,1);    % Kanal auswuerfeln
    h       = h_tmp / sqrt(sum(abs(h_tmp).^2)); % Kanal normieren
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    all_sym = [(-1 -j) (1 -j) (1 +j) (-1 +j)].'/sqrt(2);  % Symbolalphabet
    subplot(223);
    plot(all_sym+(j*1e-10),'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 16-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
    %
    % ---- Rauschen
    gamma   = 1/sqrt(2)*10^(-snr/20)/sqrt(2); % incl Eb --> Es-Korrektur
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));

    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 10;
    grafik_flag = 2;
    term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %
  case 5
    %
    % ------------ Beispiel 5 -----------
    %
    close all;
    disp('------------------ Beispiel 5 ------------------');
    disp(' ');
    disp('- BPSK');
    disp('- 8 gesendete Symbole (inkl. Terminierung)');
    disp(' ');
    %
    %
    disp('Bitte Kanalimpulsantwort eingeben')
    disp(' ');
    h0_txt = input('h0= ','s');
    h0 = str2num(h0_txt);
    disp(' ')
    h1_txt = input('h1= ','s');
    h1 = str2num(h1_txt);
    disp(' ');
    h2_txt = input('h2= ','s');
    h2 = str2num(h2_txt);
    disp(' ')

    h=[h0 h1 h2];
    h       = h / sqrt(sum(abs(h).^2)); % Kanal normieren


    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');
    %

    disp(' ');
    term_txt = input('Trellis-Terminierung: ja (1) -- nein (0) ','s');
    term_flag = str2num(term_txt);

    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    all_sym = [-1 1].';                     % Symbolalphabet
    subplot(223);
    plot(all_sym+(j*1e-10),'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 12-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal

    % ---- Rauschen
    gamma   = 1/sqrt(2) * 10^(-snr/20);
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));
    %
    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 5;
    grafik_flag = 1;
    %term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %
    %
  case 6
    %
    % ------------ Beispiel 6 -----------
    %
    close all;
    disp('------------------ Beispiel 6 ------------------');
    disp(' ');
    disp('- QPSK');
    disp('- zufaelliger Kanal 2. Ordnung');
    disp('- 8 gesendete Symbole (inkl. Terminierung)');
    disp(' ');
    %
    disp('Bitte Kanalimpulsantwort eingeben')
    disp(' ');
    h0_txt = input('h0= ','s');
    h0 = str2num(h0_txt);
    disp(' ')
    h1_txt = input('h1= ','s');
    h1 = str2num(h1_txt);
    disp(' ');
    h2_txt = input('h2= ','s');
    h2 = str2num(h2_txt);
    disp(' ')

    h=[h0 h1 h2];
    h       = h / sqrt(sum(abs(h).^2)); % Kanal normieren

    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');


    disp(' ');
    term_txt = input('Trellis-Terminierung: ja (1) -- nein (0) ','s');
    term_flag = str2num(term_txt);

    %
    %
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    all_sym = [(-1 -j) (1 -j) (1 +j) (-1 +j)].'/sqrt(2);  % Symbolalphabet
    subplot(223);
    plot(all_sym+(j*1e-10),'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 8-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
    %
    % ---- Rauschen
    gamma   = 1/sqrt(2)*10^(-snr/20)/sqrt(2); % incl Eb --> Es-Korrektur
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));

    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 10;
    grafik_flag = 1;
    %term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %



    %
  case 7
    %
    % ------------ Beispiel 7 -----------
    %
    close all;
    disp('------------------ Beispiel 7 ------------------');
    disp(' ');
    disp('- 16-QAM');
    disp('- zufaelliger Kanal 1. Ordnung');
    disp('- 8 gesendete Symbole (inkl. Terminierung)');
    disp(' ');
    %
    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');

    %
    %
    h_tmp   = randn(2,1) + j*randn(2,1);    % Kanal auswuerfeln
    h       = h_tmp / sqrt(sum(abs(h_tmp).^2)); % Kanal normieren
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    sym_tmp = repmat([-3 -1 1 3],4,1) + j * repmat([-3 -1 1 3].',1,4);
    all_sym = sym_tmp(:) / sqrt(10); % Symbolalphabet
    subplot(223);
    plot(all_sym,'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 8-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
    %
    % ---- Rauschen
    gamma  = 1/sqrt(2)*sqrt(2)*10^(-snr/20)/2; % incl Eb -> Es-Korrektur ld(M)=4
    r_sig  = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));

    %
    s_sym  = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig  = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 10;
    grafik_flag = 1;
    term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    s_ref=s_sym(1:length(r_sym));
    P_s = mean(sign(abs(s_ref-r_sym)));
    disp(sprintf('Symbolfehlerrate: %7.2e',P_s));
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %
  case 8
    %
    % ------------ Beispiel 8 -----------
    %
    close all;
    disp('------------------ Beispiel 8 ------------------');
    disp(' ');
    disp('- 16-QAM');
    disp('- zufaelliger Kanal 2. Ordnung');
    disp('- 16 gesendete Symbole (inkl. Terminierung)');
    disp('- grafische Darstellung ohne PAUSE');
    disp(' ');
    %
    %
    disp(' ');
    snr_txt = input('Bitte Eb/N0-Verhaeltnis in dB angeben --> ','s');
    snr = str2num(snr_txt);
    disp(' ');

    disp(' ');
    iv_txt = input(['Entscheidungsverzoegerung iv eingeben (gelb ',...
                    'markierter Pfad)--> '],'s');
    iv = str2num(iv_txt);
    disp(' ');
    disp(' ');

    disp(' ');
    %
    h_tmp   = randn(3,1) + j*randn(3,1);    % Kanal auswuerfeln
    h       = h_tmp / sqrt(sum(abs(h_tmp).^2)); % Kanal normieren
    %
    figure(1);
    subplot(211);
    stem([0:length(h)-1],abs(h));
    axis([-0.2 length(h)-0.8 0 1.2*max(abs(h))]);
    xlabel('i ->');
    ylabel('|h(i)|');
    title('Kanalimpulsantwort');
    grid;
    %
    %
    sym_tmp = repmat([-3 -1 1 3],4,1) + j * repmat([-3 -1 1 3].',1,4);
    all_sym = sym_tmp(:) / sqrt(10); % Symbolalphabet
    subplot(223);
    plot(all_sym,'o');
    axis([-1.5 1.5 -1.5 1.5]);
    xlabel('Re\{d\} ->');
    ylabel('Im\{d\} ->');
    title('Symbolalphabet');
    grid;
    axis('square')
    %
    %
    % ------------ Auswuerfeln der Symbole -----------
    %
    % wichtig: Einschwingen des Kanals mit Symbol Nr. 0
    %          Terminierung der Sequenz mit Symbol Nr. 0
    %
    n_data  = 16-length(h)+1;
    sym_num = [zeros(length(h)-1,1);...
               floor(length(all_sym) * rand(n_data,1));...
               zeros(length(h)-1,1)];
    s_sig   = all_sym(sym_num+1);    % Symbolfolge (Sendesignal)
    %
    r_sig   = filter(h,1,s_sig);     % Faltung mit Kanal
    %
    % ---- Rauschen
    gamma   = 1/sqrt(2)*10^(-snr/20)/2; % incl Eb --> Es-Korrektur
    r_sig   = r_sig + gamma*(randn(length(r_sig),1) + j*randn(length(r_sig),1));

    %
    s_sym   = s_sig(length(h):length(s_sig)); % Einschwingvorgang abschneiden
    r_sig   = r_sig(length(h):length(r_sig)); % Einschwingvorgang abschneiden
    %
    %
    subplot(224);
    plot(abs(s_sym),'g');
    hold on;
    plot(abs(r_sig),'r');
    hold off;
    xlabel('i ->');
    title('|Signale|');
    legend('vor Kanal','nach Kanal');
    grid;

    %
    % ------------ Viterbi-Aufruf -------------
    figure(2);
    %iv          = 10;
    grafik_flag = 2;
    term_flag   = 1;
    r_sym       = viterbi_entzerrer(r_sig,h,all_sym,term_flag,...
                                    grafik_flag,iv,s_sym);
    %
    %
    disp(' ');
    disp(' ');
    disp(' ');
    disp('zurueck zur Auswahl mit bel. Taste');
    pause
  end;

end;

% ### EOF ######################################################################
