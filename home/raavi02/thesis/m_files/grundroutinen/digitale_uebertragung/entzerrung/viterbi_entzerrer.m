% ##############################################################################
% ##  viterbi_entzerrer.m : Viterbi-Entzerrung                                ##
% ##############################################################################
%
% Aufruf:   out_sym =
%           viterbi_entzerrer(in_sym,h,all_sym,term_flag,grafik_flag,i0,source);
%
% Eingabe:  in_sym         : empfangene Sequenz
%           h              : Kanalimpulsantwort
% optional:
%           all_sym        : Vektor mit allen moeglichen Symbolen 
%                            (default = [-1,1])
%                            Anstelle eines Vektors kann auch eines der
%                            voreingestellten Alphabete durch Angabe
%                            der folgenden Nummern verwendet werden:
%                            1 : BPSK (-1 / +1)
%                            2 : QPSK (-1-j,1-j,1+j,-1+j)/sqrt(2)
%                            3 : 8-PSK (spiegelsymmetrisch)
%                            4 : 16-QAM
%            term_flag     : Dieser Schalter gibt an, ob die Sequenz terminiert
%                            war. Zur Terminierung muss die Sendesequenz mit dem
%                            ersten in 'all_sym' enthaltenen Symbol terminiert
%                            werden.
%                            (default=1)
%            grafik_flag   : Aktivierung der grafischen Trellisdarstellung
%                            (default=0). Angezeigt werden alle aktuellen Pfade
%                            und Pfadkosten. Am Ende wird der wahre Pfad gelb 
%                             eingefaerbt.
%                            'grafik_flag=1' nach jedem Takt muss eine beliebige
%                             Taste gedrueckt werden
%                            'grafik_flag=2' Abschaltung der PAUSE-Funktion
%            i0            : Einfaerbung (gelb) der i0 Takte zurueckliegenden 
%                            Pfade
%                            (default=0)
%            source        : Gesendete Symbole zur Darstellung des wahren Pfades
%                            (default=[])
%                            Der wahre Pfad wird rot eingefaerbt.
%
% Ausgabe:   out_sym       : entschiedene (entzerrte) Symbole
%
%                                       Heiko Schmidt 2001(University of Bremen)

function out_sym = viterbi_entzerrer(in_sym,h,all_sym,term_flag,grafik_flag,...
                                     i0,source);

if nargin<7            % pruefen, ob source-Vektor uebergeben wurde
  source = [];
  if nargin <6         % pruefen, ob i0 uebergeben wurde
    i0 = 0;
    if nargin <5       % pruefen ob GRAPHIC_FLAG gesetzt wurde
      grafik_flag = 0;
      if nargin < 4    % pruefen ob TERM_FLAG gesetzt wurde
        term_flag = 1;
        if nargin < 3  % pruefen, ob Symbolvektor uebergeben wurde
          all_sym = [-1 1];
        else
          if length(all_sym) == 1
            switch all_sym
            case 1
              all_sym = [-1 ; 1];
            case 2
              all_sym = [(-1 -j);(1 -j);(1 +j);(-1 +j)]/sqrt(2);
            case 3
              all_sym = exp(pi*(0.125+0.25*[0:7]'));
            case 4
              tmp = repmat([-3 -1 1 3],4,1) + j * repmat([-3 -1 1 3].',1,4);
              all_sym = tmp(:) / sqrt(10);
            end;
          end;
        end;
      end;
    end;
  end;
end;
%
%
% ---------------- Allgemeine Definitionen -------------------
%
num_data        = length(in_sym);   % Laenge des empfangenen Datenvektors
M               = length(all_sym);  % Anzahl der moeglichen Symbole
h_len           = length(h);        % Laenge der Kanalimpulsantwort
m               = h_len - 1;        % Kanalordnung
num_state       = M^m;              % Anzahl der Zustaende im Trellis
inf_met         = 10e16;            % Startpfadkosten fuer S1 - ...
%
num_state_msb   = M^(m-1);          % MSB bei Zustandsdarstellung
state_mat       = repmat([0:num_state-1]',1,M);  % temporaer benutzte 
%                                                  Zustandsmatrix
cur_sym_mat     = repmat([0:M-1],num_state,1);   % temporaer benutzte
%                                                  Symbolmatrix
sym_vec         = floor([0:num_state-1].'/num_state_msb); % Vektor mit zu jedem
%                                                        Zustand gehoerenden
%                                                        Eingangssymbolnummern
dec_sym_mat     = repmat(sym_vec,1,M);                    % Matrix zu 'sym_vec'
sym_out         = zeros(num_state,num_data);   % Zu jedem Zeitpunkt in jedem
%                                                Zustand entschiedene
%                                                Symbolnummern
trellis         = zeros(num_state,num_data);   % Entschiedene Trellispfade / 
%                                                zu jedem Zeitpunkt in jedem
%                                                Zustand --> woher
metric_vec      = inf_met * ones(num_state,1); % Aktuelle Pfadkosten aller 
%                                                Zustaende auf sehr grossen Wert
%                                                setzten
metric_vec(1)   = 0;                           % Initialisierung des Zustande S0
% ------------------------------------------------------------------------------
pre_state     = mod(state_mat * M,num_state) + (cur_sym_mat); 
% Diese Matrix enthaelt das Trellisdiagramm:
%  Aus welchem Zustand kommt man, wenn das Symbol (Spaltennummer) gesendet wurde
% ------------------------------------------------------------------------------
channel_input = double((dec2base(state_mat+(cur_sym_mat*num_state),M))-48);
channel_input = channel_input - (7*(channel_input >= M));
h_mat         = repmat(h(:).',M^(m+1),1);
z_vec         = sum((all_sym(channel_input+1) .* h_mat).');
z_mat         = zeros(num_state,M);
% ------------------------------------------------------------------------------
z_mat(:)      = z_vec; % Diese Matrix enthaelt die zu jedem Zustand 
%                        und Pfad gehoerenden Datenniveaus
% ------------------------------------------------------------------------------
%
%
if length(source) % Berechnung der Trellisstruktur zur Darstellung 
  %                 des wahren Pfades
  next_state = floor((num_state*cur_sym_mat + state_mat) /M);
  source_sym_num = zeros(length(source),1);
  for slauf = 1:length(source)
    [min_val,min_index]   = min(abs(all_sym-source(slauf)));
    source_sym_num(slauf) = min_index-1;
  end;
end;
% ------------------------------------------------------------------------------
% ---------------- Initialisierung der grafischen Ausgabe ----------------------
%
if grafik_flag
  axis off;                      % Achsen aus
  delta_x = 1/num_data;          % Skalierungsfaktor in x-Richtung
  delta_y = -1/(num_state-1);    % Skalierungsfaktor in y-Richtung
  hold on;
  for x_lauf = 0:num_data
    for y_lauf = 0:(num_state-1)
      plot(x_lauf*delta_x,y_lauf*delta_y,'o');    % Knoten zeichnen
      if x_lauf == 0
        text(-delta_x,y_lauf*delta_y,sprintf('S_{%d}',y_lauf)); 
        % Zustandsnamen darstellen
      end;

    end;
  end;

end;
% ---------------------- Ende der Initialisierung ------------------------------
%
% ------------------------------------------------------------------------------
% --------------------- Trellisdurchlauf (vorwaerts) ---------------------------
for lauf = 1:num_data
  metric               = repmat(metric_vec(:),1,M);                 
  % ver-M-fachung des Pfadkostenvektors
  metric               = metric + (abs(z_mat - in_sym(lauf)).^2);   
  % Berechnung der zu jedem Pfad gehoerenden Summenpfadkosten
  auswahl            = metric(pre_state+(dec_sym_mat*num_state)+1); 
  % Sortieren der Pfadkosten
  [path_value,min_pos] = min(auswahl.');                            
  % Entscheidung fuer geringste Summenpfadkosten
  index                = ([0:num_state-1] + (min_pos-1)*num_state); 
  % Adresse faer Pfad mit geringsten Pfadkosten
  trellis(:,lauf)      = pre_state(index+1).';                      
  % entschiedene Pfade eintragen
  sym_out(:,lauf)      = sym_vec;                                   
  % entschiedenen Symbole eintragen
  % ----------------------------------------------------------------------------
  % --------------------------- Nur fuer graphische Ausgabe --------------------
  %
  %
  %
  if grafik_flag
    pl = [];
    text_height = min(10,340/(num_state*M));    % Texthaehe bestimmen (anpassen)
    if term_flag
      plot_num_state = M^(min(m,num_data - lauf)); 
      % Anzahl der darzustellenden Zustaende (fuer Terminierung)
    else
      plot_num_state = M^m;
    end;
    if grafik_flag < 2
      for l_lauf = 0:plot_num_state-1
        plot_M = M;
        metric_txt = '';
        for s_lauf   = 0:plot_M-1
          y0=pre_state(l_lauf+1,s_lauf+1);
          pl(y0+1,l_lauf+1) = plot([(lauf-1)*delta_x lauf*delta_x],...
                                   [y0 l_lauf]*delta_y);
          if auswahl(l_lauf+1,s_lauf+1) < inf_met
            metric_txt  = sprintf('%s%6.1f',metric_txt,...
                                  auswahl(l_lauf+1,s_lauf+1));
          else
            metric_txt  = sprintf('%s   +\\infty',metric_txt);
          end;
          if s_lauf < plot_M-1
            metric_txt  = sprintf('%s\n',metric_txt);
          end;
        end;
        % Ausgabe der aktuellen Summenpfadkosten 
        % (vor der jeweiligen Entscheidung)
        if text_height > 2
          t1=text((lauf*delta_x+0.01),(l_lauf*delta_y),metric_txt);
          set(t1,'fontsize',text_height);
        end;
      end;
      if length(source)
        begin_state = 0;
        for source_lauf = 0:lauf-1
          new_state = next_state(begin_state+1,source_sym_num(source_lauf+1)+1);
          plot([source_lauf source_lauf+1]*delta_x,...
               [begin_state new_state]*delta_y,'--r');
          begin_state = new_state;
        end;
      end;
      pause
    end;

    if grafik_flag < 2 | lauf==num_data

      hold off;
      for x_lauf = 0:num_data
        for y_lauf = 0:(num_state-1)
          plot(x_lauf*delta_x,y_lauf*delta_y,'o');
          hold on;
          if x_lauf == 0
            text(-delta_x,y_lauf*delta_y,sprintf('S_{%d}',y_lauf));
          end;

        end;
      end;
      axis off;

      if lauf < num_data
        for l_lauf = 0:plot_num_state-1
          best_back = l_lauf+1;
          for ruecklauf = lauf:-1:1
            best_back_old     = best_back;
            best_back         = trellis(best_back,ruecklauf)+1;
            %
            pl=plot([(ruecklauf-1)*delta_x ruecklauf*delta_x],...
                    [(best_back-1)*delta_y (best_back_old-1)*delta_y]);
            if i0
              if (lauf-ruecklauf) >  (i0-2)
                set(pl,'linewidth',2);
                set(pl,'color','y');
              end;
            end;
          end;
        end;
      end;
    end;
  end;
  %
  %
  % ------------------------ Ende graphische Ausgabe --------------------------
  %
  %
  metric_vec = auswahl(index+1);   
  % Erzeugung des aktualisierten Metrikvektors nach der Pfadentscheidung
  %
end;
% -----------------------------
%
%
%  Darstellung aller entschiedenen Pfade, falls keine Terminierung erfolgte !!!
%
if (grafik_flag & (term_flag == 0))
  text_height = min(10,200/(num_state));    % Texthaehe bestimmen (anpassen)

  for best_back = 1:num_state

    metric_txt  = sprintf('%6.1f',metric_vec(best_back));
    if text_height > 2
      t1=text((lauf*delta_x+0.01),((best_back-1)*delta_y),metric_txt);
      set(t1,'fontsize',text_height);
      if metric_vec(best_back) == min(metric_vec)
        set(t1,'color','y');
      end;
    end;

    %   cur_txt =

    for ruecklauf = num_data:-1:1
      out_sym(ruecklauf) = all_sym(sym_out(best_back,ruecklauf)+1);   
      % entschiedenes Symbol eintragen
      best_back_old     = best_back;                                  
      % bisherigen Zustand sichern
      best_back         = trellis(best_back,ruecklauf)+1;             
      % vorhergehenden Zustand ermitteln
      %
      if grafik_flag  
        % nur fuer grafische Ausgabe: Eintragung des entschiedenen Pfades 
        pr=plot([(ruecklauf-1)*delta_x ruecklauf*delta_x],...
                [(best_back-1)*delta_y (best_back_old-1)*delta_y],'b');
        %      set(pr,'linewidth',2)
      end;            
    end;
  end;

end;
%

%
%
% ----------  Terminierung  -----------
if term_flag
  best_back = 1;  % Startzustand S_0 fuer Pfadrueckverfolgung bei Terminierung
else
  [min_metric,best_back] = min(metric_vec); % bester Startzustand
end;
% -------------------------------------
out_sym        = zeros(num_data,1);   % Initialisierung des Rueckgabevektors 
%                                       (entschiedene Daten)
% -------------------------------------
%
%
%
%
% -------------- Pfadrueckverfolgung ----------------
for ruecklauf = num_data:-1:1
  out_sym(ruecklauf) = all_sym(sym_out(best_back,ruecklauf)+1);   
  % entschiedenes Symbol eintragen
  best_back_old     = best_back;                                  
  % bisherigen Zustand sichern
  best_back         = trellis(best_back,ruecklauf)+1; 
  % vorhergehenden Zustand ermitteln
  %
  if (grafik_flag&((num_data-ruecklauf)>=i0))  
    % ------- nur fuer grafische Ausgabe: Eintragung des entschiedenen Pfades 
    pr=plot([(ruecklauf-1)*delta_x ruecklauf*delta_x],...
            [(best_back-1)*delta_y (best_back_old-1)*delta_y],'y');
    set(pr,'linewidth',2)
  end;            
end;
% ------------------------ Ende der Pfadrueckverfolgung ----------------------
%
if (length(source)&grafik_flag)        
  % ------------------------ Anzeige des wahren Pfades (optional) -------------
  begin_state = 0;
  for source_lauf = 0:num_data-1
    new_state = next_state(begin_state+1,source_sym_num(source_lauf+1)+1);
    plot([source_lauf source_lauf+1]*delta_x,...
         [begin_state new_state]*delta_y,'--r');
    begin_state = new_state;
  end;
end;                 
% ---------------------------------------------------------------------------
%
if ((term_flag == 0) & (i0 > 0))
  out_sym = out_sym(1:length(out_sym)-i0);
end;
% ------------------------------------
% Ende der Funktion viterbi_entzerrer
% ------------------------------------

% ### EOF ######################################################################
