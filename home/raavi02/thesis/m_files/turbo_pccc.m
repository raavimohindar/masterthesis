% ##############################################################################
% ##  turbo_pccc.m : Turbo-Decodierung von parallel verketteten RSC-Codes,    ##
% ##                 RSC-Codes der Rate 1/n, diverse Algorithmen,             ##
% ##                 erster Code terminiert und nicht interleaved             ##
% ##############################################################################
%
% function [L_info,L_ext] = turbo_pccc(signal,codes,P,num_it,alg)
% ------------------------------------------------------------------------------
% EINGABE:
%   signal        : Struktur mit folgenden Elementen
%      sig        : Empfangssignal
%      last_state : letzter Zustand im Trellis falls code_i.term==1
%                   (Spaltenvektor mit codes.num_codes Elementen)
%
%      codes     : Struktur mit den Elementen code_1 ... code_num_codes,
%         num_codes   : Anzahl der Codes
%         code_i      : Struktur des Codes i mit folgenden Elementen
%         trellis_out : Beschreibung der Trellisstruktur, Codewoerter fuer
%                       jeden Zweig
%         trellis_next: Folgezustand im Trellis
%         num_state   : Anzahl der Trellis Zustaende des Codes
%         block_len   : Anzahl der Codewoerter pro Rahmen
%         word_len    : Anzahl der Bit pro Codewort
%         term        : term==1 fuer terminierten Trellis, ansonsten
%                       nicht terminiert
%         P           : Vektor mit der punktierten Ausgangssequenz im codierten
%                       Bitstream
%                       (enthaelt die Positionen der uebertragenen Symbole im
%                        nichtpunktierten Strom)
%         IL          : Vektor mit der Interleaversequenz fuer Code i
%         num_it      : Anzahl der Iteration
%         alg         : String, der den Decodieralgorithmus festlegt
%                       'map','log_map','max_log_map'
%   codes             : struct containing the elements
%                       code_1 ... code_num_codes, (example for code_i)
% AUSGABE:
%   L_info         : Matrix mit LLRs der Infobit in Spalten pro Iteration
%   L_ext          : Matrix mit extrinsischer Information in den Spalten
%                    (optional)
%
%-------------------------------------------------------------------------------
% ANMERKUNGEN
%  benoetigt: demultiplex.m,
%             map, log_map, max_log_map Algorithmen
%
%-------------------------------------------------------------------------------
%% Copyright 2000 Volker Kuehn
%-------------------------------------------------------------------------------
function   [L_info,L_ext] = turbo_pccc(signal,codes,P,num_it,alg)

% Initialisierung
N_total = codes.code_1.block_len;
tail    = log2(codes.code_1.num_state);

n_i  = zeros(codes.num_codes,1);
for i=1:codes.num_codes
  eval(['n_i(i)    = codes.code_' int2str(i) '.word_len;']);
  select_syst(:,i) = [1:n_i(i):(N_total-1)*n_i(i)+1]';
  % Vektor mit den Positionen der syst. Bit
  % nur fuer Rc=1/n (Infobit sind getrennt durch n-1 Parity Bit)
end
n      = sum(n_i) - codes.num_codes + 1;
N_info = N_total - tail;


d_hat = zeros(N_info,num_it);

L_info = zeros(N_info,codes.num_codes*num_it);
if (nargout==2)
  L_ext = zeros(N_total,codes.num_codes*num_it);
end


% Am Empfaenger, Konversion vom seriellen zum parallen Strom
% um die Codeworte fuer jeden Encoder zu bekommen
data  = demultiplex(signal.sig,codes,P);


% Anfang der Decodierung

L_e   = zeros(N_total,codes.num_codes);


for it = 1:num_it
  for j = 1:codes.num_codes
    eval(['code=codes.code_' int2str(j) ';']);

    temp = sum(L_e,2) - L_e(:,j); % a priori Information generiert
    %                               durch andere Decoder

    in.sig        = data(:,j);
    in.L_a        = temp(code.IL);
    in.last_state = signal.last_state(j);

    eval(['L_xy = ' alg '(in,code);']);

    % Berechnung der extrinsischen Information
    temp = L_xy - data(select_syst(:,j),j) - in.L_a;
    L_e(code.IL,j) = temp;

    % De-Interleaven und Detektion der geschaetzen Daten pro Block
    temp(code.IL) = L_xy;
    L_info(:,j+codes.num_codes*(it-1))  = temp(1:N_info);

  end   % for j = 1:num_codes


  if (nargout==2)
    L_ext(:,(1:codes.num_codes)+codes.num_codes*(it-1)) = L_e;
  end

end    % for it = 1:iteration

% ### EOF ######################################################################
