% ##############################################################################
% ##  demultiplex.m : Demultiplexer einer codierten Sequenz,RSC Codes         ##
% ##                  der Rate 1/n                                            ##
% ##############################################################################
%
% function: dem_out = demultiplex(r, codes, puncture);
% ------------------------------------------------------------------------------
% EINGABE:
%     r:         Empfangssequenz (punktiert)
%        codes     : Struktur mit den Elementen code_1 ... code_num_codes,
%           num_codes   : Anzahl der Codes
%           code_i      : Struktur des Codes i mit folgenden Elementen
%           trellis     : Struktur des Trellis
%             out       : Codewoerter fuer jeden Zweig
%             next      : Folgezustand im Trellis
%           num_state   : Anzahl der Trellis Zustaende des Codes
%           block_len   : Anzahl der Codewoerter pro Rahmen
%           word_len    : Anzahl der Bit pro Codewort
%           term        : term==1 fuer terminierten Trellis, ansonsten
%                         nicht terminiert
%           Puncture    : Vektor mit der punktierten Ausgangssequenz im
%                         codierten Bitstream
%                         (enthaelt die Positionen der uebertragenen Symbole im
%                          nichtpunktierten Strom)
%           IL          : Vektor mit der Interleaversequenz
%           out         : out == 0, nur LLRs der Infobit werden berechnet
%                         out == 1, LLRs der Infobit und Codebit werden
%                                   berechnet
%
% AUSGABE:
%        dem_out        : Matrix der Groesse (N_total*max_n,codes.num_codes)
%                         enthaelt in jeder Spalte die Codesequenz des
%                         zugehoerigen Codes
% ------------------------------------------------------------------------------
%
function dem_out = demultiplex(r, codes, puncture);


N_total = codes.code_1.block_len;

n_i = zeros(codes.num_codes,1);
for i=1:codes.num_codes
  eval(['n_i(i) = codes.code_' int2str(i) '.word_len;']);
end
n     = sum(n_i) - codes.num_codes + 1;
max_n = max(n_i);

dem_out = zeros(N_total*max_n,codes.num_codes);
dummy   = zeros(N_total*n,1);
dummy(puncture) = r;               % Einfuegen von Dummies wegen der Punktierung

select = 1:n:(N_total-1)*n+1;
y(:,1) = dummy(select);            % Infobit mit Tailbit

for i = 1:codes.num_codes
  for j=2:n_i(i)
    select = select + 1;
    y = [y dummy(select)];         % Parity Bit von Encoder i und Generator j
  end
end


% Umordnen der Kanaldaten, sodass jeder Decodierer seine Bit erhaelt
start = 2;
for i = 1:codes.num_codes
  eval(['IL=codes.code_' int2str(i) '.IL;']);
  dem_out(:,i) = reshape([y(IL,1) y(:,start:start+n_i(i)-2)]',...
  N_total*n_i(i), 1);
  start = start + n_i(i)-1;
end

% ### EOF ######################################################################
