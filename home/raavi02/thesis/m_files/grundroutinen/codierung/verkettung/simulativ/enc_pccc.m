% ##############################################################################
% ##  enc_pccc.m : Turbo-Encoder fuer n parallel verkettete RSC Codes         ##
% ##############################################################################
%
% function: [codebits,last_states] = enc_pccc( x, N_input, codes, P)
% ------------------------------------------------------------------------------
% EINGABE:
%   x              : Vektor mit Infobit
%   N_input        : Laenge des Eingangsvektors x
%   codes          : Struktur mit den Elementen code_1 ... code_num_codes,
%     code_i       : Struktur des Codes i mit folgenden Elementen
%     g            : Generatorpolynome des Codes
%     trellis      : Struktur des Trellis
%        out       : Codewoerter fuer jeden Zweig
%        next      : Folgezustand im Trellis
%     num_state    : Anzahl der Trellis Zustaende des Codes
%     block_len    : Anzahl der Codewoerter pro Rahmen
%     word_len     : Anzahl der Bit pro Codewort
%     term         : term==1 fuer terminierten Trellis, ansonsten
%                            nicht terminiert
%     P            : Vektor mit der punktierten Ausgangssequenz im
%                    codierten Bitstream
%                    die codierte Ausgangssequenz hat eine Laenge von N_input*3
%                    Wenn P = [1...1], dann unpunktiert, ergibt eine Rate 1/n
%     IL           : Vektor mit der Interleaversequenz
%
% AUSGABE:
%     codebits     : Codeworte mit Elemente (+1/-1)
%     last_states  : Vektor mit den letzten Zustaenden im Trellis (optional)
%
% ------------------------------------------------------------------------------
function [codebits,last_states] = enc_pccc( x, N_input, codes, P )

if (nargout==2)
  last_states = zeros(codes.num_codes,1);
end

N_total = N_input + log2(codes.code_1.num_state);   % Tailbit fuer
%                                                     terminierten Code
n_i  = zeros(codes.num_codes,1);
for i=1:codes.num_codes
  eval(['n_i(i) = codes.code_' int2str(i) '.word_len;']);
end
n      = sum(n_i) - codes.num_codes + 1;

% Generierung des Codewortes, das zum ersten RSC Encoder gehoert
% (Kein Interleaven)
% end = 1, perfekt terminiert;
output = conv_encoder(x, codes.code_1.g, 1, 1); % Sequenz enthaelt
%                                               Infobit und Paritybit

% Erstellen einer Matrix mit erster Reihe Infosequenz
% naechste Reihe Pruefbit des ersten RSC Encoders usw.
y = reshape(output,n_i(1),N_total);            % Jede Spalte enthaelt Codeworte
% N_total Spalten


% Interleave fuer weitere Encoder
for dim = 2:codes.num_codes
  eval(['IL = codes.code_' int2str(dim) '.IL;']);
  eval(['g  = codes.code_' int2str(dim) '.g;']);
  input  = y(1,IL);
  [output,ls,tail] = conv_encoder(input, g, 1, -1 );
  output = reshape(output,n_i(dim),N_total);                % Umordnen: s.o.
  y      = [y; output(2:n_i(dim),:)];
  if (nargout==2)
    last_states(dim) = ls;
  end
end


% Parallel zu Seriell Umwandlung in 1-Dimensionalen Vektor und Punktierung
temp = y(:);
codebits = temp(P');

% ### EOF ######################################################################
