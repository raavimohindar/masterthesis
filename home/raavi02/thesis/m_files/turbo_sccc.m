% ##############################################################################
% ##  turbo_sccc.m : Turbo-Decodierung von zwei seriell verketteten           ##
% ##                 Faltungscodes                                            ##
% ##                 Decodieralgorithmen koennen speziell angewaehlt werden   ##
% ##                 aeusserer Code ist terminiert, Eingangssequenz nicht     ##
% ##                 interleaved                                              ##
% ##############################################################################
%
% function  [L_info[,L_code_1[,L_code_2]]] = turbo_sccc(signal,codes,num_it,alg)
%
% ------------------------------------------------------------------------------
% EINGABE:
%        signal    :  Struktur mit folgenden Elementen
%           sig    :  Empfangssignal
%
%        codes     : Struktur mit den Elementen code_1 ... code_num_codes,
%           code_i      : Struktur des Codes i mit folgenden Elementen
%           trellis_out : Beschreibung der Trellisstruktur, Codewoerter fuer
%                         jeden Zweig
%           trellis_next: Folgezustand im Trellis
%           num_state   : Anzahl der Trellis Zustaende des Codes
%           block_len   : Anzahl der Codewoerter pro Rahmen
%           word_len    : Anzahl der Bit pro Codewort
%           term        : term==1 fuer terminierten Trellis, ansonsten
%                         nicht terminiert
%           P           : Vektor mit der punktierten Ausgangssequenz im
%                         codierten Bitstream
%                         (enthaelt die Positionen der uebertragenen Symbole
%                          im nichtpunktierten Strom)
%           IL          : Vektor mit der Interleaversequenz
%           num_it      : Anzahl der Iteration
%           alg         : String, der den Decodieralgorithmus festlegt
%                         'map','log_map','max_log_map'
% ------------------------------------------------------------------------------
% AUSGABE:
%        L_info         : Matrix, die in jeder Spalte die Info LLRs
%                         des Decodierdurchgangs enthaelt
%        L_code_1       : Matrix, die in jeder Spalte die Codeworte
%                         des aeusseren Codes enthaelt
%                         nach jedem Decodiervorgang (optional)
%        L_code_2       : Matrix, die in jeder Spalte die Codeworte
%                         des inneren Codes enthaelt
%                         nach jedem Decodiervorgang (optional)
%-------------------------------------------------------------------------------
% ANMERKUNGEN:
%   - benoetigt Datei: map
%                      log_map
%                      max_log_map
%-------------------------------------------------------------------------------
%% Copyright 2001 Volker Kuehn
%-------------------------------------------------------------------------------
function   [L_info,L_code_1,L_code_2] = turbo_sccc(signal,codes,num_it,alg)


%-------------------------------------------------------------------------------
% Initialisierung

in1.L_a        = zeros(codes.code_1.block_len,1);
in1.last_state = 0;

in2.L_a        = zeros(codes.code_2.block_len,1);
in2.last_state = 0;

% Einfuegen von Dummies an den punktierten Stellen fuer den inneren Code
in2.sig                 = zeros(codes.code_2.block_len*codes.code_2.word_len,1);
in2.sig(codes.code_2.P) = signal.sig;

in1.sig                 = zeros(codes.code_1.block_len*codes.code_1.word_len,1);

L_a = zeros(codes.code_2.block_len,1);


L_info = zeros(codes.code_1.block_len,1);
if (nargout==2)
  L_code_1 = zeros(codes.code_2.block_len,1);
end
if (nargout==3)
  L_code_2 = zeros(codes.code_2.block_len*codes.code_2.block_len,1);
end


% Iterative Decodierung


for it = 1:num_it

  in2.L_a = L_a;

  % Innerer Code
  if (nargout==3)
    eval(['[L_xy,L_code_2(:,it)] = ' alg '(in2,codes.code_2);']);
  else
    eval(['L_xy = ' alg '(in2,codes.code_2);']);
  end

  % Generierung und De-Interleaving der extrinsischen Information

  tmp(codes.IL) = L_xy - in2.L_a;       % Deinterleaver
  in1.sig(codes.code_1.P) = tmp(:);     % Einfuegen von Dummies an den
  % punktierten Stellen fuer den aeusseren Code


  % aeusserer Code

  eval(['[L_info(:,it),L_xy] = ' alg '(in1,codes.code_1);']);

  if (nargout==2)
    L_code_1(:,it) = L_xy;
  end

  % Generierung und De-Interleaving der extrinsischen Information

  L_e = L_xy - in1.sig(:);
  L_e = L_e(codes.code_1.P);
  L_a = L_e(codes.IL);         % Interleaven der extrinsischen Information

end    % for it = 1:num_it

% ### EOF ######################################################################
