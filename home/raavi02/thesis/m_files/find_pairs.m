% ##############################################################################
% ## find_pairs.m : Auffinden der bevorzugten Paare von m-Sequenzen           ##
% ##############################################################################
%
% Aufruf pairs = find_pairs(fb_polys)
%
% Eingabe:
%
%   fb_polys :  Matrix mit den Feedback Polynomen
%               jede Spalte enthaelt Polynom in GF(2)
%
% Ausgabe:
%
%   pairs:  quadratische Matrix mit Einsen an den Stellen mit den bevorz. Paaren
%
% Bemerkungen:
%  benoetigt: m_seq.m, kkf_per.m
%
% Volker Kuehn, 05.07.01

function pairs = find_pairs(fb_polys)

[m,N_seq] = size(fb_polys);   % Schieberegister Laenge+1/Anzahl der Polynome 

m = m - 1;

period = 2^m-1;

seqs  = zeros(period,N_seq);
start = [zeros(m-1,1);1];

for j=1:N_seq
  seqs(:,j) = m_seq(fb_polys(:,j),start);
end

tmp   = 2^(floor((m+2)/2));
terms = [-1; -tmp-1; tmp-1];       % drei Korrelationswerte der bev. Paare

pairs = zeros(N_seq,N_seq);
cntr  = 1;

for i=1:N_seq-1
  for j=i+1:N_seq
    kkf = kkf_per(seqs(:,i),seqs(:,j));
    set = unique(round(kkf*period));
    if isempty(setdiff(set,terms))
      pairs(i,j) = 1;
      pairs(j,i) = 1;
    end
    cntr = cntr + 1;
  end
end

% ### EOF ######################################################################
