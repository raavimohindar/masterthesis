% ##############################################################################
% ## kkf_per.m : Berechnung der {per}iodischen Kreuzkorrelationsfunktion      ##
% ##############################################################################
%
% Aufruf kkf_periodic = kkf_per(x1,x2)
%
% Eingabe:
%
%   x1,x2 :  Spaltenvektoren mit den zu korrelierenden Sequenzen 
%
% Ausgabe:
%
% kkf_perdiodic:  periodische KKF 
%
% Volker Kuehn, 05.07.01

function kkf_periodic = kkf_per(x1,x2)

x1_len = length(x1);
x2_len = length(x2);
if (x1_len < x2_len)
  x  = x1;
  x1 = x2;
  x2 = x; 
  clear x;
  x1_len = length(x1);
  x2_len = length(x2);
end

kkf_periodic = zeros(x1_len+1,1);

x1 = x1(:).';
x2 = x2(:);
for i=0:x1_len
  kkf_periodic(i+1) = x1(1:x2_len) * x2;
  x1 = [x1(x1_len) x1(1:x1_len-1)];
end

% ### EOF ######################################################################
