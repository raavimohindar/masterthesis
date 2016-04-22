% ##############################################################################
% ## kkf_aper.m : Berechnung der {aper}iodischen Kreuzkorrelationsfunktion    ##
% ##############################################################################
%
% Aufruf kkf_aperiodic = kkf_aper(x1,x2)
%
% Eingabe:
%
%   x1,x2 :  Spaltenvektoren mit den zu korrelierenden Sequenzen 
%
% Ausgabe:
%
% kkf_aperdiodic:  aperiodische KKF 
%
% Volker Kuehn, 05.07.01

function kkf_aperiodic = kkf_aper(x1,x2)

x1_len = length(x1);
x2_len = length(x2);
if (x1_len ~= x2_len)
  disp('kkf_aper.m: Sequenzen haben nicht die gleiche Laenge');
  break;
end

kkf_aperiodic = zeros(2*x1_len-1,1);

x1 = x1(:).';
x2 = x2(:);
for i=-x1_len+1:0
  kkf_aperiodic(i+x1_len) = x1(1-i:x1_len) * x2(1:x2_len+i);
end
for i=1:x1_len-1
  kkf_aperiodic(i+x1_len) = x1(1:x1_len-i) * x2(1+i:x2_len);
end

% ### EOF ######################################################################
