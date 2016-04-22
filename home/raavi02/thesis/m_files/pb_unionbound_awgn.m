% ##############################################################################
% ##  pb_unionbound_awgn.m : Abschaetzung der Bitfehlerwahrscheinlichkeit     ##
% ##                         per Union Bound                                  ##
% ##############################################################################
%
% function [Pb]=pb_unionbound_awgn(EbN0db,ACwh,n,k)
%
% EINGABE:
%   EbN0db: Signal/Noise pro Infobit [db]
%           (Skalar, Spalten- oder Zeilenvektor)
%     ACwh: Input Output Weight Enumerating Function
%           w=Infowortgewicht (Matrix-Index 1 entspricht w=0 !!!)
%           h=Codewortgewicht (Matrix-Index 1 entspricgt h=0 !!!)
%           (Matrix)
%        n: Laenge der Codeworte
%           (Skalar)
%        k: Laenge der Infoworte
%           (Skalar)
%
% AUSGABE:
%   Pb: Bitfehlerwahrscheinlichkeit
%       (Spaltenvektor)
%
% ANMERKUNGEN:
%   - Ausgabe stellt eine obere Grenze der Bitfehlerwahrscheinlichkeit dar
%   - schlechte Naeherung fuer kleine EbN0db
%
% QUELLE:
%   [Fri96,S.286,(9.5.3)]
%
% AUTOR: Juergen Rinas,  31.05.1999
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pb]=pb_unionbound_awgn(EbN0db,ACwh,n,k)

EbN0=10.^(EbN0db/10);

Rc=k/n;

cd=([0:(size(ACwh,1)-1)] * ACwh);
d=(0:(length(cd)-1));

Pb=zeros(length(EbN0db),1);

for EbN0dbidx=1:length(EbN0db)
  Pb(EbN0dbidx)=0.5 * cd(d+1)/k * erfc(sqrt( d.'*Rc*EbN0(EbN0dbidx)));
end;

% ### EOF ######################################################################
