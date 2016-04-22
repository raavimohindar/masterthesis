% ##############################################################################
% ##  pb_dpsk_koh_awgn: Bitfehlerwahrsch. fuer DPSK beim AWGN Kanal,          ##
% ##                    koh. Det.                                             ##
% ##############################################################################
%
% function [Ps,Pb] = pb_dpsk_koh_awgn (EbN0db,M,method,map);
% benoetigt pmu_psk_awgn.m, pskweigth.m, pb_psk_awgn.m
% ------------------------------------------------------------------------------
% EINGABE:
%   EbN0db: Signal/Noise pro Infobit [db]
%           (Skalar, Spalten- oder Zeilenvektor)
%        M: Stufigkeit der Modulation (default=2)
%           (Skalar)
%   method: Loesungsverfahren, das fuer M>=8-DPSK angewendet werden soll
%           'approx' (default)
%             Approximation nach [Kam96,S.473,(13.4.14)]
%           'int'
%             numerische Integration ueber die Verteilungsdichte
%             nach [Ben97,S.209,(5.24)]
%           (string)
%      map: Art des Mappings
%           'gray' (default)
%             Gray-Codierung
%           'nat'
%             natuerliches Mapping
%           (string)
%
% AUSGABE:
%   Ps: Symbolfehlerwahrscheinlichkeit (Spaltenvektor)
%   Pb: Bitfehlerwahrscheinlichkeit (Spaltenvektor)
%
% ANMERKUNGEN:
%   - exakte Loesung fuer M=2
%   - Naeherungsloesungen! fuer M>=4
%   - numerische Integration als exakte Loesung fuer M>2
%   - Integration mit Romberg-Verfahren
%   - auskommentiert: simpson-Integration mit simpson.m
%
% QUELLE:
%   [Kam96,S.483,(13.4.16)] + (vgl. pb_psk_awgn.m)
%
% AUTOR: Juergen Rinas,  08.12.1998
% Erweitert um Symbolfehlerraten von Volker Kuehn, 04.01.01
% ------------------------------------------------------------------------------

function [Ps,Pb] = pb_dpsk_koh_awgn (EbN0db,M,method,map)

if (nargin<2) M=2; end;             % default: M=2
if (nargin<3) method='approx'; end; % default: method='approx'
if (nargin<4) map='gray'; end;      % default: map='gray'
K=log2(M);
if (K~=fix(K))
  disp(['# Achtung(',mfilename,') M ist keine Potenz von 2.']);
end;

EbN0=10.^(EbN0db(:)./10);

if (M==2)
  % exakte Loesung
  Pb=erfc(sqrt(EbN0)).*(1-0.5*erfc(sqrt(EbN0)));
  Ps=Pb;
else
  if ((upper(method(1))=='A') & (upper(map(1))=='G'))
    if (M==4)
      % Naeherung
      Pb=erfc( sqrt(EbN0) );
    else
      % Naeherung fuer M>=8
      Pb=2/K*erfc( sqrt(K.*EbN0) * sin(pi/M) );
    end;
    [Ps,tmp]=pb_psk_awgn(EbN0db,M,'approx');
    Ps=Ps*2;
  else
    % numerische Integration
    sgt=pskweight(K,map);

    W=toeplitz(sgt/K);

    Pb=zeros(length(EbN0db),1);
    Ps=zeros(length(EbN0db),1);
    Pmu=pmu_psk_awgn(EbN0db,M);

    for EbN0idx=1:length(EbN0db)
      Pb(EbN0idx)=[Pmu(:,EbN0idx)]'*W*Pmu(:,EbN0idx);
      Ps(EbN0idx)=1-[Pmu(:,EbN0idx)]'*Pmu(:,EbN0idx);
    end
  end;

end;

% ### EOF ######################################################################
