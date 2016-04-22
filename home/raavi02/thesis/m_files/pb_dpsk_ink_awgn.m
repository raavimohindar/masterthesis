% ##############################################################################
% ##  pb_dpsk_ink_awgn: Bitfehlerwahrsch. fuer DPSK beim AWGN Kanal,          ##
% ##                    inkoh. Det.                                           ##
% ##############################################################################

% function [Ps,Pb] = pb_dpsk_ink_awgn (EbN0db,M,method,map);
% benoetigt pb_psk_awgn.m, pskweigth.m
% ------------------------------------------------------------------------------
% EINGABE:
%   EbN0db: Signal/Noise pro Infobit [db]
%           (Skalar, Spalten- oder Zeilenvektor)
%        M: Stufigkeit der Modulation (default=2)
%           (Skalar)
%   method: Loesungsverfahren, das fuer M>=4-DPSK angewendet werden soll
%           'approx' (default)
%             Approximation nach [Kam96]
%           'int'
%             numerische Integration
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
%   - exakte Loesungen
%   - Naeherungen fuer M>=4 DPSK
%
% QUELLE:
%   2-DPSK [Kam96,S.484,(13.4.20)][Pro95,S.275,(5-2-69)]
%   4-DPSK [Pro95,S.276,(5-2-70)] auskommentiert
%   M>=4-DPSK [Kam96,S.484] Apprximation
%   [KH94]  numerische Integration
%   [PRR82] numerische Integration
%
% AUTOR: Juergen Rinas,  08.12.1998
% Erweitert um Symbolfehlerraten von Volker Kuehn, 04.01.01
% ------------------------------------------------------------------------------

function [Ps,Pb] = pb_dpsk_ink_awgn (EbN0db,M,method,map)

if (nargin<2) M=2; end;             % default: M=2
if (nargin<3) method='approx'; end; % default: method='approx'
if (nargin<4) map='gray'; end;      % default: map='gray'
K=log2(M);
if (K~=fix(K))
  disp(['# Achtung(',mfilename,') M ist keine Potenz von 2.']);
end;

EbN0=10.^(EbN0db(:)./10);

Pb=zeros(length(EbN0db),1);

if (M==2)
  Pb=0.5 * exp(- EbN0 );
  Ps=Pb;
else
  if ((upper(method(1))=='A') & (upper(map(1))=='G'))
    % Approximation
    [Ps,Pb]=pb_psk_awgn(EbN0db-10*log10(2),M);
  elseif ((upper(map(1))=='G') & (M<=32))
    if (M==4)
      Pb=F(pi/4,M,EbN0)+F(3*pi/4,M,EbN0);

      if (0) % 4-DPSK [Pro95,S.276,(5-2-70)]
        a=sqrt(2*EbN0.*(1-sqrt(.5)));
        b=sqrt(2*EbN0.*(1+sqrt(.5)));

        Pb=zeros(length(EbN0),1);

        k=0;
        while (1)
          % Schleife notwendig, da Parameter im Exponenten
          inkrement=(a./b).^k .* besseli(k,a.*b);
          Pb=Pb+inkrement;
          if (inkrement<eps)
            break;
          end;
          k=k+1;
        end;
        Pb=(Pb-0.5 .* besseli(0,a.*b)).*exp(-0.5*(a.^2+b.^2));
      end;

    elseif (M==8)
      Pb=(2/3)*(F(pi/8,M,EbN0)+F(3*pi/8,M,EbN0));
    elseif (M==16)
      Pb=0.5*(F(pi/16,M,EbN0)+F(3*pi/16,M,EbN0)...
              +F(9*pi/16,M,EbN0)-F(13*pi/16,M,EbN0));
    elseif (M==32)
      Pb=0.4*(F(pi/32,M,EbN0)+F(3*pi/32,M,EbN0)...
              +F(9*pi/32,M,EbN0)+F(17*pi/32,M,EbN0)...
              +F(19*pi/32,M,EbN0)-F(13*pi/32,M,EbN0)...
              -F(29*pi/32,M,EbN0)-F(23*pi/32,M,EbN0));
    end
  end;  % else (M>32)

  % numerische Integration ueber die Verteilungsdichte, FFT
  coeff=pskweight(K,map)/K;

  Ps=zeros(length(EbN0db),1);
  for EbN0idx=1:length(EbN0db)
    EbN0=10.^(EbN0db(EbN0idx)./10);
    EsN0=EbN0*K;

    x=0:0.001:(2*pi);
    ftheta=1/(2*pi)...
      .*(exp(-EsN0) + sqrt(pi*EsN0).*cos(x).*exp(-EsN0*sin(x).^2) ...
         .*(1+erf(sqrt(EsN0).*cos(x)))...
         );

    % Verbesserungspotential: Transformation reeller Folgen
    % nicht durchgefuehrt, da so uebersichtlicher
    % und >=64-DPSK in der Praxis uninteressant
    FCH=2*pi*ifft(ftheta,length(x));
    FCH=FCH.^2;
    fftheta=real(fft(FCH,length(x))/2/pi);

    for i=1:length(coeff)/2
      Theta1=1/M+ 2/M *(i-1);
      Theta2=1/M+ 2/M*(i);
      if Theta2>1, Theta2=1; end;

      x1=fix(length(x)*Theta1/2);
      x2=fix(length(x)*Theta2/2);
      idx=x1:x2;
      if rem(length(idx),2)==0 % wg. Simpson Laenge anpassen
        idx(1)=[];             % resultierender Fehler wird vernachlaessigt
      end;
      int=simpson(fftheta(idx),x(2)-x(1));
      if ((M>32) | (upper(map(1))=='N'))   % Pb noch nicht berechnet
        Pb(EbN0idx)=Pb(EbN0idx)+coeff(i+1)*int;
      end
      Ps(EbN0idx)=Ps(EbN0idx) + int;
      if (i<M/2)
        if ((M>32) | (upper(map(1))=='N'))   % Pb noch nicht berechnet
          Pb(EbN0idx)=Pb(EbN0idx)+coeff(M-i+1)*int;
        end
        Ps(EbN0idx)=Ps(EbN0idx) + int;
      end
    end;
    Pb(EbN0idx);
  end;  % for EbN0idx=1:length(EbN0db)

end;


function [P] =F(x,M,EbN0);
K=log2(M);

Theta=(0:0.01:1/2)*pi;

P=zeros(length(EbN0),1);
for EbN0idx=1:length(EbN0)
  EsN0=EbN0(EbN0idx)*K;  % Es/N0
  ftheta= exp(-EsN0*(1-cos(x).*cos(Theta)))...
    .* (1./(1-cos(x).*cos(Theta)));
  P(EbN0idx)=sin(x)/2/pi*simpson(ftheta,Theta(2)-Theta(1));
end;

% ### EOF ######################################################################
