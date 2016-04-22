% ##############################################################################
% ##  pb_psk_awgn: Bitfehlerwahrsch. fuer PSK beim AWGN Kanal, koh. Detektion ##
% ##############################################################################
%
% function [Ps,Pb] = pb_psk_awgn (EbN0db,M,method,map);
% benoetigt pmu_psk_awgn.m, pskweight.m, simspon.m
% ------------------------------------------------------------------------------
% EINGABE:
%   EbN0db: Signal/Noise pro Infobit [db]
%           (Skalar, Spalten- oder Zeilenvektor)
%        M: Stufigkeit der Modulation (default=2)
%           (Skalar)
%   method: Loesungsverfahren, das fuer M>8-PSK angewendet werden soll
%           'approx' (default)
%             Approximation nach [Kam96,S.473,(13.4.14)]
%           'int'
%             numerische Integration ueber die Verteilungsdichte
%             nach [IS91]
%           (string)
%      map: Art des Mappings  (string)
%           'gray' (default)
%             Gray-Codierung
%           'nat'
%             natuerliches Mapping
%
% AUSGABE:
%   Ps: Symbolfehlerwahrscheinlichkeit (Spaltenvektor)
%   Pb: Bitfehlerwahrscheinlichkeit (Spaltenvektor)
%
% ANMERKUNGEN:
%   - exakte Loesung fuer BPSK, QPSK und 8-PSK
%     hoehere Stufigkeiten mit numerischer Integration
%   - Naeherung fuer M-PSK (M>8)
%   - exakte Loesung durch Integration fuer (M>8) mit Romberg-Verfahren
%     (Routine pmu_psk_awgn.m erforderlich)
%
% QUELLE:
%   2-PSK(BPSK) [Kam96,S.473,(13.4.4c)] [Pro95,S.258,(5-2-5)]
%               [Pro95,S.271,(5-2-57)]
%   4-PSK(QPSK) [Kam96,S.473,(13.4.6)]
%   8-PSK       [IS91,S.351,(24)]
%
%   M-PSK       [Kam96,S.473,(13.4.14)] 'approx'
%   M-PSK       [Lee86,S.490] [IS91] [BBC87,S.209] 'int'
%
% AUTOR: Juergen Rinas,  08.12.1998
% Erweitert um Symbolfehlerraten von Volker Kuehn, 05.01.01
% ------------------------------------------------------------------------------

function [Ps,Pb] = pb_psk_awgn(EbN0db,M,method,map)

if (nargin<2) M=2; end;              % default: M=2
if (nargin<3) method='approx'; end;  % default: method='approx'
if (nargin<4) map='gray'; end;       % default: map='gray'
K=log2(M);
if (K~=fix(K))
  disp(['# Achtung(',mfilename,') M ist keine Potenz von 2.']);
end;

EbN0=10.^(EbN0db(:)./10);
Pmu=[];

if (M==2)
  Pb=0.5 * erfc( sqrt(EbN0) );
  Ps=0.5 * erfc( sqrt(EbN0) );
elseif (M==4)
  if (upper(map(1))=='G')
    Pb=0.5 * erfc( sqrt(EbN0) );
  else
    Pe=0.5*erfc( sqrt(EbN0) );
    Pb=3*Pe.*(1-Pe) + Pe.^2;
  end
  Ps=erfc( sqrt(EbN0) ) - (erfc( sqrt(EbN0) )/2).^2;
elseif (M>=8)
  if ((upper(method(1))=='A') & (upper(map(1))=='G'))
    % Approximation
    Ps=erfc(sqrt(K.*EbN0)*sin(pi/M));
    Pb=Ps/K;
  else
    if ((M==8) & (upper(map(1))=='G'))
      d3=sqrt(2*K*EbN0);
      Pb=2/3*Q(d3*sin(pi/8)).*(1-Q(d3*sin(3*pi/8)))+2/3*Q(d3*sin(3*pi/8));
    elseif ((M==16) & (upper(map(1))=='G'))
      d4=sqrt(2*K*EbN0);
      Pb=2/4*Q(d4*sin(pi/16)).*(1 + 2/4*Q(d4*sin(7*pi/16)))...
        +2/4*Q(d4*sin(3*pi/16)).*(1 + 2/4*Q(d4*sin(5*pi/16)))...
        -2/4*C(7,M,EbN0)-2/4*C(8,M,EbN0);
    elseif ((M==32) & (upper(map(1))=='G'))
      d5=sqrt(2*K*EbN0);
      Pb=2/5*Q(d5*sin(pi/32)).*(1 + 2.5/2*Q(d5*sin(15*pi/32)))...
        +2/5*Q(d5*sin(3*pi/32)).*(1 + 2.5/2*Q(d5*sin(13*pi/32)))...
        +1/5*Q(d5*sin(9*pi/32)).*(1 - 1.5*Q(d5*sin(7*pi/32)))...
        +1/5*Q(d5*sin(11*pi/32)).*(1 - 1.5*Q(d5*sin(5*pi/32)))...
        -1/5*Q(d5*sin(13*pi/32))-1/5*Q(d5*sin(15*pi/32))...
        +1/5*C(13,M,EbN0)+1/5*C(14,M,EbN0)-3/5*C(15,M,EbN0)-3/5*C(16,M,EbN0);
    else
      % numerische Integration aller Wahrscheinlichkeiten
      % wird verwendet, wenn M>32

      % Gewichtung der Einzelintegrale berechnen
      coeff=pskweight(K,map);

      Pmu=pmu_psk_awgn(EbN0db,M);
      Pb = coeff'*Pmu/K;
    end;  % else (M>32)

    % Exakte Berechnung der Symbolfehlerrate
    if isempty(Pmu)
      Pmu=pmu_psk_awgn(EbN0db,M);
    end
    Ps = 1 - Pmu(1,:);
  end;  % else (exakte Berechnung)
end;  % elseif (M>=8)



% zur Uebersicht wird die Q-Funktion aus [Pro95,S.40,(2-1-98)] benutzt
function [Qx]=Q(x);
Qx=0.5*erfc(x/sqrt(2));


function [P]=C(idx,M,EbN0);
m=M/2+1-idx;
K=log2(M);

P=zeros(length(EbN0),1);

if (0)
  % Integration nach Simpson

  Theta1=(1-(2*m-1)/M)*pi;
  Theta2=(1+(2*m-1)/M)*pi;
  x=(Theta1:(Theta2-Theta1)/1000:Theta2);

  for EbN0idx=1:length(EbN0),
    EsN0=EbN0(EbN0idx)*K;  % Es/N0
    ftheta=1/(2*pi)...
      .*(exp(-EsN0) + sqrt(pi*EsN0).*cos(x).*exp(-EsN0*sin(x).^2) ...
         .*(1+erf(sqrt(EsN0).*cos(x)))...
         );
    P(EbN0idx)=simpson(ftheta,x(2)-x(1));
  end;
end;

if (1)
  % Integration nach Romberg
  % Integrationsbereich
  Theta1=(1-(2*m-1)/M)*pi;
  Theta2=(1+(2*m-1)/M)*pi;
  for EbN0idx=1:length(EbN0),
    EsN0=EbN0(EbN0idx)*K;  % Es/N0

    % Integration ueber die Teilintervalle nach Romberg
    % (adaptive Erhoehung der Stuetzwerte)
    h=Theta2-Theta1;  % Schrittweite
    x=[Theta1 Theta2];
    ftheta=1/(2*pi)...
      .*(exp(-EsN0) + sqrt(pi*EsN0).*cos(x).*exp(-EsN0*sin(x).^2) ...
         .*(1+erf(sqrt(EsN0).*cos(x)))...
         );

    n=1; % Anzahl der neuen Werte im naechsten Schritt
    k=1; % Anzahl der Verfeinerungen
    T(k)=sum(ftheta)*h/2; % erster Wert laut Trapez-Formel
    while (1)
      n=n*2;
      k=k+1;  % Verfeinerungsstufe+1
      h=h/2;  % Schrittweite halbieren

      x=Theta1+h*(1:2:(n-1));
      ftheta=1/(2*pi)...
        .*(exp(-EsN0) + sqrt(pi*EsN0).*cos(x).*exp(-EsN0*sin(x).^2) ...
           .*(1+erf(sqrt(EsN0).*cos(x)))...
           );
      T(k)=T(k-1)/2+sum(ftheta)*h; % neue Werte hinzufuegen

      % Richardson Extrapolation
      p=1;
      for ii=(k-1):(-1):1,
        p=p*4;
        T(ii)=T(ii+1) + ( T(ii+1)-T(ii) )/(p-1);
      end; % for
      if (  (isnan(T(k)))... % Integrand ergibt nichts
          | (k>13)... % Schleifenabbruch wg. zu vieler Stuetzstellen
          | (abs(T(1)-T(2))<T(k)*10^-8)) % Genauigkeit erreicht
        break;
      end;
    end; % while (1)
    P(EbN0idx)=T(1);
  end;
end;

% ### EOF ######################################################################
