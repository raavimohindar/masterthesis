% ##############################################################################
% ##  pmu_psk_awgn.m : Bestimmung der paarweisen Fehlerwahrscheinlichkeiten   ##
% ##############################################################################
%
% function [Pmu] = pmu_psk_awgn (EbN0db,M);
% ------------------------------------------------------------------------------
% EINGABE:
%   EbN0db: Signal/Noise pro Infobit [db]
%           (Skalar, Spalten- oder Zeilenvektor)
%        M: Stufigkeit der Modulation (default=2)
%           (Skalar)
%
% AUSGABE:
%   Pmu: paarweise Symbolfehlerwahrscheinlichkeit (Matrix)
%        (i-te Zeile gibt Wahrscheinlichkeiten fuer Phasenfehler von
%         (i-1)*2*pi/M an,
%         jede Spalte korrespondiert mit einem Siganl-Rausch-Abstand)
%
% ANMERKUNGEN:
%   - Loesung durch Integration mit Romberg-Verfahren
%     (Simpson-Verfahren auskommentiert)
%
% AUTOR: Volker Kuehn, 04.01.01
% ------------------------------------------------------------------------------

function [Pmu] = pmu_psk_awgn(EbN0db,M)

Pmu=zeros(M,length(EbN0db));

K = log2(M);

for EbN0idx=1:length(EbN0db)
  EbN0=10.^(EbN0db(EbN0idx)./10);
  EsN0=EbN0*K;

  % Integration nach Romberg
  for m=1:M/2
    % Integrationsbereich
    Theta1=(1/M+ 2/M *(m-1))*pi;         % untere Grenze
    Theta2=(1/M+ 2/M*(m))*pi;            % obere Grenze
    if Theta2>pi, Theta2=pi; end;

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
      end;
      if (  (isnan(T(k)))... % Integrand ergibt nichts
          | (k>13)...        % zu viele Stuetzstellen
          | (abs(T(1)-T(2))<T(k)*10^-8)) % Genauigkeit erreicht
        break;
      end;
    end; % while (1)
    Pmu(m+1,EbN0idx) = T(1);
    if (m<M/2)
      Pmu(M-m+1,EbN0idx) = T(1);
    end
  end; % for m=1:M
end;  % for EbN0idx=1:length(EbN0db)
Pmu(1,:) = 1-sum(Pmu,1);

% ### EOF ######################################################################
