% ##############################################################################
% ##  cpm_mod.m : zur Erzeugung der komplexen Einhuellenden eines             ##
% ##              CPM-Signals                                                 ##
% ##############################################################################
%
% Aufruf:    y = cpm_mod(q,d,eta,M,w,isis);
%
% Beispiel:  y = cpm_mod(q,d,0.3,4,8,isis);
%
% Eingabe:   q     = Integral (Simpson) ueber den Impulsformer im
%                    Intervall [1, L*w]
%            d     = Datenblock einer Sendedatenfolge
%                    (ganzz. Vielfaches von L=length(q)/w !)
%            eta   = Modulationsindex
%            M     = Stufigkeit des gesendeten Datensignals
%            w     = Faktor der Ueberabtastung  T = w * To
%            isis  = 0 alle Daten werden moduliert; = 1 jedes L-te Datum wird
%                    moduliert (Intersymbolinterferenzfreiheit zu Testzwecken)
%
% Ausgabe:   y     = komplexe Einhuellende
%
% Bemerkungen:
%            Der Ueberabtastfaktor w ist im folgenden Programmcode mit w 
%            bezeichnet.
%
%            isis = 0 ermoeglicht es einen einzelnen Datenimpuls nach der
%            ZF-Filterung zu betrachten. Diese Intersymbolinterferenzfreiheit
%            wird erzeugt, indem nurjedes L-te Datum auf den Kanal gegeben wird.
%            Damit sinkt die Datenrate und damit die spektrale Breite des
%            Sendesignals. Die Verformung des Grundimpulses durch das ZF-Filter
%            wird so nur ungenuegend erfasst. Bei linearphasigenZF-Filtern, kann
%            aber wenigstens geprueft werden ob der Abtastzeitpunkt richtig
%            gewaehlt wurde.
%
%                                                          Benthin 9/91

function y = cpm_mod(q,d,eta,M,w,isis);

L = round(length(q)/w);
l = length(d);
jpieta = j*pi*eta;

Lw = L*w;
kr = Lw-1;
ddlen = l/L;                 % Laenge einer Teildatenfolge
xlen = ddlen*Lw + (L-1)*w;   % Es wird um (L-1)*w verlaengert
% weil die Teilfolgen spaeter verschoben
% ueberlagert werden muessen
if isis == 1,
  LL = 1;
else
  LL = L;
end;
x = zeros(L+3,xlen);
dd = zeros(L,ddlen);
% Aufteilung der Zufallsfolge d in L
for I = 1:LL                         % Teilfolgen (jeder L-te Wert wird
  dd(I,:) = d(I:L:l);                % ausgetastet).
  kshift = (I-1)*w;
  k = 1 + kshift;  % kshift sorgt dafuer, dass vorne Nullen belassen
  % werden. Dies beruecksichtigt die sonst spaeter noetige
  % Verschiebung der Teilfolgen fuer die korrekte Ueber-
  % lagerung.
  sum = 0;
  for J = 1:ddlen,
    qh = sum  + dd(I,J)*q;
    sum = sum + dd(I,J);
    x(I,k:k+kr) = qh;
    k = k + Lw;
  end;
  if I < L,
    x(I,k:xlen) = qh(Lw)*ones(1,(L-I)*w);    % Auffuellen des
    % Restvektors  mit dem "Endwert" der jeweiligen
    % Teilfolge.
  end;
end;

t = 1:1:xlen;
%plot(t,x(1,:),t,x(2,:),t,x(3,:),t,x(4,:))
%plot(t,x(1,:))

%pause(4)

y = zeros(1,xlen);
for I = 1:L,
  y = y + x(I,:);
end;
%plot(t,y);

y = exp(jpieta*y);

% ### EOF ######################################################################
