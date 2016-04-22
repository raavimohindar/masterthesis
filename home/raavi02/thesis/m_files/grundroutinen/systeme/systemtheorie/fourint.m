% ##############################################################################
% ##  fourint.m : Auswertung der Fouriertransformation des endlichen          ##
% ##              Impulses R an den diskreten Frequenzen fTb                  ##
% ##############################################################################
%
% Der Zeitimpuls R ist nur zu aequidistanten Zeitpunkten bekannt und
% wird dazwischen quadratisch interpoliert. Die Frequenzen fTb des
% Spektrums brauchen nicht aequidistant gesetzt zu werden.
%
% Aufruf:    S = fourint(k0,k2n,R,fTb,M,vv,To);
%
% Eingabe:   k0  = Untere Integralgrenze gemaess t0 = k0*To
%            k2n = obere Integralgrenze gemaess t2n = k2n*T
%            R   = Zeitimpuls vom Zeitpunkt 0 beginnend !
%            fTb = Vektor der diskreten Frequenzen f*Tb
%            M   = Stufigkeit der Daten, wird zur Festlegung von
%                  Tb in Bezug auf den Symboltakt T benoetigt.
%            vv  =  Anzahl der Stuetzstellen von R im Zeittakt T
%            To  = Zeitdifferenz zwischen zwei Abtastwerten von R  
%                  (vv * To = T)
%
% Ausgabe:   S   = Spektrum
%
% Es wird das Fourierintegral von R(t) approximativ ausgewertet.
% Diese Integration erfolgt mit einer modifizierten Simpsonformel. Dabei
% wird nur der Integrand quadratisch interpoliert und die Exponential-
% schwingung bezueglich der Integration exakt behandelt. Dies hat zur
% Folge, dass die Stuetzstellenanzahl bei wachsender Frequenz nicht
% mit wachsen muss, um etwa den Integranden angemmessen zu interpolieren.
% Der numerische Aufwand wird damit erheblich begrenzt.
%
% Der eigentliche Fehler, den man macht liegt in der quadratischen Inter-
% polation der AKF R(t), also der Funktion, die man fouriertransformieren
% moechte. Dieser Zugang erlaubt eine Interpretation des Approximations-
% fehlers im Frequenzbereich. Berechnet wird durch das beschriebene
% Vorgehen das exakte Spektrum der Funktion h(t) = R(t) + eps(t),
% eps(t) beschreibt den Interpolationsfehler und ist in grober Naeherung
% als periodische Funktion aufzufassen, die mit einem Rechteckfenster
% der Laenge L*T gefenstert wurde und die Periode 2*To besitzt. Diese
% vereinfachenden Annahmen lassen sich experimentell sehr gut nachempfinden.
% Neben dem Spektrum von R(t) findet man Stoerspektralanteile bei Viel-
% fachen von 1/(2*To), enstprechend eines periodischen Signals mit Grund-
% und Oberwellen.

function S = fourint(k0,k2n,R,fTb,M,vv,To);

zp = find(fTb==0);        % Behandlung der Frequenz w=0
if length(zp) == 1,
  fTb(zp) = 0.005;
end;

k2n_1 = k2n-1;
k1 = k0+1;
k2 = k0+2;

pi2lnvv = 2*pi*log(M) / (log(2)*vv);

expk2n   = exp(-j*pi2lnvv*fTb*k2n);
expk2n_1 = exp(-j*pi2lnvv*fTb*k2n_1);
expk1    = exp(-j*pi2lnvv*fTb*k1);
expk0    = exp(-j*pi2lnvv*fTb*k0);
R0 = R(k0+1);
R2n = R(k2n+1);

w_ = To/pi2lnvv * ( ones(1,length(fTb))./fTb );

wh = pi2lnvv * fTb;
wh_ = ones(1,length(fTb)) ./ wh;

coswh = cos(wh);
sinwh = sin(wh);

keven = k2:2:k2n-2;
EXPeven = exp( -j*pi2lnvv * (keven.' * fTb) );
Reven = R(k2+1:2:k2n-1);
REXPeven = Reven * EXPeven;

kodd = k1:2:k2n_1;
EXPodd = exp( -j*pi2lnvv * (kodd.' * fTb) );
Rodd = R(k1+1:2:k2n);
REXPodd = Rodd * EXPodd;

R0_expk0 = R0*expk0;
R0_expk1 = R0*expk1;
R2n_expk2n = R2n*expk2n;
R2n_expk2n_1 = R2n*expk2n_1;

Sp2 = j*(R2n_expk2n - R0_expk0);
Sp1 = R0_expk0 + 2*REXPeven + R2n_expk2n  ...
  + 2*coswh .* ( -2*REXPodd + 0.5*R0_expk1  ...
                + 0.5*R2n_expk2n_1 + (coswh.*REXPeven) );
Sp0 = -2*sinwh .* ( -2*REXPodd + R0_expk1 + R2n_expk2n_1  ...
                   + 2*coswh.*REXPeven );

S = w_.*( Sp2 + wh_.*( Sp1 + wh_.*Sp0) );

% ### EOF ######################################################################
