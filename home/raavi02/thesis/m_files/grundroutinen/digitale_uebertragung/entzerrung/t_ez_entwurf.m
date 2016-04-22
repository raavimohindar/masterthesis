% ##############################################################################
% ##  t_ez_entwurf.m : Zero-Forcing oder MMSE-Entwurf der Impulsantwort       ##
% ##                   eines T-Entzerrers                                     ##
% ##############################################################################
%
% Aufruf:   e = t_ez_entwurf(Typ,f,n,di0,EsN0);
%
% Eingabe:  Typ  = Loesungstyp
%                  'ZF' = Zero Forcing
%                  'MMSE' = MMSE-Entwurf
%           h    = Kanalimpulsantwort im Symboltakt
%           n    = Entzerrerordnung
%           di0  = Ablage der Verzoegerung i0 von der Mitte
%           EsN0 = Es/N0-Verhaeltnis in dB  (nur relevant fuer MMSE-Loesung)
%
% Ausgabe:  e    = Entzerrervektor (w=1)

function e=t_ez_entwurf(Typ,h,n,di0,EsN0)

sigN_sigS=10^(-EsN0/10);    % sigmaN/sigmaS fuer MMSE-Loesung
h=h(:);
H=toeplitz([h; zeros(n,1)],[h(1) zeros(1,n)]);
i=zeros(1,n+length(h));
i0=fix(length(i)/2)+di0 +1;
i(i0)=1;
R=H'*H;
li=length(R(1,:));
iz=zeros(1,li);
iz(1)=1;
I=toeplitz(iz);

if (strcmp(Typ,'ZF') == 1),

  gamma=0;
elseif (strcmp(Typ,'MMSE') == 1),

  gamma=sigN_sigS;
else
  disp('ungueltige Typ-Eingabe --> gamma=0')
  gamma=0;
end
e=inv(R+gamma*I)*H'*i';

% ### EOF ######################################################################
