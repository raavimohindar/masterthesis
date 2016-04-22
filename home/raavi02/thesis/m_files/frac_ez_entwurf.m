% ##############################################################################
% ##  frac_ez_entwurf.m : Berechnung eines T/2 Entzerrers aus der             ##
% ##                     Kanalimpulsantwort mit Nebenbedingung minimaler      ##
% ##                     Koeffizientenenergie                                 ##
% ##############################################################################
%
% Aufruf:    e=frac_ez_entwurf(h2,di0,n,gamma);
%
% Eingabe:   h2    = Kanalimpulsantwort, w=2
%            di0   = Verschiebung der "1" aus der Mittellage
%            n     = Entzerrerordnung
%            gamma = Gewichtsfaktor fuer NB e'e --> Min
%
% Ausgabe:   e     = Entzerrervektor

function e=frac_ez_entwurf(h2,di0,n,gamma);

h2=h2(:);
H=toeplitz([h2; zeros(n,1)],[h2(1) zeros(1,n)]);
[z,s]=size(H);
H2(1:round(z/2),:)=H(1:2:z,:);

clear H;
[z,s]=size(H2');
i=zeros(s,1);
im=round(s/2)+di0+1;
i(im)=1;
I=diag(ones(1,s));

e=H2'*inv(H2*H2'+gamma*I)*i;

% ### EOF ######################################################################
