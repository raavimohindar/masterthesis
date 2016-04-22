% ##############################################################################
% ##  viterbi_fehlervec.m : Fehleranalyse des Viterbi-Entzerrers fuer QPSK    ##
% ##############################################################################
%
% Aufruf:    [gamma_square,e] = viterbi_fehlervec(h);
%
% Eingabe:   h = Zeilenvektor Symboltakt-Impulsantwort der Laenge n+1
%
% Ausgabe:   gamma_square = saemtliche Werte gamma^2
%            e(1:length(gamma_square),1:3) = zugehoerige Fehlervektoren
%
% Anmerkung: Fuer vorgegebene Kanalimpulsantwort:
%            Fehlervektoren bis zur Laenge Le=3
%            sortiert nach aufsteigendem gamma^2
%            Beispiel Darstellung:
%            gamma_square(1:6)   ersten 6 gammas
%            e(1:6,:)            ersten 6 Fehlervektoren (zeilenweise)

function [gamma_square,e]=viterbi_fehlervec(h);

e=zeros(1,3);
ig=0;
h=h/sqrt(h*h');
d(1)=0;
d(2)=1;
d(3)=-1;
d(4)=j;
d(5)=-j;
d(6)=1+j;
d(7)=1-j;
d(8)=-1+j;
d(9)=-1-j;

Le=3;
h(Le+length(h)-1)=0;
z=[h(1) 0 0];
H=toeplitz(h,z);


A=H'*H;
for i1=2:9

  E(1)=d(i1);

  for i2=1:9

    E(2)=d(i2);

    for i3=1:9

      E(3)=d(i3);

      ig=ig+1;
      gamma_square(ig)=conj(E)*A*E.';
      gamma_square=real(gamma_square);
      E2(ig,:)=E;
    end
  end
end
[gamma_square,n]=sort(gamma_square);
igmax=length(gamma_square);
for i=1:igmax
  e(i,:)=E2(n(i),:);
end

% ### EOF ######################################################################
