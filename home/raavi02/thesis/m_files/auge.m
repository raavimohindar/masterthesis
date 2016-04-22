% ##############################################################################
% ## auge.m : erzeugt ein Augendiagramm                                       ##
% ##############################################################################
%
% Aufruf:      auge(x,w,o)
%
% Eingabe:     x = Vektor des abgetasteten Datensignals
%              w = Abtastintervalle pro Symbol ( T/TA ); w = gerade Zahl
%              o = Rechts/links-Verschiebung um o-Abtastintervalle; 
%                  o = -w ... +w
%                  in Verbindung mit datensig ergibt sich
%                  fuer o=0 ein symmetrisches Auge
%
%                                                      Kammeyer '00

function auge(x,w,o);

w1=round(w/8);                  % Augendarstellung um T/8 ueber T/2 hinaus
w2=2*w1;
LL=length(x);
o=-o;
if abs(o)>w
  disp('Verschiebung um mehr als +/- T (d.h. +/- w) unzulaessig!')
else
  x=x(o+1-w1+3*w/2:LL);         % Verschiebung, so dass Darstellung bei o=0
  L=floor(LL/w)-3;              % auf -(T/2+T/8) ...  +(T/2+T/8) fuehrt
  A=zeros(w,L);
  A(:)=x(1:L*w);
  for i=1:w2+1                  %  Anfuegen von w/4 Samples
    A=[A;[A(i,2:L) x(L*w+i)]];  %  zur Augendarstellung von +/- T/8
  end                           %  ueber +/- T/2 hinaus
  t=[-w/2:w/2+w2]'+o-w1;
  t=t/w;
  plot(t,A)
  xlabel('t/T \rightarrow')
  %title(' Augendiagramm')
  grid
  axis([min(t) max(t) -1.2*max(x) 1.2*max(x)])
end

% ### EOF ######################################################################
