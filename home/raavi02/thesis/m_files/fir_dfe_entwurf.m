% ##############################################################################
% ##  fir_dfe_entwurf.m : MSE-Loesung fuer FIR-DF-Entzerrer                   ##
% ##############################################################################
%
% Aufruf:    function [e, b, b0] = fir_dfe_entwurf(f, ne, nb, di0, EsN0);
%
% Eingabe:   h(0)..h(l) : Symboltaktimpulsantwort des Kanals 
%                         (m: Grad des Kanals)
%            ne         : Grad des Entzerrers e
%            nb         : Grad des DF-FIR
%            di0        : Ablage der Koeffizienten b(0) ... b(nb) von der
%                         Mittellage
%            EsN0       : Es/N0 in dB, Annahme sigma^2_D=1
%
% Ausgabe:   e(0)..e(ne) : FIR-Entzerrerkoeffizienten
%            b(1)..b(m)  : DF-Koeffizienten
%            b0          : b(0)
%
% Anmerkung: MSE(i0) = min(e,h) :  Mean square error

function [e, b, b0] = fir_dfe_entwurf(h, ne, nb, di0, EsN0);

var_N=10^(-EsN0/10);
var_D = 1;
nx=10*ne;      % links u. recht Nullen, um flexible Verschiebung zu ermoeglichen

m=length(h)-1;
i0=round((m+ne-nb)/2);    % Verschiebung des Hauptimpulses i0

c=[zeros(1,nx) h zeros(1,nx)];

rDXstern = var_D*conj(c(nx+i0+1:-1:nx+i0-ne+1));

for lambda=1:ne+1
  rxx(lambda)=0;
  for k=0:m-lambda+1
    rxx(lambda)=rxx(lambda)+conj(c(nx+k+1))*c(nx+k+lambda);
  end
  rxx(lambda)=rxx(lambda)*var_D;
end
rxx;

rxx(1)=rxx(1)+var_N;
Rxx=toeplitz(rxx);
h=h(:);
H=toeplitz([h; zeros(ne,1)],[h(1) zeros(1,ne)]);
Rxx=H'*H;
li=length(Rxx(1,:));
iz=zeros(1,li);
iz(1)=1;
I=toeplitz(iz);
Rxx=Rxx+var_N*I;

RDX = [];
for mu=1:nb
  for nu=1:ne+1
    RDX(mu,nu) = var_D*(c(nx+i0+mu-nu+2));
  end
end

e = inv(Rxx-RDX'*RDX/var_D)*rDXstern.';

b0=rDXstern*e;
b=1/var_D*RDX*e;
%e=e.';          % e wird hier als Zeilenvektor ausgegeben
%b=b.';          % b auch

if 0
  C = [];
  for mu=1:nb
    for nu=1:ne+1
      C(mu,nu) = conj(c(nx+i0+mu-nu+2));
    end
  end
  %C;
  h = C*e';
  balt = h';

  ce=conv(h,e);
  s=ce;s(i0+2:i0+nb+1)=s(i0+2:i0+nb+1)-h';

  [m1,m2]=max(abs(s));
  ss=s;ss(m2)=ss(m2)-1;
  mse=var_D*ss*ss'+var_N*e*e';
  Eb=1+b*h;
  s=ss;
end

% ### EOF ######################################################################
