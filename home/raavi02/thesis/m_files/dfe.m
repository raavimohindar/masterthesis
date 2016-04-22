% ##############################################################################
% ##  dfe.m : MSE-Loesung fuer FIR-DF-Entzerrer                               ##
% ##############################################################################
%
% Aufruf:    [estern, fstern, f0, s, Eb, mse] = dfe(ch, n, m, i0, var_n);
%
% Eingabe:   c(0)..c(l):  Impulsantwort des Kanals  (l: Grad des Kanals)
%            n         :  Grad des Entzerrers e
%            m         :  Grad des DF-FIR
%            i0        :  Verzoegerung
%            var_n     :  Rauschleistung
%
% Ausgabe:   e(0)..e(n)  :  Entzerrerkoeffizienten
%            f(1)..f(m)  :  DF-Koeffizienten
%            f(0)        :  Hauptkoeffizient
%            s           :  Gesamtimp.-f(0)...f(m)
%            Eb          :  Energie des verkuerzten Impulses
%            MSE(i0) = min(e,f):  Mean square error
%
% Annahmen:  var_d=1

function [estern, fstern,f0,s, Eb, mse] = dfe(ch, n, m, i0, var_n);

var_d = 1;
nx=10*n;

l=length(ch)-1;
c=[zeros(1,nx) ch zeros(1,nx)];

rxdstern = var_d*conj(c(nx+i0+1:-1:nx+i0-n+1));

for lambda=1:n+1
  rxx(lambda)=0;
  for k=0:l-lambda+1
    rxx(lambda)=rxx(lambda)+conj(c(nx+k+1))*c(nx+k+lambda);
  end
  rxx(lambda)=rxx(lambda)*var_d;
end
rxx;

Rxx = [];
%for mu=1:n+1
%  for nu=1:n+1
%    Rxx(mu,nu) = rxx(1+abs(mu-nu));
%  end
%end
%for mu=1:n+1
%  Rxx(mu,mu) = Rxx(mu,mu) + var_n;
%end

rxx(1)=rxx(1)+var_n;
Rxx=toeplitz(rxx);
Rxx;

Rxd = [];
for mu=1:m
  for nu=1:n+1
    Rxd(mu,nu) = var_d*conj(c(nx+i0+mu-nu+2));
  end
end
Rxd;

estern = rxdstern*inv(Rxx-Rxd'*Rxd/var_d);
f0=estern*rxdstern';
C = [];
for mu=1:m
  for nu=1:n+1
    C(mu,nu) = conj(c(nx+i0+mu-nu+2));
  end
end
C;

f = C*estern';
fstern = f';

ce=conv(ch,estern);

s=ce;s(i0+2:i0+m+1)=s(i0+2:i0+m+1)-f';

[m1,m2]=max(abs(s));
ss=s;ss(m2)=ss(m2)-1;
mse=var_d*ss*ss'+var_n*estern*estern';
Eb=1+fstern*f;
s=ss;

% ### EOF ######################################################################
