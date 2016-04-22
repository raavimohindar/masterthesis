% ##############################################################################
% ##  pr_vorcod.m : Partial-Response-Vorcodierung                             ##
% ##############################################################################
%
% Aufruf:    d_vor = pr_vorcod(b,alpha,ausg_log)
%
% Eingabe:   b        = Quelldaten in logischer Schreibweise (z.B. {0 1 1 0 1})
%            alpha    = Partial-Response Koeffizienten (z.b. {1 0 -1})
%            ausg_log = Paramter zum Ausgabeformat von d_vor (siehe Ausgabe)
%
% Ausgabe:   d_vor    = vorcodierte Daten
%                       ausg_log=1 -->  logische  Schreibweise
%                       ausg_log=0 -->  +/- Schreibweise
%
% Anmerkung: d_vor(-1)...d_vor(-n+1) werden auf 0 gesetzt
%            Ausgabelaenge von d_vor = length(b)

function d_vor=pr_vorcod(b,alpha,ausg_log);

n=length(alpha);
N=length(b);

% Erzeugung der logischen PR-Koeffizienten
alpha_log=mod(alpha,2); % alpha:   gerade --> aplha_log: 0
% alpha: ungerade --> alpha_log: 1

% Umsetzung der logischen Daten b(i) in +/- Darstellung
b=-sign(b-0.5);         % 0 --> 1  /  1 --> -1

% PR-Vorcodierung
for i=1:n-1
  d_v(i)=1;
end

d_vor=zeros(1,N);

for i=1:N
  d_vor(i)=b(i);
  for nu=1:n-1
    if alpha_log(nu+1)==1
      d_vor(i)=d_v(nu)*d_vor(i);
    else
    end
  end
  for nu=n-1:-1:2
    d_v(nu)=d_v(nu-1);
  end
  d_v(1)=d_vor(i);
end

if ausg_log==1
  d_vor=-d_vor;
  d_vor=0.5*(d_vor+1);
else
end

% ### EOF ######################################################################
