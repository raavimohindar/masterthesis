% ##############################################################################
% ##  qpuls.m : Integration ueber Impulsantwort des Impulsformers             ##
% ##############################################################################
%
%  QPULS liefert das Integral q(t) ueber die Impulsantwort des
%  Impulsformers g(t) bei CPM-Modulation.
%
% Aufruf:    q = qpuls(L,impulsform,w,versatz,f3dBT,graph);
%
% Eingabe:   L  =  Laenge der Impulsantwort in Vielfachen des Symboltakts T
%            impulsform  = 'REC', 'RC' oder 'GAUSS'
%            w           = Ueberabtastfaktor gemaess (w*To=T); To kennzeichnet
%                          Abtasttakt, f0 = 1/To die Abtastfrequenz
%            versatz     = 1 : Versatz des Impulses um -To/2 zur Kompensation
%                          des To/2 Versatzes durch den Demodulator
%                        = 0 : kein Versatz
%            <f3dBT>     = f3dB * T  : Kenngroesse des Gaussfilters
%            <graph>     = 1 -> mit Graphik
%
%
% Ausgabe:   q           = q(1*To) ... q(L*w*To) ensprechend L*w Werten
%
%                                                   Benthin 9/91

function q = qpuls(L,impulsform,w,versatz,f3dBT,graph);


f3dBTtext = '           ';
if (strcmp(impulsform,'REC') == 1),
  t = 1/(L*w):1/(L*w):1;
  q = t;

elseif (strcmp(impulsform,'RC') == 1),
  t = 1/(L*w):1/(L*w):1;
  if (versatz == 1),
    t = t + 1/(2*L*w);
  end;
  q = (t - 1/(2*pi)*sin(2*pi*t));

elseif (strcmp(impulsform,'GAUSS') == 1),
  % Erzeugung des GMSK Impulses g(t)
  aa = sqrt(2/log(2))*pi*f3dBT;
  t =  -L/2:1/(4*w):L/2+2/(4*w);   % Es wird t im 4fachen Abtasttakt
  % bestimmt um die versetzten und nicht versetzten
  % Zeitpunkte ausrechnen zu koennen.
  g = 1/2*( erf(aa*(t+1/2)) - erf(aa*(t-1/2)) );

  lq = 2*L*w+1;
  q = zeros(1,lq);

  sum = 0;                      % Integration des GMSK Impulses fuehrt auf q(t)
  kk = 0;
  for k=2+1:2:length(t),
    kk = kk + 1;
    sum = sum + g(k-2) + 4*g(k-1) + g(k);
    q(kk) = sum;
  end;

  qL = q(length(q)-1);
  q = 1/qL * q;                 % Normierung, so dass q(LT) = 1 gilt
  if (versatz == 1),
    q = q(3:2:lq);
  else
    q = q(2:2:lq);
  end;

  f3dBTtext = sprintf('%.3f',f3dBT);
  f3dBTtext = [', f3dBT*T = ' f3dBTtext];
else
  disp('Impulsform unbekannt');
end;

if exist('graph') == 1,
  t = 1/w:1/w:L;
  plot(t,q);
  title('"Integral des Impulsformers"');
  ylabel('q(t)     Magnitude');
  text = [int2str(L),' ',impulsform,f3dBTtext,'               t/T'];
  xlabel(text);
end;

% ### EOF ######################################################################
