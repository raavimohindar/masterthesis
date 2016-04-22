% ##############################################################################
% ##  wurzcos.m : Erzeugung eines Signals mit Wurzel-Cos-Roll-Off             ##
% ##              Charakteristik                                              ##
% ##############################################################################
%
% Aufruf:    [t,h]=wurzcos(r,w,L);
%
% Eingabe:   r = Roll-Off-Faktor
%            w = Anzahl von Abtastwerten pro Symbolintervall (T/TA)
%            L = Laenge von h in Symbolintervallen
%
% Ausgabe:   t = Zeitvektor auf T normiert: t --> t/T = -L/2 ... +L/2
%            h = Wurzel-Cos-Roll-Off-Impuls
%

function [t,h]=wurzcos(r,w,L);

nfg=1/(2*w);
nh=L*w+1;
kappa = (nh-1)/2;
ugr=rem(nh,2);

if ugr==0,                % gerade Filterlaenge
  for i=1:nh,
    arg=2*pi*nfg*(i-kappa-1);
    c=abs(arg*4*r/pi);
    if abs(c-1)<=0.000000000000000000000000000001,
      arg1=pi*(1+r)/(4*r);
      arg2=pi*(1-r)/(4*r);
      h(i)=-2*nfg*r*((2*cos(arg1)/pi)-cos(arg2));
    else
      hz=(arg*4*r/pi)*cos(arg*(1+r))+sin(arg*(1-r));
      hn=(1-((arg*4*r)/pi)^2)*(arg/(2*nfg));
      h(i)=hz/hn;
    end
  end
else                      % ungerade Filterlaenge
  for i=1:1:nh,
    arg=2*pi*nfg*(i-kappa-1);
    c=abs(arg*4*r/pi);
    if (i-round(kappa)-1)==0,
      h(i)=2*nfg*(1+r*((4/pi)-1));
    elseif abs(c-1)<=0.000000000000000000000000000001,
      arg1=pi*(1+r)/(4*r);
      arg2=pi*(1-r)/(4*r);
      h(i)=-2*nfg*r*((2*cos(arg1)/pi)-cos(arg2));
    else
      hz=(arg*4*r/pi)*cos(arg*(1+r))+sin(arg*(1-r));
      hn=(1-((arg*4*r)/pi)^2)*(arg/(2*nfg));
      h(i)=hz/hn;
    end
  end
end

Anorm=max(conv(h,h));
h=h/sqrt(Anorm);

LL=w*L;
t=[-LL/2:LL/2]/w;

% ### EOF ######################################################################
