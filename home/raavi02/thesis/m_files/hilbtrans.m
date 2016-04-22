% ##############################################################################
% ##  hilbtrans.m : Hilberttransformator (interner Remez-Entwurf)             ##
% ##############################################################################
%
% Aufruf:    [t,xH] = hilbtrans(t,x,del_F,n);
%
% Eingabe:   t = Zeitvektor
%            x = zu transformierender Signalvektor
%            del_F = Beginn/Ende des Durchlassbandes des
%            Hilberttransformators (del_F ... 1-del_F)
%            n = Ordnung des Hilbert-Transformators
%
% Ausgabe:   t = Zeitvektor
%            xH = hilberttransformiertes Signal

function [t,xH] = hilbtrans(t,x,del_F,n);

F=[del_F 1-del_F];
A=[1 1];
hH=remez(n,F,A,'hilbert');

hH = hH(:);
Lh = length(hH);
Lx = length(x);

if rem(Lh,2)==0
  error(['ERROR (hilbtrans.m): Der Hilberttrafo muss eine gerade Ordnung, ',...
         'd.h.ungerade Laenge haben.']);
end;

hH = [hH((Lh-1)/2+1:Lh); zeros(Lx-Lh,1); hH(1:(Lh-1)/2)];
X  = fft(x);
HH = fft(hH);
xH = ifft(X(:).*HH(:));

if size(x,2)~=1
  xH = xH.';
end;

% ### EOF ######################################################################
