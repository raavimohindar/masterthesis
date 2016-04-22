% ##############################################################################
% ##  h_trafo.m : Hilberttransformator (externer Remez-Entwurf)               ##
% ##############################################################################
%
% Aufruf:    [t,xH] = h_trafo(t,x,hH)
%
% Eingabe:   t = Zeitvektor
%            x = zu transformierender Signalvektor
%            hH = nichtkausale Impulsantwort des Hilberttransformators
%
% Ausgabe:   t = Zeitvektor
%            xH = hilberttransformiertes Signal

function [t,xH] = h_trafo(t,x,hH);

hH = hH(:);
Lh = length(hH);
Lx = length(x);

if rem(Lh,2)==0
  error(['ERROR (h_trafo.m): Das dritte Argument muss eine ungerade ',...
         'Laenge haben.']);
end;

hH = [hH((Lh-1)/2+1:Lh); zeros(Lx-Lh,1); hH(1:(Lh-1)/2)];
X  = fft(x);
HH = fft(hH);
xH = ifft(X(:).*HH(:));

if size(x,2)~=1
  xH = xH.';
end;

% ### EOF ######################################################################
