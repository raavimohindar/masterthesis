% ##############################################################################
% ##  parallelconcat.m : Parallele verk. Codes mit Uniform-Interleaver        ##
% ##                     (Laenge N)                                           ##
% ##############################################################################
%
% function ACwhp=parallelconcat(ACwh1,ACwh2,N,in_type)
% ------------------------------------------------------------------------------
% EINGABE:
%  Awh1,ACwh2: Matrix der Input Redundancy Weight Enumerating Function (IRWEF)
%              w=Infowortgewicht (Matrix-Index 1 entspricht w=0 !!!)
%              h=Redundanzgewicht (Matrix-Index 1 entspricgt h=0 !!!)
%              (Matrix)
%  N:          Interleaverlaenge
%              (Skalar)
%  in_type:    Art der Eingangsspektren ('IOWEF' oder 'IRWEF')
%
% AUSGABE:
%   ACwhp: Matrix der resultierenden IOWEF
%          (Matrix)
%
% ANMERKUNGEN:
%   - benoetigt Datei bin_coef
%   - Matrizen werden in ihrer Groesse angepasst, um Verknuepfung zu
%     ermoeglichen
%   - Routine funktioniert nur fuer systematische Codes, da von der
%     zweiten IOWEF der systematische Teil entfernt wird
%     (siehe letzte Schleife)
%
% QUELLE:
%   [BR96,S.412,(7)]
%
% AUTOR: Juergen Rinas,  31.05.1999
% Eingang jetzt IRWEF, Ausgang um systematischen Anteil ergaenzt
% von Volker Kuehn, 08.05.01
% ------------------------------------------------------------------------------

function [ACwhp]=parallelconcat(ACwh1,ACwh2,N,in_type)

if (nargin<4)
  disp('Fehler: Zu wenig Eingabeargumente!')
end


wmax1=size(ACwh1,1);
wmax2=size(ACwh2,1);

wmax=wmax1;

if (wmax1>wmax2)
  disp([mfilename,': ungleiche maximale Eingangsgewichte (wmax1>wmax2)']);
  ACwh1(wmax2+1:wmax1,:)=[];
  wmax=wmax2;
end;
if (wmax2>wmax1)
  disp([mfilename,': ungleiche maximale Eingangsgewichte (wmax2>wmax1)']);
  ACwh2(wmax1+1:wmax2,:)=[];
  wmax=wmax1;
end;

if (wmax~=N+1)
  disp([mfilename,': maximales Eingangsgewicht ungleich Interleaverlaenge ',...
        '(wmax~=N+1)']);
end;

hmax1=size(ACwh1,2);
hmax2=size(ACwh2,2);

b=bin_coef(N*ones(wmax,1),0:wmax-1);

% Ausgabematrix zur Geschwindigkeitssteigerung vorbesetzen
if strcmp(upper(in_type),'IRWEF')
  ACwhp=zeros(wmax,hmax1+hmax2+wmax);
  for w=1:wmax
    tmp = conv(ACwh1(w,:),ACwh2(w,:))/b(w);
    ACwhp(w,w:w+length(tmp)-1) = tmp;          % Infogewicht einmal addieren
  end
else
  ACwhp=zeros(wmax,hmax1+hmax2-wmax);
  for w=1:wmax
    ACwr2 = ACwh2(w,w:hmax2-wmax+w);           % Infobits nur einmal uebertragen
    ACwhp(w,:)=conv(ACwh1(w,:),ACwr2)/b(w);
  end;
end;

% ### EOF ######################################################################
