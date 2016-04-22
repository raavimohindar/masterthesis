% ##############################################################################
% ##  cpm_ini.m : Festlegung einiger charakteristischer Parameter eines       ##
% ##              CPM-Systems fuer die Datenuebertragung                      ##
% ##############################################################################
%
% Aufruf:    struct_cpm_ini = cpm_ini(L,impulsform,eta,M,w,<f3dBT>);
%
% Beispiel:  struct_cpm_ini = cpm_ini(4,'GAUSS',0.5,2,8,0.3);
%
% Eingabe:   L          = Impulslaenge L*T (T kennzeichnet die Symboldauer)
%            impulsform = 'REC','RC' oder 'GAUSS'
%            eta        = Modulationsindex
%            M          = Stufigkeit des gesendeten Datensignals
%            w          = Ueberabtastfaktor --> w*To = T
%                         --> Abtastfrequenz f0 = 1/To = w*1/T
%            <f3dBT>    = 3dB-Frequenz des Gaussfilters * T
%
% <.> sind optionale Angaben und werden mit default Einstellungen belegt, wenn
%     sie nicht explizit gesetzt werden. Argumente "von hinten" weglassen !!
%
%                                                   Benthin 9/91

function  struct_cpm_ini = cpm_ini(L,impulsform,eta,M,w,f3dBT);

if nargin < 6,
  f3dBT = 0.3;
end;

struct_cpm_ini = struct('L', L, 'impulsform', impulsform, 'eta', eta,...
                        'M', M, 'w', w, 'f3dBT', f3dBT);

% ### EOF ######################################################################
