% ##############################################################################
% ##  ofdm_auge.m : OFDM-Augendiagramm                                        ##
% ##############################################################################
%
% Aufruf:     ofdm = ofdm_auge(signal, ofdm, carrier, range);       
%
% Diese Funktion stellt ein OFDM-Augendiagramm dar. Neben den in 
% ofdm_mod.m beschriebenen Parametern, sind folgende Angaben moeglich:
%
% Eingabe:    signal        : OFDM-Zeitsignal. Dieser Vektor sollte mehrere 
%                             OFDM-Symbole (ofdm.n_sym>3) enthalten. 
%             ofdm          : Struktur mit allen OFDM-Parametern. 
%             ofdm.ch_time  : Vektor mit Kanalimpulsantwort (nur optional) 
%             ofdm.yaxis    : Minimaler und Maximaler Wert der Y-Achse 
%                             (nur optional)
%             carrier       : Optional, falls nur ein Traeger geplottet werden
%                             soll ist hier die Traegernummer anzugeben.
%                             Sonst (default=-1)          
%             range         : Vektor mit den zu plottenden x-Werten 
%                             default=-1  -->  range = [-n_fft/2:n_fft/2]; 
%
% Ausgabe:    nur graphische Ausgabe 
%             ofdm          : optional 
%
% Benoetigte m-files: ofdm_dem.m
%
% ------------------------------------------------------------------------------
%
%  Author:        Heiko Schmidt (University of Bremen)
%  Project:  
%  Date:             21-Jun-2001
%  Last modified:   07-Jun-2001 (by H. Schmidt)
%  Matlab Version: 5.3
%
% ------------------------------------------------------------------------------
%
%  Copyright (c) 2001 
%  Department of Communications Engineering / University of Bremen 
%  This software is property of the University of Bremen. 
%  Unauthorized duplication and disclosure to third parties is forbidden. 
%
%  If You have received a copy of this program from other parties, please 
%  contact us. Be sure, we have no commercial interests, but we are 
%  academically interested in your version number (last modified date) 
%  and the distribution way.  
%
% ------------------------------------------------------------------------------

function ofdm = ofdm_auge(signal,ofdm,carrier,range);

% Syntax Check
if exist('ofdm') ~=1 
  ofdm.no_dc = 0;
end;

if (isfield(ofdm,'ch_time')~=1) 
  ofdm.ch_time = 1; 
end;

if (isfield(ofdm,'yaxis')~=1) 
  ofdm.yaxis = 0; 
end;

if exist('signal') ~=1 
  singal = -1;
end;


if exist('carrier') ~=1 
  carrier = -1;
end;


if exist('range') ~=1 
  range = -1;
end;
ofdm.equal   = 1; 

r_sym_mat    = []; 

if length(range) == 1
  range = -ofdm.n_fft/2 : ofdm.n_fft/2; 
end; 
for lauf = range; 
  ofdm.t_shift  = lauf; 
  r_sym         = ofdm_dem(signal,ofdm);
  [dummy,n_sym] = size(r_sym); 
  if (carrier > 0);
    r_sym_vec = r_sym((ofdm.n_car+carrier):ofdm.n_car:((n_sym-1)*ofdm.n_car));
  else
    r_sym_vec = r_sym((ofdm.n_car+1):((n_sym-1)*ofdm.n_car));
  end;                 
  r_sym_mat     = [r_sym_mat r_sym_vec(:)]; 
end;
t_vec = range/ofdm.n_fft; 
plot(t_vec,real(r_sym_mat),'b');
grid; 
xlabel('t/T_s');
ylabel('Re\{d\}'); 
if (length(ofdm.yaxis) > 1)
  axis([min(t_vec) max(t_vec) min(ofdm.yaxis) max(ofdm.yaxis)]);
end; 

% ### EOF ######################################################################
