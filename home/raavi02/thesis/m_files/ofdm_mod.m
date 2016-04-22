% ##############################################################################
% ##  ofdm_mod.m : OFDM-Modulator                                             ##
% ##############################################################################
%
% Aufruf:    [s_sig,ofdm] = ofdm_mod(s_sym,ofdm);
%
% Eingabe:   s_sym: Quell-Symbole (komplex). Die Symbole muessen als Matrix mit
%                   den Dimensionen [n_car,n_sym] uebergeben werden, wobei
%                   n_car die Anzahl der belegten Untertraeger und n_sym die
%                   Anzahl der OFDM-Symbole bezeichnet. Falls zufaellige
%                   QPSK-Symbole benutzt werden sollen, ist hier 
%                   s_sym = -1 anzugeben. 
% 
%            ofdm:  Struktur mit allen OFDM-Parametern. 
%                     
%              ofdm.n_car   : Anzahl belegter Traeger (Der Wert wird durch die
%                             Hoehe der Matrix s_sym ueberschrieben, daher ist
%                             die Angabe nur bei s_sym=-1 erforderlich) 
%              ofdm.n_sym   : Anzahl der auzuwuerfelnden OFDM-Symbole 
%                             (Angabe nur bei s_sym=-1 erforderlich) 
%              ofdm.n_fft   : FFT-Laenge (optional, 
%                                         default: Anzahl Untertraeger) 
%              ofdm.n_guard : Guardintervall-Laenge in Abtastwerten 
%              ofdm.n_cos   : Laenge der weichen Cos-Roll-Off-Flanke 
%              ofdm.no_dc   : Flag zur Unterdrueckung des Gleichanteils 
%                             (optional, default: 0)
%              ofdm.n_ref   : Anzahl der in ofdm.ref_sym abzulegenden 
%                             OFDM-Symbole
%                             (default: alle Symbole)
%
% Ausgabe:   s_sig : Sendesignal 
%                       
%            ofdm            : wirklich verwendete Parametereinstellungen  
%            ofdm.ref_sym    : verwendete Symbolmatrix (n_car x n_sym) 
%                              --> Fuer Kanalschaetzung 
%                       
% ------------------------------------------------------------------------------
%
%  Author:  Heiko Schmidt (University of Bremen)
%  Project:  
%  Date:  07-Jun-2001
%  Last modified:   26-Jun-2001 (by H. Schmidt)
%  Matlab Version: 5.3
%
% -----------------------------------------------------------------------
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
% -----------------------------------------------------------------------

function [s_sig,ofdm] = ofdm_modulator(s_sym,ofdm);

% Syntax Check
if exist('ofdm') ~=1 
  ofdm.no_dc = 0;
end;
if exist('s_sym') ~=1 
  s_sym = -1;
end;

random_flag = 0;  

if length(s_sym)==1 
  random_flag = 1;
else 
  [n_car_s,n_sym_s] = size(s_sym); 
  % ofdm.n_car = n_car; 
  if n_car_s > ofdm.n_fft 
    if (isfield(ofdm,'n_car')== 1) 
      n_sym        = floor(n_car_s*n_sym_s/ofdm.n_car);
    else 
      n_sym        = floor(n_car_s*n_sym_s/ofdm.n_fft);
      ofdm.n_car   = ofdm.n_fft; 
    end; 
    s_sym_tmp    = zeros(ofdm.n_car,n_sym); 
    s_sym_tmp(:) = s_sym(1:(ofdm.n_car*n_sym)); 
    s_sym        = s_sym_tmp; 
    s_sym_tmp    = []; 
  else 
    n_sym = n_sym_s;
  end; 
  ofdm.n_sym = n_sym; 
end;

if (isfield(ofdm,'n_car')~= 1)
  ofdm.n_car        = 64; 
end;
if (isfield(ofdm,'n_sym')~= 1)
  ofdm.n_sym        = 100; 
end;
if (isfield(ofdm,'n_fft')~= 1)
  ofdm.n_fft        = ofdm.n_car; 
end;
if (isfield(ofdm,'n_guard')~= 1)
  ofdm.n_guard        = 0; 
end;
if (isfield(ofdm,'n_cos')~= 1)
  ofdm.n_cos        = 0; 
end;
if ofdm.n_cos > ofdm.n_guard 
  ofdm.n_cos = ofdm.n_guard;
end;

if (isfield(ofdm,'no_dc')~= 1)
  ofdm.no_dc        = 0; 
end;

%  Bei Bedarf --> zufaellige Symbole auswuerfeln 
if random_flag 
  s_sym = (1/sqrt(2))*(sign(randn(ofdm.n_car,ofdm.n_sym)) ...
                       + j * sign(randn(ofdm.n_car,ofdm.n_sym)));
end;

% Untertraeger-Mapping
n_sym     = ofdm.n_sym;
all_car   = ofdm.n_car; 
n_fft     = max(all_car,ofdm.n_fft);

n_zero    = n_fft - all_car;
if (n_zero == 0)
  index  = 1:all_car;
else
  all_car1 = round(all_car/2);
  all_car2 = n_fft - (all_car-all_car1);
  if ofdm.no_dc
    index = [(all_car2+1):n_fft 2:all_car1+1];
  else
    index = [(all_car2+1):n_fft 1:all_car1];
  end;
end;
%
if (isfield(ofdm,'n_ref') == 1)
  ofdm.ref_sym  = s_sym(:,1:min(n_sym,ofdm.n_ref)); 
else
  ofdm.ref_sym  = s_sym; 
end;

% IFFT-Eingangsmatrix bilden
w_fft               = ofdm.n_fft / ofdm.n_car;
s_sym_mat           = zeros(ofdm.n_fft,n_sym);   
s_sym_mat(index,:)  = s_sym;

% IFFT-Computation
s_sig_mat = ifft(s_sym_mat,ofdm.n_fft) * sqrt(ofdm.n_fft * w_fft);

% add guard interval and Cos-Roll-off-ramp
n_guard    = ofdm.n_guard;     % length of guard interval (incl. cos-roll-off)
n_cos      = ofdm.n_cos;       % length of cos-roll-off-ramp
g_sig_mat  = [s_sig_mat((n_fft - n_guard + 1:n_fft),:);...
              s_sig_mat;...
              s_sig_mat(1:n_cos,:)];
sym_len    = n_fft + n_guard + n_cos;

v_koeff    = [1:n_cos]';
ramp2      = 0.5*(1+cos(pi/(n_cos+1)*v_koeff));
ramp1      = ramp2(n_cos:-1:1);

m_ramp2    = repmat(ramp2,1,n_sym);
m_ramp1    = repmat(ramp1,1,n_sym);

g_sig_mat(1:n_cos,:) = m_ramp1 .* g_sig_mat(1:n_cos,:);
g_sig_mat((sym_len+1-n_cos):sym_len,:)=m_ramp2 ...
  .* g_sig_mat((sym_len+1-n_cos):sym_len,:);
g_sig_mat(1:n_cos,2:n_sym) = g_sig_mat(1:n_cos,2:n_sym) ...
  + g_sig_mat((sym_len+1-n_cos):sym_len,1:n_sym-1);

last_ramp = g_sig_mat((sym_len+1-n_cos):sym_len,n_sym);
g_sig_mat = g_sig_mat(1:(n_guard+n_fft),:);
s_sig     = [g_sig_mat(:) ; last_ramp];

% ### EOF ######################################################################
