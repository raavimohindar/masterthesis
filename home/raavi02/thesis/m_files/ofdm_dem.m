% ##############################################################################
% ##  ofdm_dem.m : OFDM-Demodulator mit Kanalschaetzung und Entzerrung        ##
% ##############################################################################
%
% Aufruf:     [rec_sym,ofdm] = ofdm_dem(r_sig,ofdm);       
%
% Eingabe:    r_sig : Basisband-Empfangssignal 
%
%             ofdm  : Struktur mit allen OFDM-Parametern. 
%         
%               ofdm.n_car      : Anzahl belegter Traeger  
%               ofdm.n_fft      : FFT-Laenge (optional, dafault: Anzahl
%                                 Untertraeger) 
%               ofdm.n_guard    : Guardintervall-Laenge in Abtastwerten 
%               ofdm.no_dc      : Flag zur Unterdrueckung des Gleichanteils 
%                                 (optional, default = 0)
%               ofdm.t_shift    : zeitliche Verschiebung des Empfangssignals 
%                                 (default = 0) 
%               ofdm.equal      : Multiplikative Entzerrung 
%                                 (optional, default: 1 = an) 
%               ofdm.ref_sym    : Matrix mit Referenzsymbolen 
%                                 (bekannte Symbole am Burstanfang) 
%               ofdm.ch_time    : Vektor mit Kanalimpulsantwort 
%               ofdm.ch_freq    : Vektor mit Kanaluebertragungsfunktion 
%                                 (Laenge: ofdm.n_fft oder ofdm.n_car) 
%               ofdm.l_max      : Laenge der Impulsantwort fuer 
%                                 Rauschreduktion 
%                                 (default: 0 == keine Rauschred.)
%
% Ausgabe:    r_sym : Sendesignal 
%                     
%             ofdm         : wirklich verwendete Parametereinstellungen 
%
%             ofdm.ch_freq : Vektor mit verwendeter Kanaluebertragungsfunktion 
%                            (Laenge: ofdm.n_fft oder ofdm.n_car) 
%
% Anmerkung:  Kanalschaetzung: Die Kanalschaetzung erfolgt durch Vergleich der
%             empfangenen Symbole mit den in ofdm.ref_sym befindlichen
%             Symbolen. Ist diese Matrix nicht vorhanden, oder weicht ihre
%             Hoehe von ofdm.n_car ab, so wird keine Kanalschaetzung
%             durchgefuehrt. Dann wird geprueft, ob der Vektor ofdm.ch_time
%             (Kanalimpulsantwort) vorhanden ist. Wenn ja, werden daraus die
%             Kanalkoeffizienten berechnet. Ist dieser Vektor auch nicht
%             vorhanden, dann wird geprueft, ob ofdm.ch_freq vorhanden ist,
%             der die Uebertragungsfunktion (FFT-Laenge) oder die
%             Untertraegerkoeffizienten (Laenge ofdm.n_car) enthalten muss. 
%             Reihenfolge:  1) ofdm.ref_sym   2) ofdm.ch_time   3) ofdm.ch_freq 
%
% ------------------------------------------------------------------------------
%
%  Author:  Heiko Schmidt (University of Bremen)
%  Project:  
%  Date:  12-Jun-2001
%  Last modified:   12-Jun-2001 (by H. Schmidt)
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


function [rec_sym,ofdm] = ofdm_dem(r_sig,ofdm);

% Syntax Check
if exist('ofdm') ~=1 
  ofdm.no_dc = 0;
end;

if exist('r_sig') ~=1 
  r_sig = -1;
end;
r_sig = r_sig(:); 

if (isfield(ofdm,'n_car')~= 1)
  ofdm.n_car        = 64; 
end;
if (isfield(ofdm,'n_fft')~= 1)
  ofdm.n_fft        = ofdm.n_car; 
end;
if (isfield(ofdm,'n_guard')~= 1)
  ofdm.n_guard        = 0; 
end;
if (isfield(ofdm,'no_dc')~= 1)
  ofdm.no_dc        = 0; 
end;
if (isfield(ofdm,'t_shift')~= 1)
  ofdm.t_shift      = 0; 
end;
if (isfield(ofdm,'l_max')~= 1)
  ofdm.l_max       = 0; 
end;
if (isfield(ofdm,'ref_sym')~= 1)
  ofdm.ce      = 0;
else
  [rsh,rsb] = size(ofdm.ref_sym); 
  if rsh ~= ofdm.n_car 
    ofdm.ce = 0; 
  else
    ofdm.ce = 1; 
  end; 
end;
if (isfield(ofdm,'equal')~= 1)
  ofdm.equal      = 1; 
end;

% Untertraeger-Mapping
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

sig_len = length(r_sig); 

if ofdm.t_shift ~=0 
  if ofdm.t_shift > 0 
    r_sig = [r_sig((1+ofdm.t_shift):sig_len);...
             r_sig(1:ofdm.t_shift)]; 
  else 
    r_sig = [r_sig((sig_len+ofdm.t_shift+1):sig_len);...
             r_sig(1:(sig_len+ofdm.t_shift))]; 
  end;
end;

% FFT-Eingangsmatrix bilden und Guardintervall entfernen
n_sps     = ofdm.n_fft + ofdm.n_guard; 
n_sym     = floor(sig_len/n_sps); 
r_sig_mat = zeros(n_sps,n_sym); 
r_sig_mat(:) = r_sig(1:(n_sps*n_sym)); 
r_sig_mat = r_sig_mat((ofdm.n_guard+1):n_sps,:);

% FFT-Eingangsmatrix bilden
w_fft     = ofdm.n_fft / ofdm.n_car;
r_sym_mat = fft(r_sig_mat,ofdm.n_fft) * sqrt(1/(ofdm.n_fft*w_fft));

r_sym   = r_sym_mat(index,:);

% Kanalschaetzung 
if ofdm.ce 
  ref_sym = ofdm.ref_sym; 
  [rsh,rsb] = size(ref_sym); 
  n_ref_sym = min(rsb,n_sym); 
  if n_ref_sym > 1
    sc_coef   = mean((r_sym(:,1:n_ref_sym) ./ ref_sym(:,1:n_ref_sym)).').'; 
  else
    sc_coef   = r_sym(:,1:n_ref_sym) ./ ref_sym(:,1:n_ref_sym);
  end; 
  %  Rauschreduktion 
  if ofdm.l_max > 0
    td_coef = ifft(sc_coef,ofdm.n_fft); 
    td_coef((ofdm.l_max+1):ofdm.n_fft) = 0; 
    sc_coef = fft(td_coef,ofdm.n_fft); 
  end; 

else 
  if isfield(ofdm,'ch_time') == 1 
    ch_time_tmp = zeros(ofdm.n_fft,1); 
    ch_time_tmp(1:length(ofdm.ch_time)) = ofdm.ch_time; 
    ch_freq = fft(ch_time_tmp(:),ofdm.n_fft); 
  elseif isfield(ofdm,'ch_freq') == 1 
    ch_freq_len = length(ofdm.ch_freq);   
    ch_freq     = ofdm.ch_freq(:); 
    if (ch_freq_len ~= ofdm.n_fft) & (ch_freq_len ~= ofdm.n_car)
      ch_freq = ones(ofdm.n_fft,1); 
    end;
  else 
    ch_freq = ones(ofdm.n_fft,1); 
  end;
  if (length(ch_freq) ~= ofdm.n_car) 
    sc_coef = ch_freq(index);  
  else
    sc_coef = ch_freq; 
  end; 
  % SC-Koeffizienten korrigieren bei Time-Shift
  if ofdm.t_shift ~=0 
    sc_twiddle_tmp = exp(j*2*pi*[0:ofdm.n_fft-1].' *(ofdm.t_shift/ofdm.n_fft)); 
    sc_twiddle     = sc_twiddle_tmp(index); 
    sc_coef        = sc_coef .* sc_twiddle; 
  end; 

end; 

if ofdm.equal
  equal_mat = repmat((1./sc_coef),1,n_sym);
  rec_sym = r_sym .* equal_mat;  
else   
  rec_sym = r_sym; 
end; 

% ### EOF ######################################################################
