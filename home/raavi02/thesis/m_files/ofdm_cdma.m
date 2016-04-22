% ##############################################################################
% ## ofdm_cdma.m : Simulation eines uncodierten {OFDM}-{CDMA} Systems, liefert##
% ##               Bitfehlerrate, perfektes Interleaving: stat. unabhaengige  ##
% ##               Kanalimpulsantworten fuer jedes OFDM Symbol,Walsh-Spreizung##
% ##############################################################################
%
% Aufruf:      ber = ofdm_cdma(ofdm,cdma,channel,sim)
%
% Eingabeparameter
%  ofdm       :    struct mit folgenden OFDM-Parametern:
%    N_c      :      Anzahl der OFDM-Untertraeger 
%    N_guard  :      Laenge des Guard-Intervalls
%    N_d      :      Anzahl Bit je OFDM-Symbol
%    IL       :      Permutationsvektor fuer Frequenzinterleaver
%  cdma       :    struct mit folgenden CDMA-Parametern:
%    G_p      :      Prozessgewinn 
%    N_user   :      Anzahl Teilnehmer
%    equalize :      string gibt Art der Entspreizung an ('MRC','ORC',EGC')
%    codes    :      Matrix mit Spreizungscodes fuer die einzelnen Teilnehmer 
%                    (spaltenweise angeordnet, G_p Zeilen, falls 
%                    Kanalschaetzung, Code in letzter Spalte fuer Pilotkanal)
%  channel    :    struct mit den Kanalparametern fuer ant_channel.m
%  sim        :    struct mit folgdenden Simulationsparametern:
%    EbN0dB   :      Signal-Rausch-Abstaende
%    N_info   :      Anzahl OFDM-Symbole mit konstantem Kanal
%    max_info :      maximale Anzahl Informationsbit (optional, default=1e4)
%    max_err  :      maximale Anzahl Bitfehler (optional, default=100)
%
% Ausgangssignal
%  ber      :      Vektor der Bitfehlerrate, Laenge entspricht der von EbN0dB
%
% Volker Kuehn, 09.09.01

function ber=ofdm_cdma(ofdm,cdma,channel,sim)

if (nargin<4)
   error('ofdm_cdma: Zu wenig Eingabeparameter!');
end

if (length(fieldnames(ofdm))<4)
   error('ofdm_cdma: Zu wenig OFDM-Parameter!');
end

if (length(fieldnames(cdma))<4)
   error('ofdm_cdma: Zu wenig CDMA-Parameter!');
end

if (length(fieldnames(sim))<4)
   error('ofdm_cdma: Zu wenig Simulationsparameter!');
end

% Speicher fuer Fehlerraten  
ber = zeros(length(sim.EbN0dB),1);

% Interleaver fuer gesamten Block aus sim.N_info OFDM-Symbolen
IL = repmat(ofdm.IL(:),1,sim.N_info+1) + ones(ofdm.N_c,1)*...
               [0:ofdm.N_c:(ofdm.N_c)*(sim.N_info)];
IL = IL(:);
  
  
% Starte Simulation ueber alle SNR
for snr = 1:length(sim.EbN0dB)
   EbN0 = 10^(sim.EbN0dB(snr)/10);
   %  Mismatching durch Guard-Intervall
   EsN0 = EbN0 * ofdm.N_c / (ofdm.N_c+ofdm.N_guard);
   
   bit_cntr = 0;
   err_cntr = 0;
   
   while ((bit_cntr<sim.max_info) & (err_cntr<sim.max_err))
      u = sign(randn(cdma.N_user,ofdm.N_d*(sim.N_info+1))); % Informationsbit
      
      X = cdma.codes * u;                         % DS-Spreizung

      X = X(IL);                                  % Frequenzinterleaving

      X = reshape(X,ofdm.N_c,sim.N_info+1);       % Aufteilung auf OFDM-Symbole

      x = ifft(X,ofdm.N_c,1);                     % Transformation in Zeitbereich 

      
      % Guardintervall voranstellen
      x = [x(ofdm.N_c-ofdm.N_guard+1:ofdm.N_c,:); x]; 

      % L-Pfad-Rayleigh-Kanal         
      [y,ch] = ant_channel(x(:),channel);         % Uebertragung ueber den Kanal
      CH = fft(ch.impulse_response,ofdm.N_c,1) / sqrt(ofdm.N_c);         

      % AWGN-Kanal
      awgn = (randn(size(y)) + j*randn(size(y))) / sqrt(2*EsN0);
      y    = y + awgn;

      % Entfernung des Guard-Intervalls
      y = reshape(y,ofdm.N_c+ofdm.N_guard,sim.N_info+1);

      y = y(ofdm.N_guard+1:ofdm.N_guard+ofdm.N_c,1:sim.N_info);         

      % Transformation zurueck in Frequenzbereich
      Y = fft(y,ofdm.N_c,1) / sqrt(ofdm.N_c);         

      CH = repmat(CH(:),1,sim.N_info);
      % Entzerrung mit
      switch upper(cdma.equalize)
      case ('MRC')                              % Maximum Ratio Combining 
         X_hat = real(Y .* conj(CH));
      case ('ORC')                              % Orthogonal Restoring Combining
         X_hat = real(Y ./ CH);     
      case ('EGC')                              % Equal Gain Combining   
         X_hat = real(Y .* conj(CH) ./ abs(CH));     
      end

      X_hat(IL(1:ofdm.N_c*sim.N_info)) = X_hat;   % Frequenzdeinterleaving
      % DS-Entspreizung 

      u_hat = sign(cdma.codes.' * reshape(X_hat,cdma.G_p,ofdm.N_d*sim.N_info));
      
      err = find(u_hat~=u(:,1:sim.N_info*ofdm.N_d));

      bit_cntr = bit_cntr + cdma.N_user*ofdm.N_d*sim.N_info;
      
      err_cntr = err_cntr + length(err);
      
   end % while
   
   ber(snr) = err_cntr / bit_cntr;
   
end % for   

% ### EOF ######################################################################
