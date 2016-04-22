% ##############################################################################
% ## sc_cdma.m : Simulation eines uncodierten SC-CDMA-Systems,                ##
% ##             liefert Bitfehlerrate                                        ##
% ##             benoetigte m-files:  ant_channel.m                           ##
% ##############################################################################
%
% Aufruf:   ber = sc_cdma(cdma,channel,sim)
%
% Eingabeparameter
%  cdma       :  struct mit folgenden CDMA-Parametern:
%    G_p      :    Prozessgewinn 
%    N_user   :    Anzahl Teilnehmer
%    codes    :    Matrix mit Spreizungscodes fuer die einzelnen Teilnehmer 
%                  (spaltenweise angeordnet (G_p Zeilen, falls 
%                   Kanalschaetzung, Code in letzter Spalte fuer Pilotkanal)
%    L_rake   :    Anzahl Rake-Finger 
%    ch_est   :    Kanalschaetzung: ==0 - ideal bekannte Impulsantwort 
%                                   >0 - Schaetzung ueber Pilotkanal und 
%                                        Mittelung ueber ch_est Werte
%                                        (!!! ch_est<=N_info !!!)
%                  Beruecksichtigung der L_rake ersten Verzoegerungen !
%  channel    :  struct mit den Kanalparametern fuer ant_channel.m
%  sim        :  struct mit folgdenden Simulationsparametern:
%    EbN0dB   :    Signal-Rausch-Abstaende
%    N_info   :    Anzahl Informationsbit mit konstantem Kanal
%    max_info :    maximale Anzahl Informationsbit (optional, default=1e4)
%    max_err  :    maximale Anzahl Bitfehler (optional, default=100)
%
% Ausgangssignal
%  ber      :    Vektor der Bitfehlerrate, Laenge entspricht der von EbN0dB
%
% Volker Kuehn, 09.09.01

function ber=sc_cdma(cdma,channel,sim)

  if (nargin<3)
      error('sc_cdma: Zu wenig Eingabeparameter!');
  end
  
  if (length(fieldnames(cdma))<5)
      error('sc_cdma: Zu wenig CDMA-Parameter!');
  end
  
  if (length(fieldnames(sim))<4)
      error('sc_cdma: Zu wenig Simulationsparameter!');
  end

  
  
% Speicher fuer Fehlerraten  
  ber = zeros(length(sim.EbN0dB),1);
  


% Starte Simulation ueber alle SNR
  for snr = 1:length(sim.EbN0dB)
      EbN0 = 10^(sim.EbN0dB(snr)/10);
      EsN0 = EbN0;
      
      bit_cntr = 0;
      err_cntr = 0;
      
      while ((bit_cntr<sim.max_info) & (err_cntr<sim.max_err))
         u = sign(randn(cdma.N_user,sim.N_info+1));     % Informationsbit
         if (cdma.ch_est>0)
             u = [u; ones(1,sim.N_info+1)];     % Pilotkanal ergaenzen
         end
         
         x = cdma.codes * u;                    % DS-Spreizung
          
% L-Pfad-Rayleigh-Kanal         
         [y,ch] = ant_channel(x,channel);       % Uebertragung ueber den Kanal
         
% AWGN-Kanal
         awgn = (randn(size(y)) + j*randn(size(y))) / sqrt(2*EsN0);
         y    = y + awgn;


% Kanalschaetzung
         if (cdma.ch_est>0)
            ch_hat = zeros(cdma.L_rake,sim.N_info);  
            for l=1:cdma.L_rake
                tmp = reshape(y(l:cdma.G_p*sim.N_info+l-1),cdma.G_p,sim.N_info);
                tmp = cdma.codes(:,cdma.N_user+1).'*tmp;
              % Mittelung ueber geschaetzte Werte
                tmp = convn(tmp,ones(1,cdma.ch_est),'valid');        
                ch_hat(l,:) = [tmp repmat(tmp(sim.N_info-cdma.ch_est+1),...
                               1,cdma.ch_est-1)]/cdma.ch_est;
            end
            
% Pilotsignal vom Empfangssignal abziehen (Interferenzreduktion)
            s_pilot = repmat(cdma.codes(:,cdma.N_user+1),1,sim.N_info);
            r_pilot = zeros(size(y));
            for l=1:cdma.L_rake
                tmp = s_pilot.*repmat(ch_hat(l,:),cdma.G_p,1);
                tmp = tmp(:);
                r_pilot(l:sim.N_info*cdma.G_p+l-1) = ...
                    r_pilot(l:sim.N_info*cdma.G_p+l-1) + tmp;
            end
            y = y - r_pilot;
         else
            ch_hat = repmat(ch.impulse_response,1,sim.N_info); 
         end

         

% Rake-Empfaenger
         r = zeros(cdma.N_user,sim.N_info);  
         for l=1:cdma.L_rake
             tmp = reshape(y(l:cdma.G_p*sim.N_info+l-1),...
                 cdma.G_p,sim.N_info);
             r = r + real((cdma.codes(:,1:cdma.N_user).'*tmp)...
                 .*conj(repmat(ch_hat(l,:),cdma.N_user,1)));
         end
         
         u_hat    = sign(r);
         
         % Bitfehler zaehlen         
         err      = find(u_hat~=u(1:cdma.N_user,1:sim.N_info));
         err_cntr = err_cntr + length(err);
         bit_cntr = bit_cntr + sim.N_info*cdma.N_user;

     end % while

     ber(snr) = err_cntr / bit_cntr;
     
  end % for   

% ### EOF ######################################################################
