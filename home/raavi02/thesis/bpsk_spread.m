%*********************BINARY PHASE SHIFT KEYING *******************
%
% Simulation of BPSK system in presence of Additive White Guassian Noise
% and Rayleigh fading and Dopper shift and using SPREAD SPECTRUM TECH
%
%******************DIRECT SEQUENCE SPREAD SEQUENCE******************

f_data=1;          % DATA FREQUENCY
f_chip=11;         % LENGTH OF CHIP SEQUENCE
fc=220;            % RELATIVE CARRIER FREQUENCY
fs=fc*3;           % SAMPLING FREQUENCY
N=fs/f_chip;       % CODING RATE
data_length=1000;  
M=2;               % BINARY LEVEL CODING

matches=0;
errors=0;
count=1;

SNRpbit=0;
SNR=SNRpbit;
rand('state',12345);    % INITIALISATION 'SEED'
randn('state',54321);

numplot=100;

msg_unspread=randsrc(data_length,1,[0:M-1]); % GENERATION OF RANDOM DATA 


PN_sequence=randsrc(11,1,[0:1]);             % GENERATION OF PN SEQUENCE

% GENERATION OF SPREADED MESSAGE
j=1;
for i = 1:data_length
    for k = j:j+f_chip
        msg_orig(k) = msg_unspread(i);
    end;
    msg_orig(j:(j+f_chip-1)) = xor(msg_orig(j:(j+f_chip-1))',PN_sequence(1:f_chip));
    j = f_chip*i+1;    
end;

len_of_orig=length(msg_orig);


% MODULATING THE SPREAD MESSAGE 
msg_tx=dmod(msg_orig,fc,f_chip,fs,'psk',M);
len_of_tx=length(msg_tx);


% % RAYLEIGH COEFFICIENTS
% magT=abs(T);
% 
% % ADDING RAYLEGH FADING COEFFICIENTS
%    k=1;
%    for i=1:len_of_tx
%       msg_tx(1,i)=msg_tx(1,i)*10*magT(1,k);
%       if(mod(i,k))
%          k=k+300;
%       end ;
%    end;    

% PERFORMANCE ANALYSIS FOR VARYING SNR'S

for SNRpbit=0:1:7
  
    SNR=SNRpbit;
    rand('state',12345);
    randn('state',54321);

   % ADDING AWGN TO THE SIGNAL 
    msg_rx_data=awgn(msg_tx,SNR-10*log10(N),'measured','dB'); 
    
    % plot(t,msg_rx(1:length(t)),'b-');
    msg_rx= msg_rx_data ;
       
   % DEMODULATING THE RECEIVED SIGNAL
    msg_demod=ddemod(msg_rx,fc,f_chip,fs,'psk',M);

   % CORRELATING WITH THE PN SEQUENCE 
    j=1;
    for i = 1:data_length
      msg_demod(j:(j+f_chip-1)) = xor(msg_demod(j:(j+f_chip-1))',PN_sequence(1:f_chip));
      j = f_chip*i+1;    
    end;

   % DESPREADING THE RECEIVED SIGNAL
    j=1;
    for i = 1:data_length
      sum=0;
      for k = j:j+f_chip
        sum=sum+msg_demod(k);
      end;
      if (sum >=7)
        msg_demod_rec(i)=1;
      else
        msg_demod_rec(i)=0;
      end;    
    j = f_chip*i+1;    
    end;
 
   % CALCULATION OF ERRORS
 
   for i=1:data_length;
     if (msg_demod_rec(1,i)== msg_unspread(i,1))
           matches=matches+1;
     else
           errors=errors+1;
     end ;
   end;     
 
  % BER_ray(count)=errors/data_length;
  BER_awgn(count)=errors/data_length;
  %BER_theo(count)=0.07*erfc(SNR^0.5);
    count=count+1;
    errors=0;
    
end;       