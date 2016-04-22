% m-file for calculating channel capacity of AWGN channel for BPSK

  clear all

  prnt = 0;

  mode = input('new calculation [y/n]?  ','s');
  
  M = 2;
  
if strcmp(lower(mode),'y')
  
% building input alphabet
  alphabet = [-1; 1];
  
% equiprobable symbols in alphabet  
  P_x = ones(M,1)/M;
  
% signal-to-noise-ratios on AWGN channel
  snr = [0:20].';
  awgn.bpsk.q0.snr = snr;
  
% computing capacity for each SNR
  awgn.bpsk.q0.C = zeros(size(snr));
  
  for run = 1:length(snr)
      awgn.bpsk.q0.C(run)  = c_awgn(alphabet,P_x,snr(run));
  
  end
% computing Eb/N0 dependant on snr, C and R0  
  awgn.bpsk.q0.EbN0_C  = awgn.bpsk.q0.snr - 10*log10(awgn.bpsk.q0.C);
  
 save Results/awgn_bpsk awgn
  
else   % load old results
  load Results/awgn_bpsk
end
  
  figure(1)
  plot(awgn.bpsk.q0.snr,awgn.bpsk.q0.C,'-o')
  grid on
  axis([0 10 0 1])
  title('capacity of AWGN and BPSK')
  xlabel('SNR in dB')
  ylabel('C in bit/s/Hz')
  legend('q=inf','q=1','q=2','q=3','q=4',4);
  if prnt
     print -depsc c_awgn_bpsk_EsN0  
  end
  
  