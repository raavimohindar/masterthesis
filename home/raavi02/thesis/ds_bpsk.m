%DS BPSK Transmitter
%Run from editor debug(F5).
%m-file for generating PRBS(Pseudo Random Binary Sequence) which modulates a
%carrier to get a DS BPSK(Direct Sequence Binary Phase Shift Keyed) transmitter.(Spread Spectrum) 
%The PRBS is exclusive ored with the message and then used to modulate the carrier. The  
%file is scaleable so should be able to run much higher frequencies if
%desired. Keep message bit duration an integral multiple of the PRBS bit or
%chip duration.
fcarr=1e3;              % Carrier frequency
N = 100;		        % Number of symbols
fmess = 1e1;		    % Message frequency
fs = 4*1e3;		        % Sampling frequency
Fn = fs/2;              % Nyquist frequency
Ts = 1/fs;	            % Sampling time = 1/fs
T = 1/100;		        % Symbol time
randn('state',0);       % Keeps PRBS from changing on reruns
t = [0:Ts:N*T-Ts];      % Time vector

symbols = sign(randn(N,1))';%generate N random binary symbols, +/- 1(notice transpose')
symbols1 = ones(T/Ts,1)*symbols;	
s2 = symbols1(:); 	%PRBS @ +/-1

s2(s2>0)=1;%use to get PRBS @ +1/0 for exclusive or(modulus 2 adder) use
s2(s2<0)=0;
%==============================================================
%generate squarewave message 
%cosine wave
%2 pi fc t is written as below
twopi_fc_t=(1:fs)*2*pi*fmess/fs; 
a=1;
phi=0;
x = a * cos(twopi_fc_t + phi);

%generate square wave
x(x>0)=1;
x(x<0)=0;
%==========================================================
%plot PRBS sequence,message and x_or out
figure(1)
subplot(3,2,1)
plot(t,s2)
axis([0 max(t) -1.2 1.2])
xlabel('                                            time');
grid on
title('PRBS')
 
 
subplot(3,2,3) 
plot(t,x)
axis([0 max(t) -1.2 1.2])
grid on
title('Message')

s2=s2';%transpose for matrix match(4000x1 not equal to 1x4000)
x_or=xor(x,s2);


subplot(3,2,5) 
plot(t,x_or)
axis([0 max(t) -1.2 1.2])
grid on
title('Xor out +/0')

%generate carrier wave
%cosine wave
%2 pi fc t is written as below
twopi_fc_t=(1:fs)*2*pi*fcarr/fs; 
a=1;
phi=0;
f_c = a * cos(twopi_fc_t + phi);

x_or = x_or-mean(x_or);%remove dc component
x_or(x_or>0)=1;%use to get x_or back to +/- 1 to remove carrier
x_or(x_or<0)=-1;

subplot(3,2,2) 
plot(t,x_or)
axis([0 max(t) -1.2 1.2])
grid on
title('Xor out +/-')

f_c=x_or.*f_c;%multiply carrier and +/- x_or out

%========================================================================
%take FFT of modulated carrier
y=f_c;
NFFY=2.^(ceil(log(length(y))/log(2)));

FFTY=fft(y,NFFY);%pad with zeros
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(y);
f1=(0:NumUniquePts-1)*2*Fn/NFFY;

%plot frequency domain
subplot(3,2,4); plot(f1,MY);xlabel('');ylabel('AMPLITUDE');
axis([500 1500 -.5 .5]);%zoom in/out
title('Frequency domain plots');
grid on;
subplot(3,2,6); plot(f1,20*log10(abs(MY).^2));xlabel('FREQUENCY(Hz)');ylabel('DB');
axis([500 1500 -60 5]);
grid on;
title('Modulated DS BPSK carrier')


