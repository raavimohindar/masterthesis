%QPSK Implemenatation
clear;
%This is the first part of the problem where we are first generating a
%random input of size 6 bits, which results in 3 QPSk sysmbols;
Input = round(rand(1,6));
%Conversion of the bits to -1 and +1;
s = 2*Input - 1;

%The frequncy of the signal along with its sampling frequency which should
%be 10 times more than the frequency of the signal as recommended;
F = input('Enter the frequency of the signal ');
Fs = input('Enter the sampling frequncy ');

%Generating the in-phase component of the first symbol
I = s(1)*sqrt(2*F)*cos(2*pi*F/Fs*[0:(Fs/F)]);
%Generating the quadrature phase component of the first symbol
Q = s(2)*sqrt(2*F)*sin(2*pi*F/Fs*[0:(Fs/F)]);

%Generating the second and the third symbol of the signal;
for i =2:3
I1 = s(2*i-1)*sqrt(2*F)*cos(2*pi*F/Fs*[0:Fs/F]);
Q1 = s(2*i)*sqrt(2*F)*sin(2*pi*F/Fs*[0:Fs/F]);
I = [I I1];
Q = [Q Q1];
end

%plot for Quadrature phase component
plot(I);
title('In-Phase phase component of QPSK');
xlabel('Sampled Time Axis');
ylabel('Amplitude');
figure;

%Plot for in-phase component
plot(Q);
title('Quadrature phase component of QPSK');
xlabel('Sampled Time Axis');
ylabel('Amplitude');
figure;
%Obtaning the QPSk Signal
x = I+Q;

title('QPSK Signal');
xlabel('Samples Time Axis');
ylabel('Amplitude');
plot(x);

figure;
a2 = s(1);
for j  = 0:F/Fs:.01
    a1 = s(1);
    a2 = [a2 a1];
end

for i = 2:3
    a = s(2*i-1);
      for j = 0:F/Fs:.01
         a = [a s(2*i-1)];   
      end
    a2 = [a2 a];
end
    plot(a2);
    title('I');
    
figure;
b2 = s(2);
for j  = 0:F/Fs:.01
    b1 = s(2);
    b2 = [b2 b1];
end

for i = 2:3
    b = s(2*i);
      for j = 0:F/Fs:.01
         b = [b s(2*i)];   
      end
    b2 = [b2 b];
end
    plot(b2);
    title('Q');


 
%Detection of QPSK signals 
N=1000;%Number of Iterations
s1=[1 0];
s2=[0 1];
s3=[-1 0];
s4=[0 -1];
No=.1507;

symbolerror=0;% initializing eror 
for i=1:N
    %generation of uniformly distributed symbols
    a=rand;
    if(a<=0.25)
        data=s1;
    elseif(a<=0.5&a>0.25)
        data=s2;
    elseif(a>0.5&a<=0.75)
        data=s3;
    else
        data=s4;
    end
    
%generation of noise

%corrupting the data with noise

    n(1)=randn*sqrt(No/2); %generating noise with variance No/2
    n(2)=randn*sqrt(No/2);
% Uncomment if Rayleigh fading
    %r=raylrnd(1)*data+n;
    r=data+n;
%detection of corrupted signal
    d1=dot(r,s1);
    d2=dot(r,s2);
    d3=dot(r,s3);
    d4=dot(r,s4);
    d=max([d1 d2 d3 d4]);
    if(d==d1)
        decision=s1;
    elseif(d==d2)
        decision=s2;
    elseif(d==d3)
        decision=s3;
    else(d==d4)
        decision=s4;
    end
    if(data~=decision)
        symbolerror=symbolerror+1;
    end
    X(i,:)=r;
end
perr=symbolerror/N;
figure;
plot(X(:,1),X(:,2),'*');%plot of the received signal points
title('Received signal')
figure;
w=[s1;s2;s3;s4];
plot(w(:,1),w(:,2),'x');
title('Transmitted signal')

 
