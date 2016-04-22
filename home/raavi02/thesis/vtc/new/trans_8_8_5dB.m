clear all
sigma_eff_cc=[0.01 0.1:0.1:5];
SNR_dB=5
SNR=10.^(SNR_dB/10);
N=2e5;
%g=[1 0 0 1 1; 1 1 0 1 1];
g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
[n,Lc] = size(g);
d=(1-sign(randn(1,N/n-Lc+1)))/2;
trellis=make_trellis(g,0);
code.trellis_out=trellis.out;
code.trellis_next=trellis.next;
code.term=1;
code.word_len=n;
code.block_len=N/n;
code.num_state=2^(Lc-1);
c=1-2*conv_encoder(d,g,0,1);
for s=1:length(sigma_eff_cc)
s
    mui=sqrt(sigma_eff_cc(s))*randn(size(c));
    x=c+mui;
    sigma_e_test(s,:)=mean(abs(mui).^3);
    signal.sig=x*4/sigma_eff_cc(s);%out(n,:);
    signal.L_a=zeros(size(x));
    signal.last_state=0;
    [L_info(:,s),L_code(:,s)]=log_map(signal,code);
    L_e(:,s) = L_code(:,s) - signal.sig;
    sigma_d_cc(s)=var(tanh(L_e(:,s)/2)-c);
    sigma_d_test(s,:)=mean(abs(tanh(L_code(:,s)/2)-c).^3);
end
figure
plot(sigma_d_test,sigma_e_test)
grid
hold on
sigma_d_ic=0:0.1:1;
alpha=7/8;
sigma_eff_ic=1/SNR+alpha*sigma_d_ic;
plot(sigma_d_ic,sigma_eff_ic,'r')
axis([0 1 0 1]);
hold off
%axis([0 1 0 4])
filename=['sigma_d_vs_sigma_eff' num2str(bi2de(g(1,:))) num2str(bi2de(g(2,:))) '.mat'];
save(filename,'SNR_dB','SNR','sigma_eff_cc','sigma_d_cc','sigma_eff_ic','sigma_d_ic')
