clear all;
close all;
EsN0_dB=3
EsN0=10.^(EsN0_dB/10);
user=4;
it_max=10;
g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
[n,Lc] = size(g);
N_u=25000-Lc+1;
N    = (N_u+Lc-1)*n;                          % number of coded bits
N_samples=200;
trellis=make_trellis(g,0);
code.trellis_out=trellis.out
code.trellis_next=trellis.next;
code.term=1;
code.word_len=n;
code.block_len=N/n;
code.num_state=2^(Lc-1);
sf=8;
I_a   = [0.001 0.05:0.05:0.95 0.999].'; %ten Brink's spacing
sigma2_a = [0.001 0.01:0.2:50].';
sigma_a  = sqrt(sigma2_a);
tmp      = zeros(length(sigma2_a),1);
for run2=1:length(sigma2_a)
    delta = 15*sigma_a(run2)/N_samples;
    x     = -5*sigma_a(run2):delta:10*sigma_a(run2);  
    tmp(run2) = sum(exp(-(x-sigma2_a(run2)/2).^2/2/sigma2_a(run2)).*(1-log2(1+exp(-x))))...
        / sqrt(2*pi*sigma2_a(run2)) * delta;
end
sigma2_a = interp1(tmp,sigma2_a,I_a);
sigma_a  = sqrt(sigma2_a);
pn_code=sign(randn(sf*N,user));
data_u=(1-sign(randn(N_u,user)))/2;
%(1-data_u(1:10,1).')/2
for n=1:user
    data(:,n)=conv_encoder(data_u(:,n),g,0,1);
end
data=1-2*data;
for n=1:user
    perm(n,:)=randperm(N);
    data_c(:,n)=data(perm(n,:),n);
end
for n=1:user
    test=repmat(data_c(:,n).',sf,1);
    data1(:,n)=test(:);
end
tx=data1.*pn_code;
tx=sum(tx,2);
sigma2   = sf/EsN0;
noise=sqrt(sigma2/2)*randn(size(tx));
tx=tx+noise;
for l=1:N
    mf(:,l)=pn_code((l-1)*sf+1:l*sf,:)'*tx((l-1)*sf+1:l*sf)/sf;
end
out=zeros(user,N);
%--------------------------------------------------------------------------
%iterative Detektion
%--------------------------------------------------------------------------
L_a=zeros(size(data_c));
for it=1:it_max
    %for na=length(sigma_a)
        %L_a=sigma2_a(na)*data_c/2+sigma_a(na)*randn(size(data_c));
        for l=1:N
            R_xx=pn_code((l-1)*sf+1:l*sf,:)'*pn_code((l-1)*sf+1:l*sf,:)/sf;
            for n=1:user
                ind=[1:n-1,n+1:user];
                out(n,l)=(mf(n,l)-R_xx(n,ind)*tanh((L_a(l,ind).')/2));
            end
        end
	sigma_eff(:,it)=var(out.'-data_c);
 	b_tilde=out;
        N_bins=40;
        mm    = max(max(abs(b_tilde)));
        delta = 2*mm/N_bins;
        bins  = -mm:delta:mm;
        number=zeros(1,length(bins));
        number1=zeros(1,length(bins));
        for n=1:user
            [a1,b1]=find(data_c(:,n)>0);
            [a,b]=find(data_c(:,n)<0);
            [numb,value]=hist(b_tilde(n,a),bins);
            number(n,:)=numb/sum(numb);
            [numb1,value1]=hist(b_tilde(n,a1),bins);
            number1(n,:)=numb1/sum(numb1);
        end
        for n=1:user
            summe(n)=0;
            idx=find(((number(n,:)+number1(n,:))~=0));
            idx1=find(((number(n,:)+number1(n,:))~=0));
            I_e_ic(n,:)=full(0.5*sum(number(n,idx).*spfun(@log2,2*number(n,idx)./(number(n,idx)+number1(n,idx))))+...
                0.5*sum(number1(n,idx1).*spfun(@log2,2*number1(n,idx1)./(number(n,idx1)+number1(n,idx1)))));
        end
        %end

        for n=1:user
            signal.sig(perm(n,:))=out(n,:)*4/sigma_eff(n,it);
            signal.L_a=zeros(size(out(n,:)));
            signal.last_state=0;
            [L_info(:,n),L_code(:,n)]=max_log_map(signal,code);
            L_a(:,n)= L_code(perm(n,:),n)-(out(n,:).'*4/sigma_eff(n,it));
        end
        b_tilde=L_a.';%-L_a.';
        N_bins=40;
        mm=max(max(abs(b_tilde)));
        delta=2*mm/N_bins;
        bins=-mm:delta:mm;
        number=zeros(1,length(bins));
        number1=zeros(1,length(bins));
        for n=1:user
            [a1,b1]=find(data_c(:,n)>0);
            [a,b]=find(data_c(:,n)<0);
            [numb,value]=hist(b_tilde(n,a),bins);
            number(n,:)=numb/sum(numb);
            [numb1,value1]=hist(b_tilde(n,a1),bins);
            number1(n,:)=numb1/sum(numb1);
        end
        for n=1:user
            summe(n)=0;
            idx=find(((number(n,:)+number1(n,:))~=0));
            idx1=find(((number(n,:)+number1(n,:))~=0));
            I_e_cc_(n,:)=full(0.5*sum(number(n,idx).*spfun(@log2,2*number(n,idx)./(number(n,idx)+number1(n,idx))))+...
                0.5*sum(number1(n,idx1).*spfun(@log2,2*number1(n,idx1)./(number(n,idx1)+number1(n,idx1)))));
        end
%--------------------------------------------------------------------------
       
    I_e_cc_m(it)=mean(I_e_cc_);
    I_e_ic_m(it)=mean(I_e_ic);
end
figure(1);
test=repmat(I_e_cc_m,2,1);
x_ax=[0;test(:)];
x_ax(end)=[];
test2=repmat(I_e_ic_m,2,1);
y_ax=test2(:);
plot(x_ax,y_ax,'-x');
hold on;
plot(x_ax,y_ax);
legend('aaaaaaaaaa','bbbbbbbbbb','location','northeastoutside');

load('exit_chart_3dB_4Nutzer_sf8_200000bit_PIC.mat')
plot(I_e_cc,I_a,I_a,I_e_pic)
hold off;
axis([0 1 0 1])
grid
xlabel('x');
ylabel('y');
title('t');
set(gcf,'paperposition',[0 0 4.5 3]);
set(gca,'fontname','times','fontsize',10);
print -deps exit_trajectories_3dB_half_psfrag.eps
%filename=['exit_chart_' num2str(EsN0_dB) 'dB_' num2str(user) 'Nutzer_sf' num2str(sf) '_' num2str(N) 'bit_PIC']
%save(filename)
