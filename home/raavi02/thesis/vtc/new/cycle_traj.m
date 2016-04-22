clear all
close all
randn('state',sum(100*clock));
EsN0_dB=1
EsN0=10.^(EsN0_dB/10);
user=8;
it_max=10;
g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
[n,Lc] = size(g);
N_u=50000-Lc+1;
N    = (N_u+Lc-1)*n;                          % number of coded bits
N_samples=200;
trellis=make_trellis(g,0);
code.trellis_out=trellis.out;
code.trellis_next=trellis.next;
code.term=1;
code.word_len=n;
code.block_len=N/n;
code.num_state=2^(Lc-1);
sf=8;
alpha=(user-1)/sf;
    pn_code=sign(randn(sf*N,user));
    data_u=(1-sign(randn(N_u,user)))/2;
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
    it
        for l=1:N
            R_xx=pn_code((l-1)*sf+1:l*sf,:)'*pn_code((l-1)*sf+1:l*sf,:)/sf;
            for n=1:user
                ind=[1:n-1,n+1:user];
                out(n,l)=(mf(n,l)-R_xx(n,ind)*tanh((L_a(l,ind).')/2));
            end
        end
        
	sigma_eff(it,:)=var(out.'-data_c);
        sigma_e_test(it,:)=mean(abs(out.'-data_c).^3);
     

        for n=1:user
	    signal.sig(perm(n,:))=out(n,:);
            signal.L_a=zeros(size(out(n,:)));
            signal.last_state=0;
            [L_info(:,n),L_code(:,n)]=log_map(signal,code);
	    
	     L_a(:,n)=L_code(perm(n,:),n);%-(out(n,:).'*4/sigma_eff(it,n));
	
            %error(it,n)=length(find((1-data_u(:,n)*2)~=sign(L_info(1:length(data_u(:,n)),n))));
        end
        sigma_d(it,:)=var(tanh(L_a/2)-data_c);
	
        sigma_d_test(it,:)=mean(abs((tanh(L_a/2)-data_c)).^3);
    end
    s_d(2:it_max+1)=mean(sigma_d.');
    s_d(1)=1;
    s_e=mean(sigma_eff.');

figure(1);
test=repmat(s_d,2,1);
x_ax=test(:);
x_ax(1)=[];
x_ax(end)=[];
test=repmat(s_e,2,1);
y_ax=test(:);
plot(x_ax,y_ax,'-x');
hold on;
plot(x_ax,y_ax);
legend('aaaaaaaaaaaa','bbbbbbbbbbbb','location','northeastoutside');
sigma_d_ic=0:0.1:1;
alpha=(user-1)/sf;
sigma_eff_ic=sigma2/(2*sf)+alpha*sigma_d_ic;
plot(sigma_d_ic,sigma_eff_ic);
hold on;
filename=['sigma_d_vs_sigma_eff' num2str(bi2de(g(1,:))) num2str(bi2de(g(2,:))) '.mat'];
load(filename)
plot(sigma_d_cc,sigma_eff_cc);
hold on;
axis([0 1 0 2]);
xlabel('x')
ylabel('y')
title('t')
grid on;
set(gca,'fontname','times','fontsize',10);
set(gcf,'paperposition',[0 0 4.5 3]);

print -deps cycle_psfrag.eps
