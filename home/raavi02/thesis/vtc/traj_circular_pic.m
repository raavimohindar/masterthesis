clear all
close all
set(gca,'fontname','times','fontsize',10);
%set(gcf,'paperposition',[0 0 4.5 3.0]);

EsN0_dB=-3
EsN0=10.^(EsN0_dB/10);
user=4;
it_max=24;
g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
[n,Lc] = size(g);
N_u=1000-Lc+1;
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
        for l=1:N
            R_xx=pn_code((l-1)*sf+1:l*sf,:)'*pn_code((l-1)*sf+1:l*sf,:)/sf;
            for n=1:user
                ind=[1:n-1,n+1:user];
                out(n,l)=(mf(n,l)-R_xx(n,ind)*tanh((L_a(l,ind).')/2));
            end
        end
        
	sigma_eff(it,:)=var(out.'-data_c);
     
        for n=1:user
	    signal.sig(perm(n,:))=out(n,:);
            signal.L_a=zeros(size(out(n,:)));
            signal.last_state=0;
            [L_info(:,n),L_code(:,n)]=max_log_map(signal,code);
	    L_a(:,n)=L_code(perm(n,:),n);

	end
if it == 1 


figure(it);

subplot(2,2,1)
plot(abs(tanh(L_a(:,1)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);


subplot(2,2,2)
plot(abs(tanh(L_a(:,2)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,3)
plot(abs(tanh(L_a(:,3)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,4)
plot(abs(tanh(L_a(:,4)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

print -deps a0_psfrag.eps

elseif it == abs(it_max/2)

figure(it);

subplot(2,2,1)
plot(abs(tanh(L_a(:,1)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,2)
plot(abs(tanh(L_a(:,2)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,3)
plot(abs(tanh(L_a(:,3)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,4)
plot(abs(tanh(L_a(:,4)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

print -deps b0_psfrag.eps

elseif it == it_max


figure(it);

subplot(2,2,1)
plot(abs(tanh(L_a(:,1)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,2)
plot(abs(tanh(L_a(:,2)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,3)
plot(abs(tanh(L_a(:,3)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

subplot(2,2,4)
plot(abs(tanh(L_a(:,4)/2).^2),'.');
set(gca,'ytick',[0 1]);
set(gca,'xtick',[0 2000]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

print -deps c0_psfrag.eps

else


end

   end
