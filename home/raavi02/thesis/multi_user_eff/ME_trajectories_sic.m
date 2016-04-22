clear all
EsN0_dB=10
EsN0=10.^(EsN0_dB/10);
decoder='max_log_map';
user=32;
it_max=20;
g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
[n,Lc] = size(g);
N_u=100-Lc+1;
N    = (N_u+Lc-1)*n;                          % number of coded bits
N_samples=200;
trellis=make_trellis(g,0);
code.trellis_out=trellis.out;
code.trellis_next=trellis.next;
code.term=1;
code.word_len=n;
code.block_len=N/n;
code.num_state=2^(Lc-1);
sf=16;
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
mu_(1,1:user)=alpha;
for it=1:it_max
    for n=1:user
        for l=1:N
            R_xx=pn_code((l-1)*sf+1:l*sf,:)'*pn_code((l-1)*sf+1:l*sf,:)/sf;
            ind=[1:n-1,n+1:user];
            out(n,l)=(mf(n,l)-R_xx(n,ind)*tanh((L_a(l,ind).')/2));
        end
        signal.sig(perm(n,:))=out(n,:);
        signal.L_a=zeros(size(out(n,:)));
        signal.last_state=0;
        [L_info(:,n),L_code(:,n)]=eval([decoder '(signal,code);']);        
        %[L_info(:,n),L_code(:,n)]=max_log_map(signal,code);
        error(it,n)=length(find((1-data_u(:,n)*2)~=sign(L_info(1:length(data_u(:,n)),n))));
        L_a(:,n)=L_code(perm(n,:),n);
        mu_(it,n)=mean(abs((data_c(:,n)-tanh((L_a(:,n))/2))).^2);
        for u=1:user
            ind=[1:u-1 u+1:user];
            nu(it,u)=1/(1+2*EsN0/sf*sum(mu_(it,ind)));
        end
    end
end
for n=1:user
    test=repmat(nu(:,n).',2,1);
    test=test(:);
    x_ax(n,:)=[0;test(1:end-1)].';
    test2=repmat(nu(:,n).',2,1);
    y_ax(n,:)=test2(:).';
end
plot(x_ax.',y_ax.')
axis([0 1 0 1])
hold on
ME_plot_pic(EsN0,user,sf,decoder)
hold off