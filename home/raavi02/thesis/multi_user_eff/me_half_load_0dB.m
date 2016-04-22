	clear all;
	close all;
	
	EsN0_dB=0
	EsN0=10.^(EsN0_dB/10);
	decoder='max_log_map';
	user=16;
	it_max=20;
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
	sf=32;
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
	tx0=data1.*pn_code;
	tx1=sum(tx0,2);
	sigma2   = sf/EsN0;
	noise=sqrt(sigma2)*randn(size(tx1));
	tx=tx1+noise;
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
	    mu(it,:)=mean(abs((data_c-tanh((L_a)/2))).^2);%mean(abs(abs(out(:))-1).^2);
	    for n=1:user
	        ind=[1:n-1 n+1:user];
	        nu(it,n)=1/(1+EsN0/sf*sum(mu(it,ind)));
	    end
	    for n=1:user
	        signal.sig(perm(n,:))=out(n,:);
	        signal.L_a=zeros(size(out(n,:)));
	        signal.last_state=0;
	        [L_info(:,n),L_code(:,n)]=eval([decoder '(signal,code);']);
	        %[L_info(:,n),L_code(:,n)]=log_map(signal,code);
	        L_a(:,n)=L_code(perm(n,:),n);
	        error(it,n)=length(find((1-data_u(:,n)*2)~=sign(L_info(1:length(data_u(:,n)),n))));
	    end
	end
	figure;
	nu_m=mean(nu.')
	test=repmat(nu_m,2,1);
	test=test(:);
	x_ax=[0;test(1:end-1)].';
	test2=repmat(nu_m,2,1);
	y_ax=test2(:).';
	plot(x_ax.',y_ax.')
	axis([0 1 0 1])
	xlabel('x')
	ylabel('y')
hold on;
title('t');
	
set(gcf,'paperposition',[0 0 4.5 3.0]);
set(gca,'fontname','times','fontsize',10);

	ME_plot_pic(EsN0,user,sf,decoder);
	hold off

print -deps me_half_load_pic_0dB_psfrag.eps
	filename=['ME_traj_' decoder '_EsN0' num2str(EsN0_dB) 'dB_user' num2str(user) '_sf' num2str(sf) '_N' num2str(N)];
	save(filename,'nu_m','EsN0','user','sf','N');
