	clear all;
	close all;
	EsN0_dB=2.5
	EsN0=10.^(EsN0_dB/10);
	user=8;
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
	decoder = 'max_log_map'
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
	it
	    for l=1:N
	         R_xx=pn_code((l-1)*sf+1:l*sf,:)'*pn_code((l-1)*sf+1:l*sf,:)/sf;
	        for n=1:user
	            ind=[1:n-1,n+1:user];
	            out(n,l)=(mf(n,l)-R_xx(n,ind)*tanh((L_a(l,ind).')/2));
	        end
	    end
	    
	    sigma_eff(it,:)=var(out.'-data_c);
	    
	    mu(it,:)=mean(abs((data_c-tanh((L_a)/2))).^2);
	    
	    for n=1:user
		    signal.sig(perm(n,:))=out(n,:)*4/sigma_eff(it,n);
	            signal.L_a=zeros(size(out(n,:)));
	            signal.last_state=0;
	            [L_info(:,n),L_code(:,n)]=log_map(signal,code);
		    L_a(:,n)=L_code(perm(n,:),n)-(out(n,:).'*4/sigma_eff(it,n));
            end
            
            for n=1:user
	        ind=[1:n-1, n+1:user];
	        nu(it,n)=1/(1+2*EsN0*(sum(mu(it,ind)))/sf);
	    end
	    
	end
	
	figure;

	nu_m=mean(nu.');

	test=repmat(nu_m,2,1);
	test=test(:);

	x_ax=[0;test(1:end-1)].';

	test2=repmat(nu_m,2,1);

	y_ax=test2(:).';
	
	grid on;
	
	plot(x_ax,y_ax,'-x','markersize',5);
	hold on;
	plot(x_ax,y_ax);
	
	legend('aaaaaa','bbbbbb','location','northeastoutside');
	
	hold on;
	
	filename=['MU_vs_SINR_CC_' decoder '_75.mat']

	load(filename);	
	
	nu1=[0:0.01:1];
	
	alph_a = (user-1)/sf

	for lauf=1:length(nu1)
	    mu_(lauf)=interp1(SINR*2,mu1,nu1(lauf)*2*EsN0);
	    nu_(lauf)=1/(1+alph_a*2*EsN0*mu_(lauf));
	end

	plot(nu1,nu_,[0,1],[0,1]);
	
	hold on
	
	axis([0 1 0 1]);
	
	xlabel('x');
	ylabel('y');
	title('t');
	
	set(gca,'fontname','times','fontsize',10);
	set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1]);
	set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
	
	set(gcf,'paperposition',[0 0 4.5 3.0]);
	
	print -deps mue_full_2p5dB_psfrag.eps
	
	
	
	
	
