	clear all
	SINR_dB=-10:1:30;
	SINR=10.^(SINR_dB/10);
	decoder='max_log_map';
	N=1e4
	g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
	[n,Lc] = size(g);
	d=(1-sign(randn(1,N/2-Lc+1)))/2;
	trellis=make_trellis(g,0);
	code.trellis_out=trellis.out;
	code.trellis_next=trellis.next;
	code.term=1;
	code.word_len=n;
	code.block_len=N/n;
	code.num_state=2^(Lc-1);
	c=1-2*conv_encoder(d,g,0,1);

	for s=1:length(SINR)
	    sigma=1/SINR(s);
	    rausch=sqrt(sigma)*randn(size(c));
	    x=c+rausch;
	    signal.sig=4*x*SINR(s);
	    signal.L_a=zeros(size(x));
	    signal.last_state=0;
	    [L_info(:,s),L_code(:,s)]=eval([decoder '(signal,code);']);
	    L_e(:,s)=L_code(:,s)-signal.sig;
	   
            mu(s)=mean( abs( c-tanh(L_e(:,s)/2) ).^2 );
	end
	mu1=mu;
	filename=['MU_vs_SINR_CC_' decoder '_' num2str(bi2de(g(1,:))) num2str(bi2de(g(2,:))) '.mat']
	save(filename,'SINR_dB','SINR','mu1')
