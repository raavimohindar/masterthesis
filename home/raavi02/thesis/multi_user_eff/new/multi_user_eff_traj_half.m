	clear all;
%	close all;
	EsN0_dB=11
	EsN0=10.^(EsN0_dB/10);
	user=8;
	decoder = 'max_log_map';
	sf=8;
	sigma2   = sf/EsN0;

	
	figure(1);

	filename=['MU_vs_SINR_CC_' decoder '_75.mat']

	load(filename);	
	
	nu1=[0:0.01:1];
	
	alph_a = (user-1)/sf

	for lauf=1:length(nu1)
	    mu_(lauf)=interp1(SINR*2,mu1,nu1(lauf)*2*EsN0);
	    nu_(lauf)=1/(1+alph_a*2*EsN0*mu_(lauf));
	end

	plot(nu1,nu_,'-s','markersize',5);
	
	hold on
	
	axis([0 1 0 1]);
	
	
