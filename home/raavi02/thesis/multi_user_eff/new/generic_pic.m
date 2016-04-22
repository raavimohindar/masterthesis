	%clear all;
	%close all;
	
	EsN0_dB=8
	EsN0=10.^(EsN0_dB/10);
	decoder='max_log_map';
	alph_a = (31+16)/16;
	
	filename=['MU_vs_SINR_CC_' decoder '_75.mat']

	load(filename);	
	
	nu1=[0:0.01:1];

	for lauf=1:length(nu1)
	    mu_(lauf)=interp1(SINR*2,mu1,nu1(lauf)*2*EsN0);
	    nu_(lauf)=1/(1+alph_a*EsN0*mu_(lauf));
	end

	plot(nu1,nu_,'-o','markersize',3)
	xlabel('x')
	ylabel('y')
	title('t');
	set(gcf,'paperposition',[0 0 4.5 3.0]);
	set(gca,'fontname','times','fontsize',10);
	axis([0 1 0 1]);
	hold on;
	
	
	
	set(gcf,'paperposition',[0 0 4.5 3])
	print -deps generic_3_me_psfrag.eps
	
	save gnereric_5