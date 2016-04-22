function dummy=ME_plot_pic(EsN0,user,sf,decoder,alpha)

filename=['MU_vs_SINR_CC_' decoder '_75.mat']

load(filename);
alpha=(user-1)/sf;
nu1=[0:0.01:1];

for lauf=1:length(nu1)
    mu_=interp1(SINR,mu1,nu1(lauf)*EsN0);
    nu_(lauf)=1/(1+alpha*EsN0*mu_);
end

plot(nu1,nu_,'-k',[0,1],[0,1],'-g')
