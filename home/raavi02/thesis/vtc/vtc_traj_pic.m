% ---- default system settings.
    clear all;
    close all;

% ---- resets it to a different state each time. 
    randn('state',sum(100*clock));
    
% ---- defining the signal to noise ratio. 
    EsN0_dB = 10;
    EsN0 = 10.^(EsN0_dB/10);

% ---- number of users.
    user = 6;
    
% ---- maximum number of iterations.    
    it_max = 20;
    
% ---- generator polynomial.    
    g = [1 1 1; 1 0 1];                             % R=1/2, L_c=3 generator polynom
    [n,Lc] = size(g);
    
% ---- number of information bits.
    N_u = 10-Lc+1;

% ---- number of coded bits.
    N = (N_u+Lc-1)*n;        
    
% ---- number of samples -- unknown parameter.    
    N_samples=200;
    
% ---- data for trellis code.
    trellis=make_trellis(g,0);
    code.trellis_out=trellis.out;
    code.trellis_next=trellis.next;
    code.term=1;
    code.word_len=n;
    code.block_len=N/n;
    code.num_state=2^(Lc-1);
    
% ---- spreading factor.    
    sf=8;

% ---- defining alpha a unknown parameter.
    alpha=(user-1)/sf;
    
% ---- PN sequence.    
    pn_code=sign(randn(sf*N,user));
    
% ---- user information bits.
     data_u=(1-sign(randn(N_u,user)))/2;
     
     
%========Convolution Encoding of information bits.     
	for n=1:user
        data(:,n)=conv_encoder(data_u(:,n),g,0,1);
	end
    
%========BPSK-Modulation.    
    data=1-2*data;
    
%========Random Interleaving.    
	for n=1:user
        perm(n,:)=randperm(N);
        data_c(:,n)=data(perm(n,:),n);
	end
    
%=======Random Spreading.    
	for n=1:user
        test=repmat(data_c(:,n).',sf,1);
        data1(:,n)=test(:);
	end

%=======Data Transmission.
	tx=data1.*pn_code;
	tx=sum(tx,2);

%=======Channel Noise.
	sigma2   = sf/EsN0;
	noise=sqrt(sigma2/2)*randn(size(tx));

%=======Transmission with Additive Noise.
    tx=tx+noise;

%=======Despreading with multiplying the transpose of the PN sequence.
	for l=1:N
        mf(:,l)=pn_code((l-1)*sf+1:l*sf,:)'*tx((l-1)*sf+1:l*sf)/sf;
	end
    
%=======Output of the Matched-Filter.    
    out=zeros(user,N);
    
%--------------------------------------------------------------------------
%===========================iterative Detektion============================
%--------------------------------------------------------------------------

    L_a=zeros(size(data_c));

    for it=1:it_max             % ---- maxmimum number of iterations. 
        it
        
        for l=1:N               % ---- number of code bits.
    
            R_xx=pn_code((l-1)*sf+1:l*sf,:)'*pn_code((l-1)*sf+1:l*sf,:)/sf; % ---- correlation props.
            
            for n=1:user        % ---- maximum number of users.
                
                ind=[1:n-1,n+1:user];
                out(n,l)=(mf(n,l)-R_xx(n,ind)*tanh((L_a(l,ind).')/2)); % ---- initial cancellation.
                
            end
            
        end
        
        sigma_eff(it,:)=var(out.'-data_c);              % ---- calculation of effective variance.
        
        sigma_e_test(it,:)=mean(abs(out.'-data_c).^3);  % ---- testing variable.
        
        out_test(it,:)=out(1,:);                        % ---- another testing variable.
        
        %==========Decoding of the matched filter output.
        
        for n=1:user                                    % ---- maximum number of users.
            
            signal.sig(perm(n,:))=out(n,:);             % ---- De-interleaving.
    
            signal.L_a=zeros(size(out(n,:)));           % ---- a-priori information bits.
    
            signal.last_state=0;                        % ---- Last state.
    
            [L_info(:,n),L_code(:,n)]=log_map(signal,code);     % ---- max_log decoding.
    
            L_a(:,n)=L_code(perm(n,:),n);               % ---- a-priori information bits.
    
                       
        end
    
        sigma_d(it,:)=var(tanh(L_a/2)-data_c);          % ---- calculation of sigma_d.
    
        sigma_d_test(it,:)=mean(abs((tanh(L_a/2)-data_c)).^3);      % ---- another testing variable.
        
end

%=======Calculation of average variance.
    s_d(2:it_max+1)=mean(sigma_d.');
    s_d(1)=1;
    s_e=mean(sigma_eff.');
    

figure;
test=repmat(s_d,2,1);
x_ax=test(:);
x_ax(1)=[];
x_ax(end)=[];
test=repmat(s_e,2,1);
y_ax=test(:);
plot(x_ax,y_ax)
grid
hold on
sigma_d_ic=0:0.1:1;
alpha=(user-1)/sf;
sigma_eff_ic=sigma2/(2*sf)+alpha*sigma_d_ic;
plot(sigma_d_ic,sigma_eff_ic,'r')
filename=['sigma_d_vs_sigma_eff' num2str(bi2de(g(1,:))) num2str(bi2de(g(2,:))) '.mat'];
load(filename)
plot(sigma_d_cc,sigma_eff_cc,'k')
hold off
xlabel('\sigma_d')
ylabel('\sigma_eff')
