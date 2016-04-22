close all;
clear all;

%  ------------------------  uplink design.     


% -------------------------  number of users.  

N_u = 4;



% ------------------------- define number of iterations.

it_max = 10;



% ---- define the spreading.

sf = 8;

% -------------------------- initilize a variable for spreaded sequence.	

x1 = [];

% -------------------------  parameters for coding  


% -----  defining the generator polynomial. 

%G = [7;5];                % [1 1 1; 1 0 1]
%G = [15;11];               % [1 1 1 1; 1 1 0 1]
%G = [19;31];              % [1 0 0 1 1; 1 1 1 1 1] 
%G = [109;79];             % [1 0 1 1 0 1 1; 1 1 1 1 0 0 1];

g = [1 1 1 1; 1 1 0 1];

% -----  additional prameters for the encoder.

r_flag = 1;                % flag indicating the feedback polynomial for RSC codes
% r_flag=0: NSC code
% r_flag=i: i-th polynomial used for feedback
term = 0;                  % flag indicating trellis truncation or termination
% term=0: truncation, term=1: termination

[n,Lc] = size(g);

% -------------------------  number of bits that each user can transmit. 

N = 10;

% -------------------------  number of coded bits.

N_c = N*n;

% ----  operating rate. 

Rc = 1/2;

% ---- spreading sequence.

pn_seq = sign(randn(sf*N_c,N_u));


% ----  Date for the APP Decoder.

trellis = make_trellis(g,r_flag);
code.trellis_out = trellis.out;
code.trellis_next = trellis.next;
code.num_state = 2^(Lc-1);
code.block_len = N+(Lc-1)*term;
code.word_len = n;
code.term = term;


% -------------------------  Start of the transmission. 

for i=1:N_u
    
    % -------------------------  info bits.
    
    d(:,i) = (1-sign(randn(N,1)))/2;
    
    [c(:,i),ls,tb] = conv_encoder(d(:,i),g,r_flag,term); 
    
    
    % BPSK
    
    c(:,i) = 1-2*c(:,i);
    
    % -------------------------  random interleaving.
    
    rand_int(:,i) = randperm(length(c))';
    
    inter_leave(:,i) = c(rand_int(:,i));
    
    %..........................  spreading. 	
    
    test = repmat(inter_leave(:,i).',sf,1);
    
    data_spread(:,i) = test(:);
    
    % ---- generate the spreading sequence. Randome spreading sequence is used.
    
    
    
end	

% ....  summing the spreaded sequence before transmitting.

spread = data_spread.*pn_seq;

spread_sum = sum(spread,2);	

% ------------------------- end of transmission. 


% ------------------------- air interface.

EbN0dB   = 10;
EsN0     = 10.^(EbN0dB/10)*Rc;
sigma2_N = 1 ./ EsN0;
sigma_N  = sqrt(sigma2_N/2);

noise = randn(size(spread_sum)) * sigma2_N; % noise

%--------------------------- end of air interface.


% ------------ downlink.


% ----  received matrix of y.

y = spread_sum+noise;	% of the format y = Ad + n

% ----  matched filter.

for u = 1:N_u
    
    for j = 1:length(c)
        
        mf(j,u) = pn_seq((j-1)*sf + [1:sf],u)' * y((j-1)*sf + [1:sf])/sf;
        
    end
    
end


L_a=zeros(size(c));

% .... number of iterations to carry on.

for iter = 1:it_max
    
    for k = 1:length(c)
        
        rxx=pn_seq((k-1)*sf + 1:k*sf,:)'*pn_seq((k-1)*sf + 1:k*sf,:)/(sf);
        
        for n = 1:N_u
            
                       
            idx = [1:n-1,n+1:N_u];
            
            out(k,n) = (mf(k,n)-rxx(n,idx)*tanh((L_a(k,idx).')/2));
            
            
        end
        
    end
    
     sigma_eff(iter,:) = var(out-inter_leave);  % sigma_eff --> after first cancellation.
    
    for n = 1:N_u
    
        signal.sig(rand_int(:,n)) = out(:,n);
        
        signal.L_a = zeros(size(out));
        
        signal.last_state = 0;
        
        [L_info(:,n),L_code(:,n)]=log_map(signal,code); 
        
        L_a(:,n)=L_code(rand_int(:,n),n);
                
    end
    
    sigma_d(iter,:) = var(tanh(L_a/2)-inter_leave);
    
end

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
alpha=(N_u-1)/sf;
sigma_eff_ic=sigma2_N/(2*sf)+alpha*sigma_d_ic;
plot(sigma_d_ic,sigma_eff_ic,'r')
filename=['sigma_d_vs_sigma_eff' num2str(bi2de(g(1,:))) num2str(bi2de(g(2,:))) '.mat'];
load(filename)
plot(sigma_d_cc,sigma_eff_cc,'k')
hold off
axis
xlabel('x');
ylabel('y');

axis([0 1 -1 4])
%set(gca,'xtick',[0 1/3 1]);
%set(gca,'ytick',[0  1]);

set(gca,'fontname','times','fontsize',10);
set(gcf,'paperposition',[0 0 4.5 3.0]);

print -deps tc_trans_pic_psfrag.eps
