% m-file defines integrand that has to be numerically integrated in order to
% calculate channel capacity of AWGN channel for a given input alphabet
%  only 1-dimensional integration for real valued symbols, e.g. ASK) 
%
% -----------------------------------------------------------------------------
% syntax: integrand = integrand_awgn_1d(re,alphabet,P_x,sigma2,show)
%
% -----------------------------------------------------------------------------
% Input Parameters:
%       re        : abscissa over which is integrated (equispaced)
%       alphabet  : discrete alphabet at AWGN channel input
%       P_x       : a priori probabilities of input alphabet
%       sigma2    : variance of AWGN noise 
%       show      : if show==1, integrand is plotted (default=0)
%
% -----------------------------------------------------------------------------
% Output Parameter:
%       integrand : represents function over which can be integrated
% -----------------------------------------------------------------------------
% Author: Volker Kuehn, 06.12.01
% -----------------------------------------------------------------------------
function integrand = integrand_awgn_1d(re,alphabet,P_x,sigma2,show)

if (nargin<5)
    show = 0;
end

num_sym = length(alphabet);

re  = re(:);
dre = abs(re(2)-re(1));

p_y_d = zeros(length(re),num_sym);
p_y   = zeros(length(re),1);
  
% conditional probability density functions p(y|d)
for symbol=1:num_sym
    p_y_d(:,symbol) = exp(-((re-alphabet(symbol)).^2) /2/sigma2) / sqrt(2*pi*sigma2) * dre;
    p_y = p_y + p_y_d(:,symbol)*P_x(symbol);  % probability density functions p(y)
end

% building function to be integrated
integrand = zeros(length(re),1);
for symbol=1:num_sym
    integrand = integrand + p_y_d(:,symbol)*P_x(symbol) .* log2(p_y_d(:,symbol)./p_y);
end
  

  
if show
   figure
   plot(re,integrand);
   xlabel('x')
   title('integrand for computing capacity')
end