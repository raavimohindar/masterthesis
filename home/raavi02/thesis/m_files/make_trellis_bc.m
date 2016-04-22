%-*- matlab -*-
%
%
% m-file for generating a trellis from the parity check matrix H
%
%
% -------------------------- input parameters -----------------------------
% H: parity check matrix of n-k rows and n columns containing 0's and 1's
% p: indicates the column to be neglected by generating the trellis 
%    (necessary for decoding with extended Battail algorithm)
%    if p==0, no column is omitted
%
% -------------------------- output parameters -----------------------------
% trellis: array with 2^(n-k) rows and n columns (n-1 columns for p>0)
%          row i represents current state i, 
%          contents of row i and column j describes next state at time j+1 for info=1
%          (info = 0 remains in the same state)
%
function trellis = make_trellis_bc(H,p)


[m,n] = size(H);   % m=n-k

state = 0:2^m-1;                        % vector containing states in decimal representation
state_bin = de2bi(state,m);
%state_bin = de2bi(state,m,'left-msb');  % matrix containing all states in binary representation
                                        % each row contains one state
dummy = ones(2^m,1);                    % vector for expanding rows of H
  
if (p == 0)
   trellis = zeros(2^m,n);
else
   trellis = zeros(2^m,n-1);
   H(:,p) = [];                    % discard p-th column
   n = n-1;
end

for i=1:n
    H_i = dummy * H(:,i)';
    trellis(:,i) = bi2de( rem( state_bin + H_i, 2 ) );
end
