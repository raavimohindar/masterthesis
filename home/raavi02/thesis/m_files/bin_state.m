function bin_state = bin_state( int_state, m )
% Copyright 1997 Yufei Wu
% for academic use only
% converts an vector of integer into a matrix; the i-th row is the binary form 
% of m bits for the i-th integer

for j = 1:length( int_state )
   for i = m:-1:1
       state(j,m-i+1) = fix( int_state(j)/ (2^(i-1)) );
       int_state(j) = int_state(j) - state(j,m-i+1)*2^(i-1);
   end
end

bin_state = state;

