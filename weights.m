function [w] = weights(nu,N)
% the weights used in l^1_nu

w = 2*nu.^(0:N);
w(1) = 1;

end