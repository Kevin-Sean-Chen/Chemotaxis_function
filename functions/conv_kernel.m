function [F] = conv_kernel(X,K)
%%%
% using time series X and kernel K and compute convolved X with K
% this takes care of padding and the time delay
% note that K is a causal vector with zero-time in the first index
%%%
    padding = ones(1,floor(length(K)/2));  %padding vector
    Fp = conv([padding, X, padding], K, 'same');  %convolvultion of vectors
    F = Fp(1:length(X));  %return the right size
end