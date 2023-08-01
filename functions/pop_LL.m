function [LL] = pop_LL(THETA, D)
%%%
% population POSITIVE log-likelihhod computation with input structure D and
% parameter THETA. The function loops through each track and sums the
% negative LL coming from nLL_kernel function

% Modified for Hessian calculation for log-likelihood

% The input D has fields such as angle dth, input dcp and dc, Basis
% function, regularization lambda, and a mask vector to remove weird or
% missing data
%%%
    LL = 0;
    ntracks = length(D);
    for nn = 1:ntracks
        Z = length(D(nn).dth);  %normalization for nLL per time!
        LL = LL + -nLL_kernel_hist2(THETA, D(nn).dth, D(nn).dcp, D(nn).dc, D(nn).Basis, D(nn).lambda, D(nn).mask);%/Z; % sum across tracks normalize by time
    end
    LL = LL / ntracks;  % normalize by number of tracks

end