function [NLL] = pop_nLL(THETA, D)
%%%
% population negative log-likelihhod computation with input structure D and
% parameter THETA. The function loops through each track and sums the
% negative LL coming from nLL_kernel function

% The input D has fields such as angle dth, input dcp and dc, Basis
% function, regularization lambda, and a mask vector to remove weird or
% missing data
%%%
    NLL = 0;
    ntracks = length(D);
    for nn = 1:ntracks
        Z = length(D(nn).dth);  %normalization for nLL per time!
        NLL = NLL + nLL_kernel_hist2(THETA, D(nn).dth, D(nn).dcp, D(nn).dc, D(nn).Basis, D(nn).lambda, D(nn).mask);%/Z; % sum across tracks normalize by time
    end
%     NLL = NLL / ntracks;  % normalize by number of tracks

end