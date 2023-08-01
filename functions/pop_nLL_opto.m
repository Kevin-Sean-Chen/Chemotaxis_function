function [NLL] = pop_nLL_opto(THETA, D)
%%%
% population negative log-likelihhod computation with input structure D and
% parameter THETA. The function loops through each track and sums the
% negative LL coming from nLL_kernel function, this is for the opto data

% The input D has fields such as angle dth, input opto, Basis
% function, regularization lambda, and a mask vector to remove weird or
% missing data
%%%
    NLL = 0;
    ntracks = length(D);
    for nn = 1:ntracks
        Z = length(D(nn).dth);  %normalization for nLL per time!
        NLL = NLL + nLL_kernel_opto(THETA, D(nn).dth, D(nn).dc, D(nn).opto, D(nn).Basis, D(nn).lambda, D(nn).mask);%/Z; % sum across tracks normalize by time
    end
%     NLL = NLL / ntracks;  % normalize by number of tracks

end