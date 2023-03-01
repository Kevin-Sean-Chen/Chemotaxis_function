%Negative Log-likelihood for angular random walk
% this negLL function is the negative control for continuous random walk without
% chemotaxis components. Should be the lower bound for model comparison
function [NLL] = nLL_randomwalk(THETA, dth, dcp, dc, Basis, lambda, mask)
    
    %%% handling default mask and regularizer
    if nargin < 6
        lambda = 0;
        mask = ones(1,length(dth));
    end
    
    %%% Assume we parameterize in such way first
    mu = THETA(2);  % mean of the angular change
    kappa = THETA(1)^0.5;  % concentration for von Mises angles
    gamma = THETA(13);  % weights on random turns
    
    %%% von Mises probability
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa^2));
    VM = C* exp(kappa^2* cos((dth-mu)*d2r));
    
    %%% marginal probability
    marginalP = (1-gamma)*1/(2*pi) + gamma*VM(2:end);  % marginal LL, indexed to match time step delay
    NLL = -nansum( mask(2:end).* ( log(marginalP + 1*1e-20) ) );
end
