%Negative Log-likelihood for angular random walk
% this negLL function is the negative control for continuous random walk without
% chemotaxis components. Should be the lower bound for model comparison
% two parameters: kappa (wv concentration) and gamma (weight on making a turn)
function [NLL] = nLL_randomwalk(THETA, dth, dcp, dc, Basis, lambda, mask)
    
    %%% handling default mask and regularizer
    if nargin < 6
        lambda = 0;
        mask = ones(1,length(dth));
    end
    
    %%% Assume we parameterize in such way first
    Pturn = THETA(2);  % mean of the angular change
    kappa = THETA(1)^0.5;  % concentration for von Mises angles
    kappa_turn = THETA(12)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(13);        % weight for uniform angle in the turn
    
    %%% von Mises probability for run
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa^2));
    VM = C* exp(kappa^2* cos((dth-0)*d2r));
    
    %%% von Mises and uniform for turns
    Ct = 1/(2*pi*besseli(0,kappa_turn^2));
    VM_turn = Ct* exp(kappa_turn^2* cos((dth-pi)*d2r));
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    
    %%% marginal probability
    Pturn = Pturn*14/5;  % temporal step!
    marginalP = Pturn*VM_turn + (1-Pturn)*VM;  % marginal LL, indexed to match time step delay
    NLL = -nansum( mask.* ( log(marginalP + 1*1e-20) ) );
end
