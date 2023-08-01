%Negative Log-likelihood for binomial model for chemotaxs
% this negLL function is the benchmark comparison to our statistical model
% it takes the 'Data' stucture in, counts how many tracks with .dc goes up,
% then estimates and returns the probability (chemotaxis index) p_hat

% CI = up-down/(up+down) = 2up/(up+down) = 2p_hat - 1; -> same as
% esitimating p_hat: probability of going up gradient

function [NLL] = pop_nLL_binomial(THETA, Data)
    
    p_hat = THETA;
    n = length(Data);
    y = 0;
    for nn = 1:n
        dc_i = Data(nn).dc;
        if dc_i(end) > dc_i(1)
            y = y+1;   % counts up gradient
        end
    end
    NLL = y*log(p_hat) + (n-y)*log(1-p_hat);  % nll of Binomial distribution
    NLL = -NLL;
    
end