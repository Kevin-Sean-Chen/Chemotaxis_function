%Negative Log-likelihood for chemotaxis with kernels
function [NLL] = nLL_Gamma_Mix(THETA, dis, mask)

    %%% loading parameters
    alpha = THETA(1);
    k1 = THETA(2);
    theta1 = THETA(3);
    k2 = THETA(4);
    theta2 = THETA(5);
    
    %%% log-likelihood
    N = length(dis);
    dis = mask(1:end).*dis;  % remove masked elements
    L1 = (k1-1)*nansum(log(dis)) - nansum(dis)/theta1 - N*k1*log(theta1) - N*log(gamma(k1));
    L2 = (k2-1)*nansum(log(dis)) - nansum(dis)/theta2 - N*k2*log(theta2) - N*log(gamma(k2));
    %%% log_likelihood = sum((k-1)*log(data) - data/theta - k*log(theta) - gammaln(k));
    
    marginalP = alpha*L1 + (1-alpha)*L2;
    NLL = -marginalP;
end