function [logpp] = logPrior_opto(THETA, Basis, lambda)
%%%
% log-prior of parameters for the parameters for concentration and opto
% input, used for mcmc sampling for posterior of these kernels
%%%

%%% information for maximum likelihood
% LB = [1e-0, 1e-1, ones(1,nB*2)*-inf, 0, 1., 0.1];
% UB = [200, 1., ones(1,nB*2)*inf, 0.1, 20, 1.];
% prs0 = [50, 0.1, randn(1,nB*2)*1, 0.01, 5, .5];
% K_ = x(1); A_ = x(2); K_c = x(3:6); K_opto = x(7:10); C_ = x(11); K2_ = x(12); gamma = x(13);
%%%

%%% Assume we parameterize in such way first
    kappa_wv = THETA(1)^0.5;      % variance of von Mises
    A_ = THETA(2);            % max turn probability
    alpha_dc = THETA(3:6);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    alpha_opto = THETA(7:10);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    C_ = THETA(11);            % baseline turning rate
    kappa_turn = THETA(12)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(13);        % weight for uniform angle in the turn

% bounded parameters
logp_bounds = (kappa_wv>0)&&(kappa_wv<200) && (kappa_turn>0)&&(kappa_turn<200) && (A_-C_>0)&&(A_-C_<=1) && (C_>0)&&(C_<.5) && (gamma>0)&&(gamma<1);
if logp_bounds==0
    logp_bounds = -10^30;
end
% kernel parameters
K_dc = alpha_dc' * Basis';
logp_kernel_c = -sum((sum(K_dc))^2)*lambda;
K_opto = alpha_opto' * Basis';
logp_kernel_opto = -sum((sum(K_dc))^2)*lambda;

logpp = logp_bounds + logp_kernel_c + logp_kernel_opto;
end