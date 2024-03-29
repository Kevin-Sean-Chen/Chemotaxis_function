function [logpp] = logPrior_mVM(THETA, Basis, lambda)
%%%
% log-prior of parameters for the parameters for mixture of von-Mises
% distribution, used for mcmc sampling for posterior
%%%

%%% information for maximum likelihood
% LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0 -inf, 1e-0, -inf, 1e-1, 1e-0*2, 0.1    -inf -180];
% UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1    inf 180];
% prs0 = [5, 0.01, randn(1,nB)*10, 0.01, 10, 25, 10, 25, 5, 1.    0 10]; 
% K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7)/1; Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); sb = x(12); K2_ = x(13);  gamma = x(14);
%%%

%%% Assume we parameterize in such way first
    kappa_wv = THETA(1)^0.5;      % variance of von Mises
    A_ = THETA(2);            % max turn probability
    alpha_dc = THETA(3:6);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    C_ = THETA(7);            % baseline turning rate
    Amp_dcp = THETA(8);       % amplitude for kernel for dcp normal concentration difference
    tau_dcp = THETA(9);       % time scale for dcp kernel
    Amp_h = THETA(10);        % amplitude of turning history kernel
    tau_h = THETA(11);        % time scale for dth history kernel
    kappa_turn = THETA(12)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(13);        % weight for uniform angle in the turn
    
    base_dc = THETA(14);      % baseline for dc probability
    base_dcp = THETA(15);     % baseline for dcp probability

% bounded parameters
logp_bounds = (kappa_wv>0)&&(kappa_wv<200) && (kappa_turn>0)&&(kappa_turn<200) && (A_-C_>0)&&(A_-C_<=1) && (C_>0)&&(C_<.5) && (tau_h>0)&&(tau_dcp>0) && (gamma>0)&&(gamma<1);
if logp_bounds==0
    logp_bounds = -10^30;
end
% kernel parameters
K_dc = alpha_dc' * Basis';
logp_kernel = -sum((sum(K_dc))^2)*lambda;

logpp = logp_bounds + logp_kernel;
end