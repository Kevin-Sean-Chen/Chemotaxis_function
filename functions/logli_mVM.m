function [logli] = logli_mVM(mm, xx, yy, mask)

% Compute log-likelihood term under a mixture of von Mesis model
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%      .wts   [1 len(param) K] - per-state parameters for the model
%      .basis [len(alpha) len(kernel)] - basis functions used for the kernel
%      .lambda -scalar for regularization of the logistic fit
%    xx [2 T] - inputs (time series of dc,dcp, and dth)
%    yy [1 T] - outputs (dth angles)
%
% Output
% ------
%  logpy [T K] - loglikelihood for each observation for each state
    
    %%% loading parameters
    THETA = mm.wts;
    Basis = mm.basis;
    lambda = mm.lambda;
    dth = yy;
    dc = xx(1,:);
    dcp = xx(2,:);
    
    K = size(THETA,3);
    lt = length(yy);
    logli = zeros(lt, K);
    for k = 1:K
        %%% Assume we parameterize in such way first
        kappa_wv = THETA(1,1,k)^0.5;      % variance of von Mises
        alpha_dc = THETA(1,2:5,k);    % kernel for dC transitional concentration difference (weights on kerenl basis)
        Amp_dcp = THETA(1,6,k);       % amplitude for kernel for dcp normal concentration difference
        tau_dcp = THETA(1,7,k);       % time scale for dcp kernel
        Amp_h = THETA(1,8,k);         % amplitude of turning history kernel
        tau_h = THETA(1,9,k);         % time scale for dth history kernel
        kappa_turn = THETA(1,10,k)^0.5;   % vairance of the sharp turn von Mises
        gamma = THETA(1,11,k);        % weight for uniform angle in the turn

        %%% turning decision
        %[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(5, [0, 10], 1.5);
        K_dc = alpha_dc * Basis';  %reconstruct with basis
        K_win = 1:length(K_dc);
        h_win = 1:length(K_dc)*1;
        K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel (make longer~)
        filt_dth = conv_kernel(abs(dth), K_h);
        filt_dc = conv_kernel(dc, K_dc);
        P = 1 ./ (1 + exp( -(filt_dc + filt_dth + 0))) +0;  %sigmoid(A_,B_,dc); 

        %%% weathervaning part
        d2r = pi/180;
        C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
        K_dcp = Amp_dcp * exp(-K_win/tau_dcp);    % dcp kernel
        filt_dcp = conv_kernel(dcp, K_dcp);
        VM = C * exp(kappa_wv^2*cos(( dth - filt_dcp )*d2r));  %von Mises distribution

        %%% turning analge model
        VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth*d2r - pi)));  %test for non-uniform turns (sharp turns)
        VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;

        %%% marginal probability
        marginalP = (1-P).*VM + VM_turn.*P;
        logli(:,k) = ( mask.* ( log(marginalP + 1*1e-10) ) ) + lambda(k)/lt*(1*sum((K_dc - 0).^2));% + 0.1*sum((E_ - 0).^2) + 0*C_^2);  % adding slope l2 regularization
    end
end