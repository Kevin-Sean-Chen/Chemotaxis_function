function [logli] = logli_staPAW(mm, xx, yy, mask)

% Compute log-likelihood term under a mixture of von Mesis angle + gamma
% speed model; both under a shared sigmoid
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%      .wts   [1 len(param) K] - per-state parameters for the model
%      .basis [len(alpha) len(kernel)] - basis functions used for the kernel
%      .lambda -scalar for regularization of the logistic fit
%    xx [2 T] - inputs (time series of dc,dcp, and dth)
%    yy [2 T] - outputs (dth angles and dis speed displacement)
%
% Output
% ------
%  logpy [T K] - loglikelihood for each observation for each state
    
    %%% loading parameters
    THETA = mm.wts;
    Basis = mm.basis;
    lambda = mm.lambda;
    dth = yy(1,:);
    dv = yy(2,:);
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
        A_ = THETA(1,12,k);           % max turning probability
        B_ = THETA(1,13,k);           % baseline turning probability
        gamm_shapes = THETA(1,14:15,k);       % shape parameters for speed distribution
        gamm_scales = THETA(1,16:17,k);       % scale parameters for speed distribution
%         gamm_AR = THETA(1,18:19,k);           % AR weight for speed distribution
%         base_dc = THETA(1,20,k);      % bias of turning
%         base_dcp = THETA(1,21,k);     % bias of curving
        gamm_AR = [0,0];
        base_dc = THETA(1,18,k);      % bias of turning
        base_dcp = THETA(1,19,k);     % bias of curving
        

        %%% turning decision
        %[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(5, [0, 10], 1.5);
        K_dc = alpha_dc * Basis';  %reconstruct with basis
        K_win = 1:length(K_dc);
        h_win = 1:length(K_dc)*1;
        K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel (make longer~)
        filt_dth = conv_kernel(abs(dth(1:end-1)), K_h);
        filt_dc = conv_kernel(dc(2:end), K_dc);
        P = (A_-B_) ./ (1 + exp( -(filt_dc + filt_dth + base_dc))) +B_;  %sigmoid(A_,B_,dc); 
%         P = P*14/5;
        
        %%% weathervaning part
        d2r = pi/180;
        C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
        K_dcp = Amp_dcp * exp(-K_win/tau_dcp);    % dcp kernel
        filt_dcp = conv_kernel(dcp(2:end), K_dcp);
        VM = C * exp(kappa_wv^2*cos(wrapTo180( dth(2:end) - filt_dcp - base_dcp)*d2r));  %von Mises distribution
%         VM = C * exp(kappa_wv^2*cos(( dth(2:end) - filt_dcp - base_dcp)*d2r));
        gamm_wv = gampdf(dv(2:end), gamm_shapes(2) + gamm_AR(2)*dv(1:end-1), gamm_scales(2));
        
        %%% turning analge model
        VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
        VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
        gamm_pr = gampdf(dv(2:end), gamm_shapes(1) + gamm_AR(1)*dv(1:end-1), gamm_scales(1));

        %%% marginal probability
        marginalP = (1-P).*VM.*gamm_wv + P.*VM_turn.*gamm_pr;
%         pos = (~isinf(log(dv)));
%         ll_gamm = (gamm_shape-1)*log(dv.*pos) - dv/gamm_scale - gamm_shape*log(gamm_scale) - gammaln(gamm_shape);
        if K>1
            logli(2:end,k) = ( mask(2:end).* ( log(marginalP + 1*1e-10) ) ) + lambda(k)*(1*sum((K_dc - 0).^2))/lt;
        elseif K==1
            logli(2:end) = (( mask(2:end).* ( log(marginalP + 1*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2))/lt)';
        end
        
    end
end