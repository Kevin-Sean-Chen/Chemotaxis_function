%Negative Log-likelihood for chemotaxis with kernels
%%% negative control without using Kc kernel
function [NLL] = nLL_wo_kc(THETA, dth, dcp, dc, Basis, lambda, mask)

    if nargin < 6
        lambda = 0;
        mask = ones(1,length(dth));
    end
% mask = ones(1,length(dth));
    %%% Assume we parameterize in such way first
    kappa_wv = THETA(1)^0.5;      % variance of von Mises
    A_ = THETA(2);            % max turn probability
    alpha_dc = THETA(3:6);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    C_ = THETA(7);            % baseline turning rate
    Amp_dcp = THETA(8);       % amplitude for kernel for dcp normal concentration difference
    tau_dcp = THETA(9);       % time scale for dcp kernel
    Amp_h = THETA(10);        % amplitude of turning history kernel
    tau_h = THETA(11);        % time scale for dth history kernel
%     sb = THETA(12);         % baseline in the nonlinear function
    kappa_turn = THETA(12)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(13);        % weight for uniform angle in the turn
    
    base_dc = 0;%THETA(3);      % baseline for dc probability
    base_dcp = 0;%THETA(4);     % baseline for dcp probability
    
    %%% turning decision
    try
        K_dc = alpha_dc * Basis';  %reconstruct with basis
    catch
        K_dc = alpha_dc' * Basis';
    end
    lk = size(Basis,1); 
    K_win = 1:length(lk)-1;
    h_win = 1:length(lk)-1;
    K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel
%     X_dth = (convmtx(abs(dth),length(E_)));
%     filt_dth = X_dth'*E_';
    filt_dth = conv_kernel([abs(dth(1:end-1))], K_h);  %(1:end-1), indexed for history
    filt_dc = conv_kernel([dc(2:end)], K_dc);
    P = (A_-C_)*1 ./ (1 + exp( -(filt_dc(randperm(length(filt_dc)))*0 + filt_dth + base_dc) )) + C_;  %sigmoid(A_,B_,dc); 
    P = P*14/5;  % temporal step!
%     P = A_ ./ (1 + exp( -(filt_dc + filt_dth + 0))) + C_;
    
    %%% weathervaning part
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    K_dcp = Amp_dcp * exp(-K_win/tau_dcp);    % dcp kernel
    filt_dcp = conv_kernel(dcp(2:end), K_dcp);  %
    VM = C * exp(kappa_wv^2*cos(( dth(2:end) - filt_dcp - base_dcp )*d2r));  %von Mises distribution
%     VM = C * exp(kappa_wv^2*cos(( dth - filt_dcp - base_dcp )*d2r));
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
%     VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth*d2r - pi)));
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    
    %%% marginal probability
    marginalP = (1-P(1:end)).*VM(1:end) + VM_turn(1:end).*P(1:end);  % marginal LL, indexed to match time step delay
%     lambda = 10;
    NLL = -nansum( mask(2:end).* ( log(marginalP + 1*1e-20) ) ) + 1*lambda*(1*sum((K_dc - 0).^2));% + 0.1*sum((E_ - 0).^2) + 0*C_^2);  % adding slope l2 regularization
%     NLL = -nansum( mask.* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));%
end