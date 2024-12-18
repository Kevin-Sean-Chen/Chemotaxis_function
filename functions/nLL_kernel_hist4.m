% Negative Log-likelihood for chemotaxis with kernels, all with basis functions
function [NLL] = nLL_kernel_hist4(THETA, dth, dcp, dc, Basis, lambda, mask)
    
    %%% regularization
    if nargin < 6
        lambda = 0;
    end

    %%% Assume we parameterize in such way first
    alpha_h = THETA(1:4);
    alpha_dc = THETA(5:8);
    alpha_dcp = THETA(9:12);
    kappa_turn = THETA(13)^0.5;
    kappa_wv = THETA(14)^0.5;
    gamma = THETA(15);
    A = THETA(16);
    B = THETA(17);
    
    %%% kernel with basis
    K_h = (alpha_h*Basis');  % dth kernel
    K_dc = (alpha_dc*Basis');  % dC kernel
    dcp_win = 1:length(K_dc)-1;
    K_dcp = THETA(11) * exp(-dcp_win/THETA(12)); 
    
    %%% turning decision
    d2r = pi/180;  
    filt_dth = conv_kernel(abs(dth(1:end-1))*1,K_h);
    filt_dc = conv_kernel(dc(2:end),K_dc);
    P = (A-B)./(1+exp(-(filt_dth + filt_dc))) + B;
    
    %%% weathervaning part
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    filt_dcp = conv_kernel(dcp(2:end),K_dcp);
    VM = C * exp(kappa_wv^2*cos(( dth(2:end) - filt_dcp )*d2r));  %von Mises distribution
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
    VM_turn = gamma*1/(2*pi) + (1-gamma)*VM_turn;  %%% revisit mixture inference !!!
    
    marginalP = (1-P).*VM + VM_turn.*P;
    
    NLL = -nansum(mask(2:end).*log(marginalP + 1*1e-20)) + lambda*sum(K_dc.^2);  % adding slope l2 regularization
end
