% Negative Log-likelihood for chemotaxis with kernel for optogenetic
% stimuli, which only acts on the biased-random walk strategy
function [NLL] = nLL_kernel_opto(THETA, dth, dc, opto, Basis, lambda, mask)
    
    if nargin < 6
        lambda = 0;
        mask = ones(1,length(dth));
    end

    %%% loading parameters
    kappa_wv = THETA(1)^0.5;  % variance of von Mises
    A_ = THETA(2);            % max turn probability
    alpha_dc = THETA(3:6);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    alpha_opto = THETA(7:10); % kernel for opto input
    C_ = THETA(11);            % baseline turning rate
%     Amp_h = THETA(8);        % amplitude of turning history kernel
%     tau_h = THETA(9);        % time scale for dth history kernel
    kappa_turn = THETA(12)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(13);        % weight for uniform angle in the turn
    
    %%% turning decision
    try
        K_dc = alpha_dc * Basis';  %reconstruct with basis
        K_opto = alpha_opto * Basis'; 
    catch
        K_dc = alpha_dc' * Basis';
        K_opto = alpha_opto' * Basis'; 
    end
%     h_win = 1:length(K_dc)-1;
%     K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel
%     filt_dth = conv_kernel([abs(dth(1:end-1))], K_h);  %(1:end-1), indexed for history
    filt_dc = conv_kernel([dc(2:end)], K_dc);
    filt_opto = conv_kernel([opto(2:end)], K_opto);
    P = (A_-C_)*1 ./ (1 + exp( -(filt_dc + filt_opto) )) + C_;  %sigmoid(A_,B_,dc); 
    P = P*14/5;  % temporal step!
    
    %%% weathervaning part
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    VM = C * exp(kappa_wv^2*cos(( dth(2:end) - 0 )*d2r));  %von Mises distribution
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
    VM_turn = (1-gamma)*1/(2*pi) + (gamma)*VM_turn;  %%% revisit mixture inference !!!
    
    marginalP = (1-P).*VM + VM_turn.*P;
    
    NLL = -nansum(mask(2:end).*log(marginalP + 1*1e-20)) + lambda*sum(K_dc.^2 + K_opto.^2);  % adding slope l2 regularization
end
