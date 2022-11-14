%Negative Log-likelihood for chemotaxis with kernels
function [NLL] = nLL_kernel_hist2(THETA, dth, dcp, dc, Basis, lambda, mask)

    if nargin < 6
        lambda = 50;
    end

    %%% Assume we parameterize in such way first
    kappa_wv = THETA(1);      % variance of von Mises
    A_ = THETA(2);            % max turn probability
    alpha_dc = THETA(3:6);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    C_ = THETA(7);            % baseline turning rate
    Amp_dcp = THETA(8);       % amplitude for kernel for dcp normal concentration difference
    tau_dcp = THETA(9);       % time scale for dcp kernel
    Amp_h = THETA(10);        % amplitude of turning history kernel
    tau_h = THETA(11);        % time scale for dth history kernel
%     sb = THETA(12);         % baseline in the nonlinear function
    kappa_turn = THETA(12);   % vairance of the sharp turn von Mises
    gamma = THETA(13);        % weight for uniform angle in the turn
    
    base_dc = 0;%THETA(14);      % baseline for dc probability
    base_dcp = 0;%THETA(15);     % baseline for dcp probability
    
    %%% turning decision
    K_dc = alpha_dc * Basis';  %reconstruct with basis
    K_win = 1:length(K_dc);
    h_win = 1:length(K_dc)*1;
    K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel
%     X_dth = (convmtx(abs(dth),length(E_)));
%     filt_dth = X_dth'*E_';
    filt_dth = conv_kernel(abs([dth(1:end)]), K_h);  %(1:end-1), indexed for history
    filt_dc = conv_kernel([dc(1:end)], K_dc);
    P = A_ ./ (1 + exp( -(filt_dc + filt_dth + base_dc) )) + C_;  %sigmoid(A_,B_,dc); 
%     P = A_ ./ (1 + exp( -(filt_dc + filt_dth + 0))) + C_;
    
    %%% weathervaning part
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    K_dcp = Amp_dcp * exp(-K_win/tau_dcp);    % dcp kernel
    filt_dcp = conv_kernel(dcp, K_dcp);  %(1:end-1)
%     VM = C * exp(kappa_wv^2*cos(( dth(2:end) - filt_dcp - base_dcp )*d2r));  %von Mises distribution
    VM = C * exp(kappa_wv^2*cos(( dth - filt_dcp - base_dcp )*d2r));
    
    %%% turning analge model
%     VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth*d2r - pi)));
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    
    %%% marginal probability
    marginalP = (1-P(1:end)).*VM(1:end) + VM_turn(1:end).*P(1:end);  % marginal LL, indexed to match time step delay
%     lambda = 10;
    NLL = -nansum( mask(1:end).* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));% + 0.1*sum((E_ - 0).^2) + 0*C_^2);  % adding slope l2 regularization
%     NLL = -nansum( mask.* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));%
end

% def nLL(THETA, dth,dcp,dc):
%     """
%     negative log-likelihood objective function for fitting
%     THETA includes parameter to be inferred and dth, dcp, dc are from recorded data
%     """
%     #a_, k_, A_, B_, C_, D_ = THETA  #inferred paramter
%     k_, A_, B_, Amp, tau = THETA[0], THETA[1], THETA[2:8], THETA[8], THETA[9] #, THETA[8]#, THETA[9] #Kappa,A,Kdc,Kdcp,dc_amp,dcp_amp
%     #B_ = np.dot(B_,RaisedCosine_basis(len(K_win),5))  #test with basis function
%     #B_ = 100* B_/np.linalg.norm(B_)
%     #P = sigmoid(A_, B_, C_, D_, dcp)
%     P = sigmoid2(A_,B_,dc)
%     #VM = np.exp(k_*np.cos((dth-a_*dcp)*d2r)) / (2*np.pi*iv(0,k_))#von Mises distribution
%     #vm_par = vonmises.fit((dth-a_*dcp)*d2r, scale=1)
%     rv = vonmises(k_)#(vm_par[0])
%     C_ = -Amp *np.exp(-K_win/tau)  #sign change due to the way simulated above
%     VM = rv.pdf((dth-np.dot(dcp,C_))*d2r)
%     marginalP = np.multiply((1-P), VM) + (1/(2*np.pi))*P
%     nll = -np.sum(np.log(marginalP+1e-9))#, axis=1)
%     #fst = np.einsum('ij,ij->i', 1-P, VM)
%     #snd = np.sum(1/np.pi*P, axis=1)
%     return np.sum(nll)