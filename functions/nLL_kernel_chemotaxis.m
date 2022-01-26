%Negative Log-likelihood for chemotaxis with kernels
function [NLL] = nLL_kernel_chemotaxis(THETA, dth, dcp, dc, Basis, lambda)

    if nargin < 6
        lambda = 50;
    end


    %%% Assume we parameterize in such way first
    K_ = THETA(1);  %variance of von Mises
    A_ = THETA(2);  %max turn probability
    B_ = THETA(3:6);  %kernel for dC transitional concentration difference (weights on kerenl basis)
    C_ = THETA(7);  %baseline turning rate
    Amp = THETA(8);  %amplitude for kernel for dcp normal concentration difference
    tau = THETA(9);  %time scale for dcp kernel
    K2_ = THETA(10);  %vairance of the sharp turn von Mises
    
    %[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(5, [0, 10], 1.5);
    B_ = B_ * Basis';  %reconstruct with basis
    B_ = B_ - mean(B_);  %zero-mean to remove dc baseline
    P = A_ ./ (1 + exp( conv( dc, B_, 'same' ))) + C_;  %sigmoid(A_,B_,dc); 
    
    K_win = 1:length(B_);
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,K_^2));  %normalize for von Mises
    D_ = Amp * exp(-K_win/tau);  %dcp kernel
    acdc = conv(dcp, D_, 'same');
    VM = C * exp(K_^2*cos(( dth - acdc )*d2r));  %von Mises distribution
    %%%plot(dth,VM,'o') #check VM shape
%     marginalP = (1-P).*VM + (1/(2*pi))*P;  %marginalized probability
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,K2_^2)) * exp(K2_^2*cos((dth*d2r - pi)));  %test for non-uniform turns (sharp turns)
    gamma = 0.75;
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    
    marginalP = (1-P).*VM + VM_turn.*P;
    
%     lambda = 10;
    NLL = -nansum(log(marginalP + 0*1e-10))  + lambda*(sum((B_ - 0).^2) + C_^2);  % adding slope l2 regularization
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