%Negative Log-likelihood for chemotaxis with kernels
function [NLL] = nLL_kernel_chemotaxis(THETA, dth, dcp, dc, Basis)

    %%% Assume we parameterize in such way first
    K_ = THETA(1);
    A_ = THETA(2);
    B_ = THETA(3:8);
    Amp = THETA(9);
    tau = THETA(10);
    
    %[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(5, [0, 10], 1.5);
    B_ = B_ * Basis';
    P = A_ ./ (1 + exp( conv( dc, B_, 'same' )));  %sigmoid(A_,B_,dc); 
    
    K_win = 1:length(B_);
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,K_));
    C_ = Amp * exp(-K_win/tau);
    acdc = conv(dcp, C_, 'same');
    VM = C * exp(K_*cos(( dth - acdc )*d2r));  %von Mises distribution
    %%%plot(dth,VM,'o') #check VM shape
    marginalP = (1-P).*VM + (1/(2*pi))*P;  %marginalized probability
    NLL = -nansum(log(marginalP + 1e-10));
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