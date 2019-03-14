%Negative Log-likelihood for chemotaxis
function [NLL] = nLL_chemotaxis(THETA, dth, dcp, dc)

end

% def nLL(THETA, dth,dcp,dc):
%     #a_, k_, A_, B_, C_, D_ = THETA  #inferred paramter
%     a_,k_,A_,B_ = THETA
%     #P = sigmoid(A_, B_, C_, D_, dcp)
%     P = sigmoid2(A_,B_,dc)
%     #VM = np.exp(k_*np.cos((dth-a_*dcp)*d2r)) / (2*np.pi*iv(0,k_))#von Mises distribution
%     #vm_par = vonmises.fit((dth-a_*dcp)*d2r, scale=1)
%     rv = vonmises(k_)#(vm_par[0])
%     VM = rv.pdf((dth-a_*dcp)*d2r)
%     marginalP = np.multiply((1-P), VM) + (1/(2*np.pi))*P
%     nll = -np.sum(np.log(marginalP+1e-7))#, axis=1)
%     #fst = np.einsum('ij,ij->i', 1-P, VM)
%     #snd = np.sum(1/np.pi*P, axis=1)
% return np.sum(nll)