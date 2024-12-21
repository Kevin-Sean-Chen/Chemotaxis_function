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
%     B_ = B_ - mean(B_);  %zero-mean to remove dc baseline
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
    gamma = .5;
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    
    marginalP = (1-P).*VM + VM_turn.*P;
    
%     lambda = 10;
    NLL = -nansum(log(marginalP + 0*1e-10))  + lambda*(sum((B_ - 0).^2) + C_^2);  % adding slope l2 regularization
end
