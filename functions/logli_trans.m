function [logli] = logli_trans(mm,xx,yy,mask)
%%%
% Inputs:
% -------
% mm: structure with parameter fields and objective functions
% xx: [2 x T] input with C and dC^p
% yy: [2 x T] ouput with dtheta and dr observations
% mask: [1 x T] logic for likelihood calculations
%
% Output:
%-------
% returns the logli [state x state x T] from solftmax calculation
%%%
    dc = xx(1,:);  % input driving the state transitions
    nStates = size(mm.A,2);
    T = length(xx);
    logli = zeros(nStates, nStates, T);
    alpha_ij = mm.wts_state;
    Basis = mm.basis;
    Kij = zeros(nStates, nStates, size(Basis,1));  % transition filters
    Pij = mm.A;
    Pij = Pij./sum(Pij, 2);  % transition matrix
    for ii = 1:nStates
        for jj = 1:nStates
            if ii ~=jj
                K_ij_dc = (squeeze(alpha_ij(ii,jj,:))'*Basis');  % reconstruct kernel
                Kij(ii,jj,:) = K_ij_dc;
                logli(ii,jj,:) = conv_kernel(dc, K_ij_dc) + log(Pij(ii,jj));
            else
                logli(ii,jj,:) = 0 + log(Pij(ii,jj));  % remove diagonal components
            end
        end
        logli(ii,:,:) = logli(ii,:,:) - logsumexp(logli(ii,:,:));
    end
end
