function [logli] = logli_trans_opto(mm,xx,yy,mask)
%%%
% Inputs:
% -------
% mm: structure with parameter fields and objective functions
% xx: [3 x T] input with C and dC^p and opto
% yy: [2 x T] ouput with dtheta and dr observations
% mask: [1 x T] logic for likelihood calculations
%
% Output:
%-------
% returns the logli [state x state x T] from solftmax calculation
%%%
    dc = xx(1,:);  % input driving the state transitions
    opto = xx(3,:);  % opto input
    nStates = size(mm.A,2);
    T = length(xx);
    logli = zeros(nStates, nStates, T);
    alpha_ij = mm.wts_state;
    alpha_opto = mm.w_state_opto;
    Basis = mm.basis;
    Kij = zeros(nStates, nStates, size(Basis,1));  % transition filters
    Pij = mm.A;
%     Pij = Pij./sum(Pij, 2);  % transition matrix
    for ii = 1:nStates
        for jj = 1:nStates
            if ii ~=jj
                K_ij_dc = (squeeze(alpha_ij(ii,jj,:))'*Basis');  % reconstruct kernel
                Kij(ii,jj,:) = K_ij_dc;
                K_ij_opto = (squeeze(alpha_opto(ii,jj,:))'*Basis');  % reconstruct kernel
                logli(ii,jj,:) = log(exp(conv_kernel([dc(2:end) dc(1)], K_ij_dc)) + alter_sigmoid(conv_kernel([opto(2:end) dc(1)], K_ij_opto)) + 0*(Pij(ii,jj)));  % test with remove baseline
            else
                logli(ii,jj,:) = 0 + log(Pij(ii,jj));  % remove diagonal components
            end
        end
        logli(ii,:,:) = logli(ii,:,:) - logsumexp(logli(ii,:,:));
    end
end

function nl = alter_sigmoid(x)
    nl = 0.5 + 0.5*x./(1+abs(x));
end