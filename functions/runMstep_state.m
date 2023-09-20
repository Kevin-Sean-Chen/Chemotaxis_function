function mm = runMstep_state(mm, xx, yy, xis, mask)
%%%
% M-step to optimize input-driven state transitions in HMM
% take the structure mm, input xx, output yy, expected joints E(z',z) xis,
% and mask logic as input.
% Return updated mm structure with new parameters
%%%
    
    dc = xx(1,:);  % input driving the state transitions
    nStates = size(mm.wts_state,1);
    Basis = mm.basis;
    Pij_ = mm.A;  % baseline transitions
    lfun = @(x)nll_MLR_(x, dc, xis, Basis, mask);  % multi-nomial logistic regression loss function
    prs0 = [reshape(mm.wts_state,1,[])    reshape(Pij_',1,[])];  % concatenate flatten parameters
    LB = [zeros(1,length(prs0))-100    zeros(1,nStates^2)];
    UB = [zeros(1,length(prs0))+100    ones(1,nStates^2)];
    Aeq = zeros(nStates, length(prs0));   % transition probability sums to one!
    lw = length(reshape(mm.wts_state,1,[]));
    for ii = 1:nStates
        Aeq(ii,lw+(ii-1)*nStates+1:lw+ii*nStates) = ones(1,nStates);  % "smart" parameterizing with the last nStates^2 being state transitions
    end
    beq = ones(nStates,1);  % this all makes the constrain: Ax = b
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],Aeq,beq,LB,UB,[]);
    mm.wts_state = reshape(x(1:end-nStates^2), nStates, nStates, []);   % update weights
    mm.A = reshape(x(end-nStates^2+1:end), nStates, nStates)';  % update baseline, note the transpose here, at flatten, and corresponding to the constaints, important!!
    
end

function nll = nll_MLR_(x, dc, xis, Basis, mask)
%%%
% Multinomial logistic regession, with parameter vector x, input dc,
% expected joint xis, Basis vectors, and mask logic as input
% return the negative log-liklihood for optimization
% -sum_t(sum_state(sum_state(xi(s,s,t)*log(alpha(s,s,t))))), where alpha is
% the soft-max transition probability
%%% 
    logP_ = xis*0;  % reconstructed transitional probability (with driven and baseline)
    nStates = size(xis,1);
    Kij = zeros(nStates, nStates, size(Basis,1));
    alpha_ij = reshape(x(1:end-nStates^2), nStates, nStates,[]);  % weights on basis for kernels
    Pij = reshape(x(end-nStates^2+1:end), nStates, nStates)';  % basline probability, important to transpose given the way its flattened!!
%     Pij = Pij + kappa*eye(nStates);  % stickiness
    Pij = Pij./sum(Pij, 2);  % transition matrix
    for ii = 1:nStates
        for jj = 1:nStates
            if ii ~= jj
                K_ij_dc = (squeeze(alpha_ij(ii,jj,:))'*Basis');  % reconstruct kernel
                Kij(ii,jj,:) = K_ij_dc;
                logP_(ii,jj,:) = conv_kernel(dc, K_ij_dc) + log(Pij(ii,jj));
            else
                logP_(ii,jj,:) = 0 + log(Pij(ii,jj));
            end
        end
        logP_(ii,:,:) = logP_(ii,:,:) - logsumexp(logP_(ii,:,:));
    end
    
    nll = -squeeze(sum(sum(xis(:,:,1:end) .* logP_(:,:,1:end))))'*mask(1:end)';
    
end