function [logp,gams,xisum,logcs] = runFB_GLMHMM(mm,xx,yy,mask)
% [logp,gams,xisum,logcs] = runFB_GLMHMM(mm,xx,yy,mask)
% 
% Run forward-backward algorithm (E-Step) for GLM-HMM with discrete
% observations using logSumExp trick (to avoid underflow)
%
% - uses formulas for hat-alpha and hat-Beta from Bishop (pg 627-629)
% 
% - Uses logSumExp trick (which is slightly slower) to avoid underflow
% errors
%
%
% Inputs:
% -------
%   mm - parameter struct
%      .A    [k x k] - transition matrix (k states)
%      .mus  [1 x k] - mean for each state
%      .vars [1 x k] - variance for each state
%   rr [1 x T] - residual for each time bin
%
% Outputs:
% -------
%    logp [1 x 1] - log marginal probability: log p(X|theta)
%    gams [k x T] - marginal posterior over state at each time: p(z_t|X)
%  xissum [k x k] - summed marginal over pairs of states: p(z_t, z_{t+1}|X)
%   logCs [1 x T] - log conditional marginal li: log p(y_t | y_{1:t-1}) 
%
% Note: requires MATLB v2016b or later (for expansion operations)


% Extract sizes
nStates = size(mm.A,1);  
nT = size(yy,2);

if nargin < 4
    mask = true(1,nT);
end

% Set initial state probability
logpi0 = log(ones(1,nStates)/nStates); % assume uniform prior over initial state

% Compute log-likelihood for each observation under each state
logpy = mm.loglifun(mm,xx,yy,mask);

%% Forward pass
logaa = zeros(nStates,nT); % forward probabilities p(z_t | y_{1:t})
logbb = zeros(nStates,nT); % backward: p(y_{t+1:T} | z_t)/p(y_{t+1:T} | y_{1:T})
logcs = zeros(1,nT); % forward marginal likelihoods: p(y_t|y_{1:t-1})

% First bin
if mask(1)    
    logpyz = logpi0 + logpy(1,:); % log joint: log P(y_1,z_1)
    logcs(1) = logsumexp(logpyz,2); % log-normalizer
    logaa(:,1) = logpyz' - logcs(1); % log conditional log P(z_1|y_1)
else
    logcs(1) = 0;
    logaa(:,1) = logpi0;
end

% Remaining time bins
for jj = 2:nT
    logaaprior = log(exp(logaa(:,jj-1))'*mm.A);  % propagate uncertainty forward
    if mask(jj) % Include likelihood term
        logpyz = logaaprior + logpy(jj,:); % joint log P(y_{1:t},z_t)
        logcs(jj) = logsumexp(logpyz,2);  % conditional log P(y_t|y_{1:t-1})
        logaa(:,jj) = logpyz' - logcs(jj); % conditional log P(z_t|y_{1:t})
    else  % ignore this likelihood term
        logcs(jj)=0;
        logaa(:,jj) = logaaprior;
    end
end

%% Backward pass 
for jj = nT-1:-1:1
    if mask(jj+1) % Include likelihood term
        logbbpy = logbb(:,jj+1) + logpy(jj+1,:)';
        logbb(:,jj) = log(mm.A*exp(logbbpy)) - logcs(jj+1);
    else  % ignore this likelihood term
        logbb(:,jj) = log(mm.A*exp(logbb(:,jj+1))) - logcs(jj+1);
    end
end

%% Compute outputs

% marginal likelihood
logp = sum(logcs);

% marginal posterior over states p(z_t | Data)
loggams = logaa + logbb; 
gams = exp(loggams);

% Compute only sum of xis 
logA = log(mm.A);
if nargout > 2
     xisum = zeros(nStates,nStates); % P(z_n-1,z_n|y_1:T)
     for jj = 1:nT-1
         
         if mask(jj+1)
             % % version without logSumExp trick:
             % xisum = xisum + (aa(:,jj) * (bb(:,jj+1)'.*Py(jj+1,:,ydat(:,jj+1))) .*A)/cs(jj+1);

             % Using logSumExp:
             xisum = xisum + exp(logA + (logaa(:,jj)-logcs(jj+1)) + (logbb(:,jj+1)' + logpy(jj+1,:)));
         else
             % % version without logSumExp trick:
             % xisum = xisum + ((aa(:,jj)*bb(:,jj+1)').*A)/cs(jj+1);

             % Using logSumExp:
             xisum = xisum + exp(logA + (logaa(:,jj)-logcs(jj+1)) + (logbb(:,jj+1)'));
         end
     end

end
