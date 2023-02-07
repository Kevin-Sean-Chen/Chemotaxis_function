function [mmhat, logpTrace] = runEM_mVMHMM(Data, xx, yy, mask, options)
%%%
% function for running EM algorithm to infer mVM-HMM model
% Input: Data structure with track data
%        xx input vector [2 x N*t]
%        yy output vector [1 x N*T]
%        mask logic vector [1 x N*T]
%        options structure with inference parameters
%%%

% unpack options
nStates = options.nStates;
alpha = options.alpha; % Dirichlet shape parameter as a prior
kappa = options.kappa; % upweighting self-transition for stickiness
lambs = options.lambs;
maxiter = options.maxiter;
dlogp_tol = options.dlogp_tol;

% default function and parameters
nX = 11;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
loglifun = @logli_mVM;  % log-likelihood function
alpha_diag = 30;  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
% A0 = [0.99,0.01; 0.01,0.99];
nB = 4;  % number of basis used for kernel
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);% basis function

% Set linear weights & output noise variances
wts0 = rand(nY, nX, 5); % parameters for the mixture-VonMesis behavioral model (put in 5 as we likely won't use up to 5 states anyway)
wts0(1,:,1) = [20, randn(1,nB)*10, 50, 25, 10, 25, 5, 1.]; %[20, randn(1,nB)*10, 10, 25, 10, 25, 1, 1.]; %single mGLM
wts0(1,:,2) = [10,  randn(1,nB)*10, -50, 25, -10, 25, 20,.5];
wts0(1,:,3) = [1,  randn(1,nB)*10, -10, 25, -10, 25, 20,.5];
wts0 = wts0(:,:,1:nStates);

% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'loglifun',loglifun,'basis',cosBasis,'lambda',lambs);

% The main EM loop
EMdisplay = 3;
logpTrace = zeros(maxiter,1); % trace of log-likelihood
dlogp = inf; % change in log-likelihood
logpPrev = -inf; % prev value of log-likelihood
jj = 1; % counter

while (jj <= maxiter) && (dlogp > dlogp_tol)
    
    % --- run E step  -------
    [logp, gams, xisum] = runFB_GLMHMM(mmhat,xx,yy,mask); % run forward-backward
    logpTrace(jj) = logp;
   
    % --- run M step  -------
    
    % Update transition matrix (with stickiness)
    normed_xisum = xisum ./ sum(xisum,2);
    unormedA = kappa*eye(nStates) + (alpha-1)*ones(nStates,nStates) + normed_xisum;
    mmhat.A = unormedA ./ sum(unormedA,2);
    
    % Update model params
    mmhat = runMstep_mVM(mmhat, xx, yy, gams, mask);
    
    % ---  Display progress ----
    if mod(jj,EMdisplay)==0
        fprintf('EM iter %d:  logli = %-.6g\n',jj,logp);
    end
    
    % Update counter and log-likelihood change
    jj = jj+1;  
    dlogp = logp-logpPrev; % change in log-likelihood
    logpPrev = logp; % previous log-likelihood

    if dlogp<-1e-6
        warning('Log-likelihood decreased during EM!');
        fprintf('dlogp = %.5g\n', dlogp);
    end

end
jj = jj-1;


end