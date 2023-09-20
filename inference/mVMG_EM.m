%% mVMG_EM
% mixture of von Mesis and Gamma model

%% load some test data
wind = 1:100000;%500000:length(allas);
yy = [allas(wind)*1;
      alldis(wind)];
xx = [alldC(wind); 
      alldcp(wind)];

mask = true(1,length(yy));
mask(isnan(alltrials(wind))) = false;%

%%% for LL0 control
% yy = yy(randperm(length(yy)));

%% Observation and input
% Set parameters: transition matrix and emission matrix
nStates = 4; % number of latent states
nX = 17;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_mVM_gam;  % log-likelihood function

% Set transition matrix by sampling from Dirichlet distr
alpha_diag = 30;  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
% A0 = [0.99,0.01; 0.01,0.99];

% sticky priors
alpha = 1.;  % Dirichlet shape parameter as a prior
kappa = .9;%.5;  % upweighting self-transition for stickiness

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
% wts0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.]; 
wts0 = rand(nY,nX,nStates); % parameters for the mixture-VonMesis behavioral model
wts0(1,:,1) = [20, randn(1,nB)*10, 50, 25, 10, 25, 5, 1.  0.1 0  1 1 0 0]; %[20, randn(1,nB)*10, 10, 25, 10, 25, 1, 1.]; %single mGLM
wts0(1,:,2) = [10,  randn(1,nB)*10, -50, 25, -10, 25, 20,.5  0.1 0  1 1 0 0];
wts0(1,:,3) = [10,  randn(1,nB)*10, -10, 25, -10, 25, 20,.5 0.1 0  1 1 0 0];
wts0(1,:,4) = [10,  randn(1,nB)*10, -20, 25, -20, 25, 20,.5 0.1 0  1 1 0 0];
% wts0(1,:,5) = [3,  randn(1,nB)*10, -20, 25, -20, 25, 20,.5 0.1 0  1 1 0 0];
%[9.9763  -0.5343  -0.0776   0.1238  -0.0529   0.5335   7.7254  367.3817  0.1990  1.0000  0.1000];


% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'loglifun',loglifun,'basis',cosBasis,'lambda',zeros(1,nStates)+0.);

%% Set up variables for EM
maxiter = 50;
EMdisplay = 2;
logpTrace = zeros(maxiter,1); % trace of log-likelihood
dlogp = inf; % change in log-likelihood
logpPrev = -inf; % prev value of log-likelihood
jj = 1; % counter

while (jj <= maxiter) && (dlogp>1e-3)
    
    % --- run E step  -------
    [logp,gams,xisum] = runFB_GLMHMM(mmhat,xx,yy,mask); % run forward-backward
    logpTrace(jj) = logp;
   
    % --- run M step  -------
    
    % Update transition matrix (with stickiness)
    normed_xisum = xisum ./ sum(xisum,2);
    unormedA = kappa*eye(nStates) + (alpha-1)*ones(nStates,nStates) + normed_xisum;
    mmhat.A = unormedA ./ sum(unormedA,2);
%     mmhat.A = (alpha-1 + xisum) ./ (nStates*(alpha-1) + sum(xisum,2)); % normalize each row to sum to 1
    
    %%% sticky one
%     expected_joints = sum([np.sum(Ezzp1, axis=0) for _, Ezzp1, _ in expectations]) + 1e-16
%     expected_joints += self.kappa * np.eye(self.K) + (self.alpha-1) * np.ones((self.K, self.K))
%     P = (expected_joints / expected_joints.sum(axis=1, keepdims=True)) + 1e-16
    
    %%% place-holder for input-driven transitions
    % logp = logp0;
    % logp = logp + w'*u;
    % logp = logp - logsumexp(logp);
    % obj = sum(xis.*logp);l  %(xis is xisum before summation)
    % mmhat = runMstep_trans(mmhat, xx, yy, gams, xis, mask);  % write
    % another function for optimization of driven transition
    
    % Update model params
    mmhat = runMstep_mVM_gam(mmhat, xx, yy, gams, mask);
    %mmhat.Mstepfun(mmhat,xx(:,mask),yy(:,mask),gams(:,mask));
    
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

%% testing
wind_test = max(wind):length(allas);
yy = [allas(wind_test);
      alldis(wind_test)];
xx = [alldC(wind_test); 
      alldcp(wind_test)];
mask = true(1,length(yy));
mask(isnan(alltrials(wind_test))) = false;%
[logp_test,gams_test,xisum_test] = runFB_GLMHMM(mmhat,xx,yy,mask);
logp_test/length(wind_test)

%% sampling tests
n_samp = 10;  % repeats
l_samp = 20000;  % length
test_lls = zeros(1, n_samp);
for nn = 1:n_samp
    randomInteger = randi([max(wind), length(allas)-l_samp]);
    wind_i = randomInteger:randomInteger+l_samp;
    yy = [allas(wind_i)*1;
      alldis(wind_i)];
    xx = [alldC(wind_i); 
          alldcp(wind_i)];
    mask = true(1,length(yy));
    mask(isnan(alltrials(wind_i))) = false;
    [logp_test,gams_test,xisum_test] = runFB_GLMHMM(mmhat,xx,yy,mask);
    test_lls(nn) = logp_test/length(wind_i);
end
test_lls