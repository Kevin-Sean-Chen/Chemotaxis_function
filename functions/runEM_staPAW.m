function [mmhat, logpTrace] = runEM_staPAW(Data_train, K_states)

%% runEM_staPAW
% input training Data structure and number of states pre-defined
% run EM for the mixture of von Mesis and Gamma model, with input-driven states!
% dr and dtheta are both under a sigmoid decision!
% the output is the optimized parameter and logP trace
% this function is called by scan_staPAW to cross-validate and scan through states

%% load training data
[xxf, yyf, alltrials, time] = data2xy(Data_train);  %Data
alldis = extractfield(Data_train, 'dis');  % Data

% creating new xx,yy for fitting if needed
yyf = [yyf; alldis];
% xxf = [xxf; allopto];
maskf = alltrials;  %((alltrials)==0) = false;%

%%% not picking a window
xx = xxf;
yy = yyf;
mask = maskf;

%% Observation and input
% Set parameters: transition matrix and emission matrix
nStates = K_states; % number of latent states
nX = 17 + 2 + 0;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_staPAW;  % log-likelihood function
loglitrans = @logli_trans;  % log-likelihood of soft-max state transitions

% Set transition matrix by sampling from Dirichlet distr
alpha_diag = 60;%60  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;%5  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1

% sticky priors
alpha = 1.;  % Dirichlet shape parameter as a prior
kappa = .5;  % .5 upweighting self-transition for stickiness

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
wts0 = rand(nY, nX, 5); % parameters for the mixture-VonMesis behavioral model, assuming we never check above 5 states...
%%% kappa_wv, alpha_Kc, a_dc, tau_dc, a_h, tau_h, gamma, kappa_brw, A, B, k1, k2, theta1, theta2, phi1, phi2, base_c, base_dcp
wts0(1,:,1) = [50,  randn(1,nB)*10,  10, 25,  10, 25,  5,   1.  0.1 0  1 1 1 1  0 0];
wts0(1,:,2) = [10,  randn(1,nB)*10, -10, 25, -10, 25, 20,  .5   0.1 0  1 1 1 1  0 0];
wts0(1,:,3) = [10,  randn(1,nB)*10, -10, 25, -10, 25, 20,  .5   0.1 0  1 1 1 1  0 0];
wts0(1,:,4) = [10,  randn(1,nB)*10, -20, 25, -20, 25, 20,  .5   0.1 0  1 1 1 1  0 0];
wts0(1,:,5) = [10,  randn(1,nB)*10, -20, 25, -20, 25, 20,  .5   0.1 0  1 1 1 1  0 0];
wts0 = wts0(:,:,1:nStates);

%%% transition kernels
alpha_ij = randn(nStates,nStates, nB)*.1;  % state x state x weights
lk = size(cosBasis,1);
K_ij = zeros(nStates, nStates, lk);
for ii = 1:nStates
    for jj = 1:nStates
        if ii~=jj
            K_ij(ii,jj,:) = fliplr(squeeze(alpha_ij(ii,jj,:))'*cosBasis');  % build state-transition kernels, without self-kernels
        else
            alpha_ij(ii,jj) = alpha_ij(ii,jj) *0;
        end
    end
end
wts_state0 = alpha_ij;

% Build struct for initial params
A0 = A0*50;  % test with non-normalized initial weights
mmhat = struct('A',A0,'wts',wts0,'wts_state',wts_state0,'loglifun',loglifun,'loglitrans',loglitrans,'basis',cosBasis,'lambda',zeros(1,nStates)+0.0);

%% Set up variables for EM
maxiter = 50;
EMdisplay = 2;
logpTrace = zeros(maxiter,1); % trace of log-likelihood
dlogp = inf; % change in log-likelihood
logpPrev = -inf; % prev value of log-likelihood
jj = 1; % counter

while (jj <= maxiter) && (dlogp>1e-3)
    
    % --- run E step  -------
    [logp,gams,xis,xisum] = runFB_GLMHMM_xi(mmhat,xx,yy,mask); % run forward-backward
    logpTrace(jj) = logp;
   
    % --- run M step  -------
    
    % Update transition matrix (with stickiness)
%     normed_xisum = xisum ./ sum(xisum,2);
%     unormedA = kappa*eye(nStates) + (alpha-1)*ones(nStates,nStates) + normed_xisum;
%     mmhat.A = unormedA ./ sum(unormedA,2);
%     mmhat.A = (alpha-1 + xisum) ./ (nStates*(alpha-1) + sum(xisum,2)); % normalize each row to sum to 1
    
    % Update model params
    mmhat = runMstep_staPAW(mmhat, xx, yy, gams, mask);
    % Update for input-driven transitions
    if K_states > 1   % only run state-update if there is HMM (K>1)
        mmhat = runMstep_state(mmhat, xx, yy, xis, mask); 
    elseif K_states==1
        break
    end
    
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
% wind_test = max(wind):max(wind)+100000;%length(alldis)-30000:length(alldis)-10000; %length(alldis);%1:max(wind); %
% yy = yyf(:,wind_test);
% xx = xxf(:,wind_test);
% mask = maskf(wind_test);
% [logp_test,gams_test,xis_test,xisum_test] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
% logp_test/length(wind_test)

end