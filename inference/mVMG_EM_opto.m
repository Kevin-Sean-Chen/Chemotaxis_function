%% mVMG_EM_opto
% mixture of von Mesis and Gamma model, with input-driven states!
% now include opto input
% Note that now the loaded data includes opto input, and the log-likelihood
% as well as M-step functions are different

%% load pre-processed Data structure with the proper fields
% run or load the Data structure
% load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_nai_fit.mat')
% [Data] = track2data(Tracks, cand, M);

% extracting typical data
[xxf, yyf, alltrials, time] = data2xy(Data);
% opto
allopto = extractfield(Data, 'opto');
% displacement
alldis = extractfield(Data, 'dis');

% creating new xx,yy for fitting
yyf = [yyf; alldis];
xxf = [xxf; allopto];

maskf = true(1,length(yyf));
maskf = alltrials;  %((alltrials)==0) = false;%

%%% pick a window
wind = 1:300000;
xx = xxf(:,wind);
yy = yyf(:,wind);
mask = maskf(wind);

%% Observation and input
% Set parameters: transition matrix and emission matrix
nStates = 2; % number of latent states
nX = 17 + 4;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_mVMG_opto;  % log-likelihood function, with opto input
loglitrans = @logli_trans_opto;  % log-likelihood of soft-max state transitions, with opto driven

% Set transition matrix by sampling from Dirichlet distr
alpha_diag = 30;  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1

% sticky priors
alpha = 1.;  % Dirichlet shape parameter as a prior
kappa = .5;  % upweighting self-transition for stickiness

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
% wts0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.]; 
wts0 = rand(nY,nX,nStates); % parameters for the mixture-VonMesis behavioral model
wts0(1,:,1) = [50,  randn(1,nB)*1,  10, 25,  10, 25,  5,   1.  0.1 0  1 1 0 0  randn(1,nB)*.1,]; %[20, randn(1,nB)*10, 10, 25, 10, 25, 1, 1.]; %single mGLM
wts0(1,:,2) = [10,  randn(1,nB)*1, -10, 25, -10, 25, 20,  .5   0.1 0  1 1 0 0  randn(1,nB)*.1,];
% wts0(1,:,3) = [10,  randn(1,nB)*10, -10, 25, -10, 25, 20,.5 0.1 0  1 1 0 0];
% wts0(1,:,4) = [10,  randn(1,nB)*10, -20, 25, -20, 25, 20,.5 0.1 0  1 1 0 0];

%%% transition kernels
alpha_ij = randn(nStates,nStates, nB)*.1;  % state x state x weights
alpha_opto = randn(nStates,nStates, nB)*.1;  % state x state x weights
lk = size(cosBasis,1);
K_ij = zeros(nStates, nStates, lk);
Kij_opto = K_ij*1;
for ii = 1:nStates
    for jj = 1:nStates
        if ii~=jj
            K_ij(ii,jj,:) = fliplr(squeeze(alpha_ij(ii,jj,:))'*cosBasis');  % build state-transition kernels, without self-kernels
            Kij_opto(ii,jj,:) = fliplr(squeeze(alpha_opto(ii,jj,:))'*cosBasis');
        else
            alpha_ij(ii,jj) = alpha_ij(ii,jj) *0;
        end
    end
end
wts_state0 = alpha_ij;
w_state_opto0 = alpha_opto;

% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'wts_state',wts_state0,'w_state_opto',w_state_opto0,'loglifun',loglifun,'loglitrans',loglitrans,'basis',cosBasis,'lambda',zeros(1,nStates)+0.0);

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
    mmhat = runMstep_mGVM_opto(mmhat, xx, yy, gams, mask);
    % Update for input-driven transitions
    mmhat = runMstep_state_opto(mmhat, xx, yy, xis, mask); 
    
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
wind_test = max(wind):length(alldis);%1:max(wind); %
yy = yyf(:,wind_test);
xx = xxf(:,wind_test);
mask = maskf(wind_test);
[logp_test,gams_test,xis_test,xisum_test] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
logp_test/length(wind_test)

%% sampling tests
n_samp = 10;  % repeats
l_samp = 20000;  % length
test_lls = zeros(1, n_samp);
for nn = 1:n_samp
    randomInteger = randi([max(wind), length(yyf)-l_samp]);%randi([1,min(wind)]);%
    wind_i = randomInteger:randomInteger+l_samp;
    yy = yyf(:, wind_i);
    xx = xxf(:, wind_i);
    mask = maskf(wind_i);
    [logp_test,gams_test,xis_test, xisum_test] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
    test_lls(nn) = logp_test/length(wind_i);
end
test_lls
