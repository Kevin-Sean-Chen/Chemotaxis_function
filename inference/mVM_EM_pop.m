%% mVM_EM_pop
%% load some test data
% assuming we have conditioned tracks and analysis in the structure 'Data' and concatenated data all...
% tracks
ntracks = length(Data);

% make test and train data tracks
indices = crossvalind('Kfold',[1:ntracks], 3); 
test_set = (indices==ki);
train_set = ~test_set;
% concatenated data vectors

[xx, yy, mask] = data2xy(Data(train_set));

% wind = 1:length(allas);
% yy = allas(wind);
% xx = [alldC(wind); 
%       alldcp(wind)];
% mask = true(1,length(yy));
% mask(isnan(alltrials(wind))) = false;

%% Observation and input
% Set parameters: transition matrix and emission matrix
nStates = 2; % number of latent states
nX = 11;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_mVM;  % log-likelihood function

% Set transition matrix by sampling from Dirichlet distr
alpha_diag = 30;  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
% A0 = [0.99,0.01; 0.01,0.99];

alpha = 1.;  % Dirichlet shape parameter as a prior
kappa = 0.5;  % upweighting self-transition for stickiness

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
% wts0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.]; 
wts0 = rand(nY,nX,nStates); % parameters for the mixture-VonMesis behavioral model
wts0(1,:,1) = [20, randn(1,nB)*10, -50, 25, 10, 25, 5, 1.]; %[20, randn(1,nB)*10, 10, 25, 10, 25, 1, 1.]; %single mGLM
wts0(1,:,2) = [10,  randn(1,nB)*10, -50, 25, -10, 25, 20,.5];
% wts0(1,:,3) = [1,  randn(1,nB)*10, -10, 25, -10, 25, 20,.5];
%[9.9763  -0.5343  -0.0776   0.1238  -0.0529   0.5335   7.7254  367.3817  0.1990  1.0000  0.1000];


% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'loglifun',loglifun,'basis',cosBasis,'lambda',[0.1,.0]);

%% Set up variables for EM with tracks-based likelihood
maxiter = 50;
EMdisplay = 2;
logpTrace = zeros(maxiter,1); % trace of log-likelihood
dlogp = inf; % change in log-likelihood
logpPrev = -inf; % prev value of log-likelihood
jj = 1; % counter

train_id = find(train_set==1);
while (jj <= maxiter) && (dlogp>1e-3)
    
    gams = [];
    xisum = zeros(nStates,nStates);
    %%% load tracks
    for k = 1:length(train_id) %ntracks
        
        % load track observation
%         y_k = Data(k).dth;
%         x_k = [Data(k).dc;
%               Data(k).dcp];
%         mask_k = Data(k).mask;
        [x_k, y_k, mask_k] = data2xy(Data(train_id(k)));
        mask_k(isnan(mask_k)) = 0;
        
        % --- run E-step ---
        [logp,gams_k,xisum_k] = runFB_GLMHMM(mmhat,x_k,y_k,mask_k); % run forward-backward
        logpTrace(jj) = logpTrace(jj) + logp;
        
        % append or adding
        gams = [gams gams_k];
        xisum = xisum + xisum_k;
    end
    
    % --- run M-step ---
    % Update transition matrix
%     mmhat.A = (alpha-1 + xisum) ./ (nStates*(alpha-1) + sum(xisum,2)); % normalize each row to sum to 1
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

%% testing
[xx_test, yy_test, mask_test] = data2xy(Data(test_set));
% wind_test = max(wind):length(allas);
% yy = allas(wind_test);
% xx = [alldC(wind_test); 
%       alldcp(wind_test)];
% mask = true(1,length(yy));
% mask(isnan(alltrials(wind_test))) = false;%
[logp_test,gams_test,xisum_test] = runFB_GLMHMM(mmhat,xx_test,yy_test,mask_test);
logp_test/length(mask_test)
