% mGLMHMM
%%% serve as ground-truth to varify the inference procedure
%% generate hidden states
lt = 50000;  % length of simulation
nStates = 2;  % two hidden states for now
% Set transition matrix by sampling from Dirichlet distr
alpha_diag = 25;  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
A0 = [0.99,0.01; 0.02,0.98];
A0 % Markov transition
mc = dtmc(A0);
X = simulate(mc, lt);
alpha = 1;

%% define variables
nB = 4;  % number of basis function for the kernel
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3); % basis function
%%% true params
beta = 1; % nonlinear parameter
A = 0.5;  % max probability of turning
alpha_h1 = [-4:-1]*.01;  % angle history kernel coefficient
alpha_dc1 = [1,4,-2,-1]*5;  % dC kernel coefficient
alpha_dcp1 = [-4:-1]*.01;  % dCp kernel coefficient
base = 0;  %baseline
kappa_turn1 = 5;  % turning angle variance
kappa_wv1 = 20;  % weather-vaning angle variance

alpha_h2 = [-4:-1]*.1;  % angle history kernel coefficient
alpha_dc2 = [1,4,-2,-1]*1;  % dC kernel coefficient
alpha_dcp2 = [-4:-1]*-.01;  % dCp kernel coefficient
kappa_turn2 = 20;  % turning angle variance
kappa_wv2 = 5;  % weather-vaning angle variance

%%% construct as kernel
K_h1 = fliplr(alpha_h1*cosBasis');  % dth kernel
K_dc1 = fliplr(alpha_dc1*cosBasis');  % dC kernel
K_dcp1 = fliplr(alpha_dcp1*cosBasis');  % dCp kernel
K_h2 = fliplr(alpha_h2*cosBasis');  % dth kernel
K_dc2 = fliplr(alpha_dc2*cosBasis');  % dC kernel
K_dcp2 = fliplr(alpha_dcp2*cosBasis');  % dCp kernel

%% generate data
r2d = 1;%180/pi;
lc = 10;  %length of smoothing
dC = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dC stimulus vector
dCp = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dCp stimulus vector
dth = zeros(1,lt);
turns = zeros(1,lt);
F = dth*0;
pad = length(K_h1);
for tt=pad+1:lt
    if X(tt) == 1
        K_h = K_h1; K_dc = K_dc1; K_dcp = K_dcp1; kappa_turn = kappa_turn1; kappa_wv = kappa_wv1;
    elseif X(tt) == 2
        K_h = K_h2; K_dc = K_dc2; K_dcp = K_dcp2; kappa_turn = kappa_turn2; kappa_wv = kappa_wv2;
    end
    F(tt) = dC(tt-pad+1:tt)*K_dc' + abs(dth(tt-pad:tt-1))*r2d*K_h';  % linear filtering
    turns(tt) = choice(NL(F(tt)+base,A));  % nonlinearity and binary choice
    dth(tt) = turns(tt)*circ_vmrnd(pi,kappa_turn,1) + (1-turns(tt))*circ_vmrnd(dCp(tt-pad+1:tt)*K_dcp',kappa_wv,1);  % angle drawn from mixture of von Mesis
end

allas = dth*r2d;
alldCs = dC;
alldCps = dCp;

figure; hist(allas,100)
%% load training data
wind = 1:lt;
yy = allas(wind);
xx = [alldCs(wind); 
      alldCps(wind)];

mask = true(1,length(yy));
% mask(isnan(alltrials(wind))) = false;%

%% call inference procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Observation and input
% Set parameters: transition matrix and emission matrix
nStates = 2; % number of latent states
nX = 14;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_GTmVM;  % ground truth model log-likelihood function

% Set transition matrix by sampling from Dirichlet distr
% alpha_diag = 25;  % concentration param added to diagonal (higher makes more diagonal-dominant)
% alpha_full = 10;  % concentration param for other entries (higher makes more uniform)
% G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
% A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
A0 = A0;

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
% wts0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.]; 
wts0 = rand(nY,nX,nStates); % parameters for the mixture-VonMesis behavioral model
wts0(1,:,1) = [alpha_h1, alpha_dc1, alpha_dcp1, kappa_turn1, kappa_wv1]; %single mGLM
wts0(1,:,2) = [alpha_h2, alpha_dc2, alpha_dcp2, kappa_turn2, kappa_wv2];
% wts0 = wts0 + rand(nY, nX, nStates);

% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'loglifun',loglifun,'basis',cosBasis,'lambda',.0, 'Mstepfun',@runMstep_GTmVM);

%% debug
% [logli] = logli_GTmVM(mmhat, xx, yy, mask);

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
    
    % Update transition matrix
    mmhat.A = (alpha-1 + xisum) ./ (nStates*(alpha-1) + sum(xisum,2)); % normalize each row to sum to 1
    
    % Update model params
    mmhat = runMstep_GTmVM(mmhat, xx, yy, gams, mask);
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


%% compare to gound-truth
stateK = 2;
x = squeeze(mmhat.wts(:,:,stateK));
alpha_h_ = x(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
alpha_dc_ = x(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
alpha_dcp_ = x(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
kappa_turn_ = x(13);     % vairance of the sharp turn von Mises
kappa_wv_ = x(14);       % variance of von Mises

K_dcp_rec = alpha_dcp_*cosBasis';
K_h_rec = alpha_h_*cosBasis';
K_dc_rec = alpha_dc_*cosBasis';

figure()
subplot(131)
plot(fliplr(K_h_rec)); hold on; plot(K_h1,'--'); plot(K_h2,'--');
subplot(132)
plot(fliplr(K_dc_rec)); hold on; plot(K_dc1,'--'); plot(K_dc2,'--');
subplot(133)
plot(fliplr(K_dcp_rec)); hold on; plot(K_dcp1,'--'); plot(K_dcp2,'--');

[aa,bb] = max( gams ,[], 1 );
figure;
plot(allas(:)); hold on
% plot(smooth((bb(wind)-1)*100,10))
plot((bb(:)-1)*3)

%% functions
%%% nonlinear function
function [P] = NL(F, A)
    beta = 1;
    P = A./(1+exp(-beta*F));
end

%%% stochastic choice
function [b] = choice(P)
    pp = rand();
    if pp<P
        b = 1;
    else
        b = 0;
    end
end

%% inference ll and M-step
function [logli] = logli_GTmVM(mm, xx, yy, mask)

% Compute log-likelihood term under a mixture of von Mesis model
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%      .wts   [1 len(param) K] - per-state parameters for the model
%      .basis [len(alpha) len(kernel)] - basis functions used for the kernel
%      .lambda -scalar for regularization of the logistic fit
%    xx [2 T] - inputs (time series of dc,dcp, and dth)
%    yy [1 T] - outputs (dth angles)
%
% Output
% ------
%  logpy [T K] - loglikelihood for each observation for each state
    
    %%% loading parameters
    THETA = mm.wts;
    Basis = mm.basis;
    lambda = mm.lambda;
    dth = yy;
    dc = xx(1,:);
    dcp = xx(2,:);
    
    K = size(THETA,3);
    lt = length(yy);
    logli = zeros(lt, K);
    for k = 1:K
        %%% Assume we parameterize in such way first
        alpha_h = THETA(1,1:4,k);       % kernel for dth angle history kernel (weights on kerenl basis)
        alpha_dc = THETA(1,5:8,k);      % kernel for dC transitional concentration difference (weights on kerenl basis)
        alpha_dcp = THETA(1,9:12,k);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
        kappa_turn = THETA(1,13,k)^0.5;     % vairance of the sharp turn von Mises
        kappa_wv = THETA(1,14,k)^0.5;       % variance of von Mises

        %%% construct as kernel
        K_h = (alpha_h * Basis');  % dth kernel
        K_dc = (alpha_dc * Basis');  % dC kernel
        K_dcp = (alpha_dcp * Basis');  % dCp kernel

        %%% turning decision
        filt_dth = conv_kernel(abs(dth(1:end-1)), K_h);
        filt_dc = conv_kernel(dc(2:end), K_dc);
        P = NL(filt_dc + filt_dth, 0.5); %1 ./ (1 + exp( -(filt_dc + filt_dth + 0))) +0;  %sigmoid(A_,B_,dc); 

        %%% weathervaning part
        d2r = 1;%pi/180;
        C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
        filt_dcp = conv_kernel(dcp(2:end), K_dcp);
        VM = C * exp(kappa_wv^2*cos(( filt_dcp - dth(2:end) )*d2r));  %von Mises distribution

        %%% turning analge model
        VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)

        %%% marginal probability
        marginalP = (1-P).*VM + VM_turn.*P;
        logli(2:end,k) = ( mask(2:end).* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));% + 0.1*sum((E_ - 0).^2) + 0*C_^2);  % adding slope l2 regularization
    end
end

function mm = runMstep_GTmVM(mm, xx, yy, gams, mask)
% mm = runMstep_LinGauss(mm,xx,yy,gams)
%
% Run m-step updates for Gaussian observation model
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%        .wts  [1 K] - per-state slopes
%        .vars [1 K] - per-state variances
%    xx [d T] - input (design matrix)
%    yy [1 T] - outputs
%  gams [K T] - log marginal sate probabilities given data
%
% Output
% ------
%  mmnew - new model struct

% normalize the gammas to sum to 1
gamnrm = gams./(sum(gams,2)+1e-10);  
nStates = size(gams,1);

%%% loading parameters
Basis = mm.basis;
lambda = mm.lambda;
dth = yy;
dc = xx(1,:);
dcp = xx(2,:);
nB = size(Basis,2);

% setup fmincon
% lfun = @(x)nll_mVM(x, dth, dcp, dc, gams, Basis, lambda, mask);
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [ones(1,12)*-10, 0, 0.];
UB = [ones(1,12)*10, 20, 20.];
% prs0 = rand(1,10);
% prs0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.] ;
% prs0 = [9.9763  -0.5343  -0.0776   0.1238  -0.0529   0.5335   7.7254  367.3817  0.1990  1.0000  0.1000]; %single mGLM
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

for jj = 1:nStates
    lfun = @(x)nll_mVM(x, dth, dcp, dc, gamnrm(jj,:), Basis, lambda, mask);
%     prs0 = prs0 + prs0.*randn(1,length(UB))*0.5;
    prs0 = mm.wts(:,:,jj); %+ mm.wts(:,:,jj).*(2*(rand(1,length(UB))-0.5))*0.5;  %from last time!
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
%     x = fminunc(lfun, prs0);
    mm.wts(:,:,jj) = x; % weighted optimization
end


end

function [nll] = nll_mVM(THETA, dth, dcp, dc, gams, Basis, lambda, mask)
    
    %%% Assume we parameterize in such way first
    alpha_h = THETA(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
    alpha_dc = THETA(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
    alpha_dcp = THETA(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
    kappa_turn = THETA(13)^0.5;     % vairance of the sharp turn von Mises
    kappa_wv = THETA(14)^0.5;       % variance of von Mises

    %%% construct as kernel
    K_h = (alpha_h * Basis');  % dth kernel
    K_dc = (alpha_dc * Basis');  % dC kernel
    K_dcp = (alpha_dcp * Basis');  % dCp kernel
    
    filt_dth = conv_kernel(abs(dth(1:end-1)), K_h);
    filt_dc = conv_kernel(dc(2:end), K_dc);
    P = NL(filt_dc + filt_dth, 0.5);  %1 ./ (1 + exp( -(filt_dc + filt_dth + 0))) +0;  %sigmoid(A_,B_,dc); 

    %%% weathervaning part
    d2r = 1;%pi/180;
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    filt_dcp = conv_kernel(dcp(2:end), K_dcp);
    VM = C * exp(kappa_wv^2*cos(( filt_dcp - dth(2:end) )*d2r));  %von Mises distribution

    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)

    %%% marginal probability
    marginalP = (1-P).*VM + VM_turn.*P;
    nll = -nansum( mask(2:end) .* gams(2:end) .* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));% + 0.1*sum((
end
