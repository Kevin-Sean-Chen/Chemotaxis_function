% mGLMHMM_gamma
%%% serve as ground-truth to varify the inference procedure
%%% this version contains speed distribution with gamma, together with the
%%% angular mixture model; now also a mixture of gamma!

%% generate hidden states
lt = 50000;  % length of simulation
nStates = 2;  % two hidden states for now
% Set transition matrix by sampling from Dirichlet distr
alpha_diag = 25;  % concentration param added to diagonal (higher makes more diagonal-dominant)
alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
A0 = [0.99,0.01; 0.01,0.99];
A0 % Markov transition
mc = dtmc(A0);
X = simulate(mc, lt);
alpha = 1.;
kappa = 0.;

%% define variables
nB = 4;  % number of basis function for the kernel
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3); % basis function
%%% true params
beta = 1; % nonlinear parameter
alpha_h1 = [-4:-1]*.05;  % angle history kernel coefficient
alpha_dc1 = [1,4,-2,-1]*2;  % dC kernel coefficient
alpha_dcp1 = [-4:-1]*.01;  % dCp kernel coefficient
base = 0;  %baseline
A = 0.5;  % max turning probability
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

%% run distributions
alphas = [1, 10;
          2,  5];  % shape (state x mixture)
betas = [1, 2;
         5, 7];  % scale
w_gam = [0.1, 0.9];  % weight on one another

%% generate data
r2d = 1;%180/pi;
lc = 10;  %length of smoothing
dC = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dC stimulus vector
dCp = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dCp stimulus vector
dth = zeros(1,lt);
dv = zeros(1,lt);
turns = zeros(1,lt);
F = dth*0;
pad = length(K_h1);
for tt=pad+1:lt
    %%% loading states
    if X(tt) == 1
        K_h = K_h1; K_dc = K_dc1; K_dcp = K_dcp1; kappa_turn = kappa_turn1; kappa_wv = kappa_wv1; aa1 = alphas(1,1); bb1 = betas(1,1); aa2 = alphas(1,2); bb2 = betas(1,2); wg=w_gam(1);
    elseif X(tt) == 2
        K_h = K_h2; K_dc = K_dc2; K_dcp = K_dcp2; kappa_turn = kappa_turn2; kappa_wv = kappa_wv2; aa1 = alphas(2,1); bb1 = betas(2,1); aa2 = alphas(2,2); bb2 = betas(2,2); wg=w_gam(2);
    end
    %%% angular decisions
    F(tt) = dC(tt-pad+1:tt)*K_dc' + abs(dth(tt-pad:tt-1))*r2d*K_h';  % linear filtering
    turns(tt) = choice(NL(F(tt)+base,A));  % nonlinearity and binary choice
    dth(tt) = turns(tt)*circ_vmrnd(pi,kappa_turn,1) + (1-turns(tt))*circ_vmrnd(dCp(tt-pad+1:tt)*K_dcp',kappa_wv,1);  % angle drawn from mixture of von Mesis
    %%% mixture velocity
    if rand()<wg
        dv(tt) = gamrnd(aa1, bb1, 1, 1);  % gamma random variable
    else
        dv(tt) = gamrnd(aa2, bb2, 1, 1);  % gamma random variable
    end
end

allas = dth*r2d;
alldis = dv;
alldCs = dC;
alldCps = dCp;

figure; subplot(121); hist(allas,100); subplot(122); hist(alldis, 100);

%% load training data
wind = 1:lt;
yy = [allas(wind);
      alldis(wind)];
xx = [alldCs(wind); 
      alldCps(wind)];

mask = true(1,length(yy));
mask(1:50) = false; %%% important!!!! log(v) would have undefined numbers without masking, be careful!
% mask(isnan(alltrials(wind))) = false;%

%% call inference procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Observation and input
% Set parameters: transition matrix and emission matrix
nStates = 2; % number of latent states
n_mix_gam = 2;  % specify number of mixtures for velocity
nX = 14;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_GTmVM;  % ground truth model log-likelihood function

% % Set transition matrix by sampling from Dirichlet distr
% alpha_diag = 30;  % concentration param added to diagonal (higher makes more diagonal-dominant)
% alpha_full = 5;  % concentration param for other entries (higher makes more uniform)
% G = gamrnd(alpha_full*ones(nStates) + alpha_diag*eye(nStates),1); % sample gamma random variables
% A0 = G./repmat(sum(G,2),1,nStates); % normalize so rows sum to 1
% % A0 = A0;

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
% wts0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.]; 
wts0 = rand(nY,nX,nStates); % parameters for the mixture-VonMesis behavioral model
wts0(1,:,1) = [alpha_h1, alpha_dc1, alpha_dcp1, kappa_turn1, kappa_wv1];%, alphas(1,1), betas(1,1), alphas(1,2), betas(1,2), w_gam(1)]; %single mGLM
wts0(1,:,2) = [alpha_h2, alpha_dc2, alpha_dcp2, kappa_turn2, kappa_wv2];%, alphas(2,1), betas(2,1), alphas(2,2), betas(2,2), w_gam(2)];
wts0 = wts0 + wts0.*rand(nY,nX,nStates)*0.1; % small perturbation

mix_gamm0 = zeros(nStates, n_mix_gam, 3);  % three for shape, scale, and weight (probably easier for joint likelihood computation)
mix_gamm0(1,1,:) = [alphas(1,1), betas(1,1), w_gam(1)];
mix_gamm0(1,2,:) = [alphas(1,2), betas(1,2), 1-w_gam(1)];
mix_gamm0(2,1,:) = [alphas(2,1), betas(2,1), w_gam(2)];
mix_gamm0(2,2,:) = [alphas(2,2), betas(2,2), 1-w_gam(2)];
% mix_gamm0(:,1) = [alphas(1,1), betas(1,1), alphas(1,2), betas(1,2), w_gam(1) ];


% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'mix_gamma',mix_gamm0,'loglifun',loglifun,'basis',cosBasis,'lambda',.0, 'Mstepfun',@runMstep_GT_mVM_mG);

%%
% [logp,gams,xisum] = runFB_GLMHMM(mmhat,xx,yy,mask)
% test = runMstep_GTmVM(mmhat, xx, yy, gams, mask)
% test2 = logli_GTmVM(mmhat, xx, yy, mask)
% -nLL_gamma([10, 2], dv, mask)

%% Set up variables for EM
maxiter = 50;
EMdisplay = 2;
logpTrace = zeros(maxiter,1); % trace of log-likelihood
dlogp = inf; % change in log-likelihood
logpPrev = -inf; % prev value of log-likelihood
jj = 1; % counter
% lamb_dk = ones(nStates, n_mix_gam);
% lamb_dk = lamb_dk./sum(lamb_dk,2);  % P(mixture|state), for constrained hierarchical model
lamb_dk = [w_gam(1)  1-w_gam(1);
           w_gam(2)  1-w_gam(2)];   % initualize with true weights on the mixtures (state x mixture)

while (jj <= maxiter) && (dlogp>1e-3)
    
    % --- run E step  -------
    [logp,gams,xisum] = runFB_GLMHMM(mmhat,xx,yy,mask); % run forward-backward
    logpTrace(jj) = logp;
   
    [delta_tdk] = run_E_mix_gamma(mmhat, xx, yy, mask, lamb_dk);  % E-step for mixture gamma: P(mixture | state, observed)
    
    % --- run M step  -------
    
    % Update transition matrix
    mmhat.A = (alpha-1 + xisum) ./ (nStates*(alpha-1) + sum(xisum,2)); % normalize each row to sum to 1
    % Update transition matrix (with stickiness)
%     normed_xisum = xisum ./ sum(xisum,2);
%     unormedA = kappa*eye(nStates) + (alpha-1)*ones(nStates,nStates) + normed_xisum;
%     mmhat.A = unormedA ./ sum(unormedA,2);
    
    % Update model params
    mmhat = runMstep_GTmVM(mmhat, xx, yy, gams, mask);
    [mmhat, lamb_dk] = runMstep_mix_gamma(mmhat, xx, yy, gams, delta_tdk, mask);
    % jointed update seems better?
%     [mmhat, lamb_dk] = runMstep_GT_mVM_mG(mmhat, xx, yy, gams, delta_tdk, mask);
     
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
figure()
for kk = 1:stateK
x = squeeze(mmhat.wts(:,:,kk));
alpha_h_ = x(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
alpha_dc_ = x(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
alpha_dcp_ = x(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
kappa_turn_ = x(13);     % vairance of the sharp turn von Mises
kappa_wv_ = x(14);       % variance of von Mises

K_dcp_rec = alpha_dcp_*cosBasis';
K_h_rec = alpha_h_*cosBasis';
K_dc_rec = alpha_dc_*cosBasis';

subplot(131)
plot(fliplr(K_h_rec)); hold on;
subplot(132)
plot(fliplr(K_dc_rec)); hold on;
subplot(133)
plot(fliplr(K_dcp_rec)); hold on;

end
subplot(131)
plot(K_h1,'--'); plot(K_h2,'--');
subplot(132)
plot(K_dc1,'--'); plot(K_dc2,'--');
subplot(133)
plot(K_dcp1,'--'); plot(K_dcp2,'--');

[aa,bb] = max( gams ,[], 1 );
figure;
subplot(121)
plot(allas(:)); hold on
% plot(smooth((bb(wind)-1)*100,10))
plot((bb(:)-1)*3)
subplot(122)
plot(alldis); hold on; plot(bb*10)

%% evaluate density
figure;
kk = 1;
x = squeeze(mmhat.wts(:,:,kk));
alpha_h_ = x(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
alpha_dc_ = x(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
alpha_dcp_ = x(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
kappa_turn_ = x(13);     % vairance of the sharp turn von Mises
kappa_wv_ = x(14);       % variance of von Mises

K_dcp_rec = alpha_dcp_*cosBasis';
K_h_rec = alpha_h_*cosBasis';
K_dc_rec = alpha_dc_*cosBasis';

filt_ddc = conv_kernel(dC, fliplr(K_dc_rec));
filt_dth = conv_kernel(abs(dth), fliplr(K_h_rec));
dc_dth = filt_ddc + filt_dth;
Pturns = NL(dc_dth,A);
n_brw = sum(Pturns)*1;
n_wv = sum(1-Pturns);
p_z = n_brw + n_wv;
p_brw = n_brw/p_z;
p_wv = n_wv/p_z;
filt_dcp = conv_kernel(dCp, fliplr(K_dcp_rec));
figure;
pos = find(X(1:end-1)==kk);
nbins = 100;
hh = histogram(dth(pos), nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1);
% [aa,bb] = hist((dth - filt_dcp)*1 , 1000);
scal = sum(1/(2*pi*besseli(0,kappa_wv_)) * exp(kappa_wv_*cos( bb )) * p_wv + (1/(2*pi*besseli(0,kappa_turn_)) * exp(kappa_turn_*cos( bb-pi ))) * p_brw);
plot( bb, 1/(2*pi*besseli(0,kappa_wv_)) * exp(kappa_wv_*cos( bb )) * p_wv / scal , 'b'); hold on
plot( bb, ((1/(2*pi*besseli(0,kappa_turn_)) * exp(kappa_turn_*cos( bb-pi ))) * p_brw)/ scal ,'r'); xlim([-pi,pi])
title('von Mises for \delta C^{\perp}')

%% density for speed
figure
for kk = 1:nStates
    x = reshape(mmhat.mix_gamma(kk,:,:),1, []);
%     x = reshape(mix_gamm0(kk,:,:),1, []);
    gam_shape_1 = x(1);      % shape parameter for gamma1
    gam_shape_2 = x(2);      % shape parameter for gamma2
    gam_scale_1 = x(3);      % scale parameter for gamma1
    gam_scale_2 = x(4);      % scale parameter for gamma2
    gam_weight1 = x(5);       % weight on gamma mixture1
    gam_weight2 = x(6);       % weight on gamma mixture2
    
    pos = find(X(1:end-1)==kk);
    hh = histogram(dv(pos), nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
    bb = hh.BinEdges(1:end-1);
    test = hh.Values;
    scal = 1;%length(find(X==kk))/ length(X);
    gam_p = gam_weight1*gampdf(bb, gam_shape_1, gam_scale_1) + (gam_weight2)*gampdf(bb, gam_shape_2, gam_scale_2);
    plot(bb, scal*gam_p/sum(gam_p))
end

%% debugging
%%
% gamnrm = gams./(sum(gams,2)+1e-8);  
% nStates = size(gams,1);
% nMix = size(delta_tdk,3);
% 
% %%% loading parameters
% Basis = mmhat.basis;
% lambda = mmhat.lambda;
% dth = yy(1,:);
% dv = yy(2,:);
% dc = xx(1,:);
% dcp = xx(2,:);
% jj = 1;
% x = mmhat.wts(:,:,jj);
% % nll_mVM_mG(x, dth, dcp, dc, dv, gamnrm(jj,:), squeeze(delta_tdk(:,jj,:)), Basis, lambda, mask, mmhat, jj)
% nll_mi_gamma( dv, gamnrm(jj,:), squeeze(delta_tdk(:,jj,:)), mask, mmhat, jj)

%%
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
    THETA_gam = mm.mix_gamma;
    Basis = mm.basis;
    lambda = mm.lambda;
    dth = yy(1,:);
    dv = yy(2,:);
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
        gamm_shape1 = THETA_gam(k,1,1);     % shape of gamma1
        gamm_scale1 = THETA_gam(k,1,2);     % scale of gamma1
        gamm_weight1 = THETA_gam(k,1,3);    % weight on gamma mixture1
        gamm_shape2 = THETA_gam(k,2,1);     % shape of gamma2
        gamm_scale2 = THETA_gam(k,2,2);     % scale of gamma2
        gamm_weight2 = THETA_gam(k,2,3);    % weight on gamma mixture2
        
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
        pos = (~isinf(log(dv)));
%         ll_gamm = ((gamm_shape1-1)*log(dv.*pos) - dv/gamm_scale1 - gamm_shape1*log(gamm_scale1) - gammaln(gamm_shape1)) + log(gamm_weight1) + ...
%                   ((gamm_shape2-1)*log(dv.*pos) - dv/gamm_scale2 - gamm_shape2*log(gamm_scale2) - gammaln(gamm_shape2)) + log(gamm_weight2);
        %%% log of weighted pro
        ll_gamm = log(gamm_weight1*gampdf(dv, gamm_shape1, gamm_scale1) + gamm_weight2*gampdf(dv, gamm_shape2, gamm_scale2)).*pos;
        logli(2:end,k) = ( mask(2:end).* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2))  + mask(2:end).*ll_gamm(2:end);
        % - nLL_gamma([gamm_shape, gamm_scale], dv, mask);
        %+ log(gampdf(dv,gamm_shape, gamm_scale));
        %
    end
end

function [mm, lamb_dk] = runMstep_GT_mVM_mG(mm, xx, yy, gams, delta_tdk, mask)

gamnrm = gams./(sum(gams,2)+1e-8);  
nStates = size(gams,1);
nMix = size(delta_tdk,3);

%%% loading parameters
Basis = mm.basis;
lambda = mm.lambda;
dth = yy(1,:);
dv = yy(2,:);
dc = xx(1,:);
dcp = xx(2,:);

opts = optimset('display','iter');
LB = [ones(1,12)*-10, 0, 0.,  zeros(1,4),  0  0];
UB = [ones(1,12)*10, 20, 20., zeros(1,4)+100, 1,1];

for jj = 1:nStates
    
    %%% optimize weights for angles
    lfun = @(x)nll_mVM_mG(x, dth, dcp, dc, dv, gamnrm(jj,:), squeeze(delta_tdk(:,jj,:)), Basis, lambda, mask, mm, jj);
%     prs0 = [mm.wts(:,:,jj)   reshape(mm.mix_gamma(jj,:,:),1,[]) ]; % combining parameters!
    prs0 = [mm.wts(:,:,jj) ]; 
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    mm.wts(:,:,jj) = x;%(1:14); % weighted optimization
    
    %%% optimize parameters for gamma mixture
%     mm.mix_gamma(jj,:,:) = reshape(x(15:end), nMix,[]);
end

d_states = size(mm.mix_gamma,2);
k_mixtures = size(mm.wts,3);
%%% update lambda_dk: P(mixture | state)
lamb_dk = zeros(d_states, k_mixtures);
for d = 1:d_states
    for k = 1:k_mixtures
        lamb_dk(d,k) = gams(d,:)*delta_tdk(:,d,k);
    end
end
lamb_dk = lamb_dk ./ (nansum(gams,2)+1e-8);
% mm.mix_gamma(:,:,end) = lamb_dk;  % updating the weights directly!

end

%%

function mm = runMstep_GTmVM(mm, xx, yy, gams, mask)
% mm = runMstep_LinGauss(mm,xx,yy,gams)
%
% Run m-step updates for Gaussian observation model
%
% Inputs
% ------
%    mm [struct] - model structure with params 
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
gamnrm = gams./(sum(gams,2)+1e-8);  
nStates = size(gams,1);

%%% loading parameters
Basis = mm.basis;
lambda = mm.lambda;
dth = yy(1,:);
dv = yy(2,:);
dc = xx(1,:);
dcp = xx(2,:);
nB = size(Basis,2);

% setup fmincon
% lfun = @(x)nll_mVM(x, dth, dcp, dc, gams, Basis, lambda, mask);
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [ones(1,12)*-10, 0, 0.,  0.,  0.  0.,  0., 0];
UB = [ones(1,12)*10, 20, 20., 100, 100, 100, 100, 1];
% prs0 = rand(1,10);
% prs0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.] ;
% prs0 = [9.9763  -0.5343  -0.0776   0.1238  -0.0529   0.5335   7.7254  367.3817  0.1990  1.0000  0.1000]; %single mGLM
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

for jj = 1:nStates
    
    %%% optimize weights for angles
    lfun = @(x)nll_mVM(x, dth, dcp, dc, dv, gamnrm(jj,:), Basis, lambda, mask);
%     prs0 = prs0 + prs0.*randn(1,length(UB))*0.5;
    prs0 = mm.wts(:,:,jj); %+ mm.wts(:,:,jj).*(2*(rand(1,length(UB))-0.5))*0.5;  %from last time!
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
%     x = fminunc(lfun, prs0);
    mm.wts(:,:,jj) = x; % weighted optimization
    
%     %%% optimize parameters for gamma mixture
%     prs0_gam = mm.mix_gamma(:,jj);
%     xg = mix_gam_EM(prs0_gam, dv, gamnrm(jj,:), );
%     mm.mix_gamma(:,jj) = xg;
end


end

function [delta_tdk] = run_E_mix_gamma(mm, xx, yy, mask, lamb_dk)
    dv = yy(2,:);  % load observable
    d_states = size(mm.mix_gamma,1);
    k_mixtures = size(mm.mix_gamma,2);
    lt = length(xx);
    delta_tdk = ones(lt, d_states, k_mixtures);
    mix_gamma_params = mm.mix_gamma;
    pos = (~isnan((dv))).*(~isinf((dv))).*mask;
    % mix_gamm0(:,1) = [alphas(1,1), betas(1,1), alphas(1,2), betas(1,2), w_gam(1) ];
    for d = 1:d_states
        for k = 1:k_mixtures
            delta_tdk(:, d, k) = pos.*  mix_gamma_params(d, k, 3) .* gampdf(dv, mix_gamma_params(d, k, 1), mix_gamma_params(d, k, 2));
        end
    end
    delta_tdk = delta_tdk ./ (nansum(delta_tdk, 3)+1e-8);
end

function [mm, lamb_dk] = runMstep_mix_gamma(mm, xx, yy, gams, delta_tdk, mask);
    dv = yy(2,:);  % load observable
    d_states = size(mm.mix_gamma,2);
    k_mixtures = size(mm.wts,3);
    %%% update lambda_dk: P(mixture | state)
    lamb_dk = zeros(d_states, k_mixtures);
    mask_ = (~isnan(delta_tdk));
    for d = 1:d_states
        for k = 1:k_mixtures
            lamb_dk(d,k) = gams(d,:)*(delta_tdk(:,d,k).*mask_(:,d,k));
        end
    end
    lamb_dk = lamb_dk ./ (nansum(gams,2)+1e-8);
    mm.mix_gamma(:,:,end) = lamb_dk;  % updating the weights directly!
    
    %%% update mixture parameters!
%     opts = optimset('display','iter');
%     LB = zeros(1,4);
%     UB = zeros(1,4)+100;
%     lfun = @(x)nll_mix_gamma(x, dv, gams, delta_tdk, mask);
%     prs0 = reshape(mm.mix_gamma(:,:,1:2),1,[]);  % initialize from the last parameter estimation
%     [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
%     mm.mix_gamma(:,:,1:2) = reshape(x, d_states, k_mixtures, 2);
    
    %%% closed-form update!
    N = length(dv);
    pos = (~isinf(log(dv))).*mask;
    gamm_params = mm.mix_gamma*1; 
    for d = 1:d_states
        tau_di = gams(d,:);  % P(state_i | x)
        z_nk = squeeze(delta_tdk(:,d,:));  %P(mixture | state_i) 
        for k = 1:k_mixtures
            sum_z = nansum(tau_di'.*z_nk(:,k));
            sum_zx = nansum((tau_di'.*z_nk(:,k))'.*dv.*pos);
            sum_zxlx = nansum((tau_di'.*z_nk(:,k))'.*dv.*log(dv).*pos);
            sum_zlx = nansum((tau_di'.*z_nk(:,k))'.*log(dv).*pos);
            gamm_params(d,k,1) = sum_z*sum_zx / (sum_z*sum_zxlx - sum_zlx*sum_zx);
            gamm_params(d,k,2) = (sum_z*sum_zxlx - sum_zlx*sum_zx) / sum_z^2;
            
%             sum_z = nansum(z_nk(:,k));
%             sum_zx = nansum((z_nk(:,k))'.*dv.*pos);
%             sum_zxlx = nansum((z_nk(:,k))'.*dv.*log(dv).*pos);
%             sum_zlx = nansum((z_nk(:,k))'.*log(dv).*pos);
%             gamm_params(d,k,1) = sum_z*sum_zx / (sum_z*sum_zxlx - sum_zlx*sum_zx);
%             gamm_params(d,k,2) = (sum_z*sum_zxlx - sum_zlx*sum_zx) / sum_z^2;
            
        end
        
%     gamm_params(d,:,3) = sum(z_nk,1) / N;
        
    end
    mm.mix_gamma = gamm_params;  % place back for iterations
    
end

function [nll] = nll_mVM_mG(x, dth, dcp, dc, dv, gams_t, delta_t_k, Basis, lambda, mask, mm, jj)
    nll_dth = nll_mVM(x(1:14), dth, dcp, dc, dv, gams_t, Basis, lambda, mask);
%     nll_dr = nll_mi_gamma(x(15:end), dv, gams, delta_t_k, mask);
    nll_dr = nll_mi_gamma( dv, gams_t, delta_t_k, mask, mm, jj);
    nll = nll_dth + nll_dr;
end

function [nll] = nll_mi_gamma(dv, gams, delta_t_k, mask, mm, jj)
%     nk = size(delta_t_k,2);
%     dk_params = reshape(x,nk,[]);  % 2 for shape and scale
%     
%     ll = 0;
%     pos = (~isinf(log(dv)));
%     for kk = 1:nk
%         logp = (dk_params(kk,1)-1)*log(dv).*pos - dv/dk_params(kk,2) - dk_params(kk,1)*log(dk_params(kk,2)) - gammaln(dk_params(kk,1));
%         ll = ll + nansum(gams.*delta_t_k(:,kk)'.*logp  .*mask);
%     end
%     nll = -ll;  % return negative loglikelihood
    
    %%%%%%%
    %%% TRY the analytic method here!!!
    %%%%%%%
    %%% closed-form update!
    N = length(dv);
    pos = (~isinf(log(dv))).*mask;
    nk = size(delta_t_k,2);
    gamm_params = mm.mix_gamma(jj,:,:);  %reshape(x, nk, []) ;
    tau_di = gams*1;  % P(state_i | x)
    z_nk = delta_t_k*1; %squeeze(delta_tdk(:,d,:));  %P(mixture | state_i) 
    for k = 1:nk
        sum_z = nansum(tau_di'.*z_nk(:,k));
        sum_zx = nansum((tau_di'.*z_nk(:,k))'.*dv.*pos);
        sum_zxlx = nansum((tau_di'.*z_nk(:,k))'.*dv.*log(dv).*pos);
        sum_zlx = nansum((tau_di'.*z_nk(:,k))'.*log(dv).*pos);
        gamm_params(1,k,1) = sum_z*sum_zx / (sum_z*sum_zxlx - sum_zlx*sum_zx);
        gamm_params(1,k,2) = (sum_z*sum_zxlx - sum_zlx*sum_zx) / sum_z^2;
    end
        
%     gamm_params(1,:,3) = sum(z_nk,1) / N;

    mm.mix_gamma(jj,:,:) = gamm_params;
    pos = (~isinf(log(z_nk)));
    nll = -nansum(nansum(log(z_nk).*pos,1));
    
end

function [nll] = nll_mVM(THETA, dth, dcp, dc, dv, gams, Basis, lambda, mask)
    
    %%% Assume we parameterize in such way first
    alpha_h = THETA(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
    alpha_dc = THETA(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
    alpha_dcp = THETA(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
    kappa_turn = THETA(13)^0.5;     % vairance of the sharp turn von Mises
    kappa_wv = THETA(14)^0.5;       % variance of von Mises
%     gam_shape1 = THETA(15);  % shape and scale for gamma
%     gam_scale1 = THETA(16);
%     gam_shape2 = THETA(17);
%     gam_scale2 = THETA(18);
%     gam_weight = THETA(19);

    %%% construct as kernel
    K_h = (alpha_h * Basis');  % dth kernel
    K_dc = (alpha_dc * Basis');  % dC kernel
    K_dcp = (alpha_dcp * Basis');  % dCp kernel
    
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
    nll_dth = -nansum( mask(2:end) .* gams(2:end) .* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));
%     pos = (~isinf(log(dv)));
%     ll_gam = gam_weight*(gam_shape1-1)*nansum(gams(2:end).*log(dv(2:end)).*pos(2:end)) - nansum(gams(2:end).*dv(2:end))/gam_scale1 - nansum(gams(2:end).*gam_shape1*log(gam_scale1) + gams(2:end).*gammaln(gam_shape1)) ...
%             +(1-gam_weight)*(gam_shape2-1)*nansum(gams(2:end).*log(dv(2:end)).*pos(2:end)) - nansum(gams(2:end).*dv(2:end))/gam_scale2 - nansum(gams(2:end).*gam_shape2*log(gam_scale2) + gams(2:end).*gammaln(gam_shape2));
    nll = nll_dth;% - ll_gam;
end

function [nll] = nll_mix_gamma(x, dv, gams, delta_tdk, mask)
    nd = size(delta_tdk,2);
    nk = size(delta_tdk,3);
    dk_params = reshape(x,nd,nk,2);  % 2 for shape and scale
    
    ll = 0;
    pos = (~isinf(log(dv)));
    for dd = 1:nd
        for kk = 1:nk
            logp = (dk_params(dd,kk,1)-1)*log(dv).*pos - dv/dk_params(dd,kk,2) - dk_params(dd,kk,1)*log(dk_params(dd,kk,2)) - gammaln(dk_params(dd,kk,1));
            ll = ll + nansum(gams(dd,:).*delta_tdk(:,dd,kk)'.*logp  .*mask);
        end
    end
    nll = -ll;  % return negative loglikelihood
end

function [NLL] = analytical_mix_gam(THETA, dis, mask)

    n_mixtures = 3;  % Number of mixture components
% alpha_est = [1 7 ];  % Initial shape parameters
% beta_est = [.5 .2 ];  % Initial rate parameters
alpha_est = [1. 5 6];
beta_est = [.1 .1 .1];
lamb_est = ones(1,n_mixtures)/n_mixtures;  % Initial component probabilities
z_nk = zeros(N, n_mixtures);  % latent z
% EM algorithm
maxIter = 100;  % Maximum number of iterations
log_likelihood = zeros(maxIter, 1);
for iter = 1:maxIter
    % Expectation step (E-step)
    for k = 1:n_mixtures
        z_nk(:, k) = lamb_est(k) * gampdf(x, alpha_est(k), beta_est(k));
    end
    z_nk = z_nk ./ sum(z_nk, 2);
    
    % Maximization step (M-step)
    lamb_est = sum(z_nk,1) / N;
    
        %%% closed-form method(?)
    for k = 1:n_mixtures
        sum_z = sum(z_nk(:,k));
        sum_zx = sum(z_nk(:,k)'.*x);
        sum_zxlx = sum(z_nk(:,k)'.*x.*log(x));
        sum_zlx = sum(z_nk(:,k)'.*log(x));
        alpha_est(k) = sum_z*sum_zx / (sum_z*sum_zxlx - sum_zlx*sum_zx);
        beta_est(k) = (sum_z*sum_zxlx - sum_zlx*sum_zx) / sum_z^2;
    end
    
%     alpha_est = [param_est(1), param_est(3)];
%     beta_est = [param_est(2), param_est(4)];
    
    % Calculate the log-likelihood
    log_likelihood(iter) = nansum(nansum(log(z_nk),1));  %sum(log(sum(z_nk, 2))); %
    
    % Check convergence
    if iter > 1 && abs(log_likelihood(iter) - log_likelihood(iter-1)) < 1e-7
        break;
    end
    disp(['iteration',num2str(iter)])
end
end
