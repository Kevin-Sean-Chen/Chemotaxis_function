% staPAW_gamma
%%% serve as ground-truth to varify the inference procedure
%%% this version contains speed distribution with gamma, together with the
%%% angular mixture model; in addition, the state transitions are driven by
%%% concentration input!
%%% different from mGLMHMM setups, the gamma speed is under the beta
%%% decision in the original dPAW model, so the pirouette decision
%%% corresponds to a speed profile too (if this does not work, we would
%%% consider mixture-gamma or gamma-GLM frameworks)

%% generate hidden states, as priors
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
alpha = 1.5;  % Direchelt prior
kappa = 0.;  % stickiness
Pij = A0*1;  % as baseline transition probability

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

%%% transition kernels
alpha_ij = randn(nStates,nStates, nB)*1;  % state x state x weights
lk = length(K_h1);
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

%% run distributions
alphas = [1, 10;
          4,  5];  % shape
betas = [5, 2;
         3, 8];  % scale
phis = [.1, .5;
        .8, .2]*.5;  % autoregressive
    
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
    %%% compute states (currently only driven by dC)
    Pu_ij = zeros(nStates,nStates);
    for ii = 1:nStates
        for jj = 1:nStates
            if ii ~= jj
                Pu_ij(ii,jj) = exp(dC(tt-pad+1:tt)*squeeze(K_ij(ii,jj,:))) + (Pij(ii,jj));  % input-driven state transitions
            else
                Pu_ij(ii,jj) = 0 + Pij(ii,jj);
            end
        end
    end
    Pu_ij = Pu_ij ./ (sum(Pu_ij,2)+1e-8);  % normalize row, conditional probability
    
    temp_state = find(rand < cumsum(Pu_ij(X(tt-1),:)));  % Markov transition
    X(tt) = temp_state(1);
    
    if X(tt) == 1
        K_h = K_h1; K_dc = K_dc1; K_dcp = K_dcp1; kappa_turn = kappa_turn1; kappa_wv = kappa_wv1; aa = alphas(1,:); bb = betas(1,:); phi = phis(1,:);
    elseif X(tt) == 2
        K_h = K_h2; K_dc = K_dc2; K_dcp = K_dcp2; kappa_turn = kappa_turn2; kappa_wv = kappa_wv2; aa = alphas(2,:); bb = betas(2,:); phi = phis(2,:);
    end
    %%% compute emissions
    F(tt) = dC(tt-pad+1:tt)*K_dc' + abs(dth(tt-pad:tt-1))*r2d*K_h';  % linear filtering
    turns(tt) = choice(NL(F(tt)+base,A));  % nonlinearity and binary choice
    dth(tt) = turns(tt)*circ_vmrnd(pi,kappa_turn,1) + (1-turns(tt))*circ_vmrnd(dCp(tt-pad+1:tt)*K_dcp',kappa_wv,1);  % angle drawn from mixture of von Mesis
    dv(tt) = turns(tt)*gamrnd(aa(1) + phi(1)*dv(tt-1), bb(1), 1, 1) + (1-turns(tt))*gamrnd(aa(2) + phi(2)*dv(tt-1), bb(2), 1, 1);
    %gamrnd(aa, bb, 1, 1);  % gamma random variable
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
nX = 14+2+2 +2;  % number of input dimensions (i.e., dimensions of regressor)
nY = 1;  % number of output dimensions 
nT = length(yy); % number of time bins
loglifun = @logli_GTmVMG;  % ground truth model log-likelihood function
loglitrans = @logli_trans;  % log-likelihood of soft-max state transitions

% basis function
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);

% Set linear weights & output noise variances
% wts0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.]; 
wts0 = rand(nY,nX,nStates); % parameters for the mixture-VonMesis behavioral model
wts0(1,:,1) = [alpha_h1, alpha_dc1, alpha_dcp1, kappa_turn1, kappa_wv1, alphas(1,:), betas(1,:), phis(1,:)]; %single mGLM
wts0(1,:,2) = [alpha_h2, alpha_dc2, alpha_dcp2, kappa_turn2, kappa_wv2, alphas(2,:), betas(2,:), phis(2,:)];
wts0 = wts0 + wts0.*rand(nY,nX,nStates)*0.1; % small perturbation

wts_state0 = alpha_ij;

% Build struct for initial params
mmhat = struct('A',A0,'wts',wts0,'wts_state',wts_state0,'loglifun',loglifun,'loglitrans',loglitrans,'basis',cosBasis,'lambda',.0, 'Mstepfun',@runMstep_GTmVMG);

%%
% [logp,gams,xisum] = runFB_GLMHMM(mmhat,xx,yy,mask)
% test = runMstep_GTmVM(mmhat, xx, yy, gams, mask)
% test2 = logli_GTmVM(mmhat, xx, yy, mask)
% -nLL_gamma([10, 2], dv, mask)
% [logli] = logli_trans(mmhat,xx,yy,mask);
% [logli] = logli_GTmVM(mmhat, xx, yy, mask);
% [logp,gams, xis, xisum] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
% gams

%% Set up variables for EM
maxiter = 50;
EMdisplay = 2;
logpTrace = zeros(maxiter,1); % trace of log-likelihood
dlogp = inf; % change in log-likelihood
logpPrev = -inf; % prev value of log-likelihood
jj = 1; % counter

while (jj <= maxiter) && (dlogp>1e-3)
    
    % --- run E step  -------
    [logp,gams, xis, xisum] = runFB_GLMHMM_xi(mmhat,xx,yy,mask); % run forward-backward, returning the full xi series
    logpTrace(jj) = logp;
   
    % --- run M step  -------
    
    % Update transition matrix
%     mmhat.A = (alpha-1 + xisum) ./ (nStates*(alpha-1) + sum(xisum,2)); % normalize each row to sum to 1
    
    % Update transition matrix (with stickiness)
%     normed_xisum = xisum ./ sum(xisum,2);
%     unormedA = kappa*eye(nStates) + (alpha-1)*ones(nStates,nStates) + normed_xisum;
%     mmhat.A = unormedA ./ sum(unormedA,2);
    
    % Update model params
    mmhat = runMstep_GTmVMG(mmhat, xx, yy, gams, mask);
    % Update for input-driven transitions
    mmhat = runMstep_trans(mmhat, xx, yy, xis, mask);  % Figure out if its alphas_ or xis that should be used!!!!!!
    % logp = logp0;
    % logp = logp + w'*u;
    % logp = logp - logsumexp(logp);
    % obj = sum(xis.*logp);l  %(xis is xisum before summation)
    % 
    % another function for optimization of driven transition
    
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


%% compare to gound-truth kernels
stateK = 1;
figure()
for kk = 1:stateK
x = squeeze(mmhat.wts(:,:,kk));
alpha_h_ = x(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
alpha_dc_ = x(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
alpha_dcp_ = x(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)

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
nbins = 50;
hh = histogram(dth(pos), nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1);
% [aa,bb] = hist((dth - filt_dcp)*1 , 1000);
scal = sum(1/(2*pi*besseli(0,kappa_wv_)) * exp(kappa_wv_*cos( bb )) * p_wv + (1/(2*pi*besseli(0,kappa_turn_)) * exp(kappa_turn_*cos( bb-pi ))) * p_brw);
plot( bb, 1/(2*pi*besseli(0,kappa_wv_)) * exp(kappa_wv_*cos( bb )) * p_wv / scal , 'b'); hold on
plot( bb, ((1/(2*pi*besseli(0,kappa_turn_)) * exp(kappa_turn_*cos( bb-pi ))) * p_brw)/ scal ,'r'); xlim([-pi,pi])
title('von Mises for \delta C^{\perp}')

%% density for speed
% kk = 1;
pos = find(X(1:end-1)==kk);
x = squeeze(mmhat.wts(:,:,kk));
gam_shapes_ = x(15:16);      % shape parameter for gamma
gam_scales_ = x(17:18);      % scale parameter for gamma
gam_phis_ = x(19:20);        % autoregressive for gamma
figure;
dv_k = dv(pos);
p_turn_k = Pturns(pos);
hh = histogram(dv_k, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end)+.0;

binIndices = discretize(dv_k, bb);
prob_dr_wv = zeros(1,length(bb));
prob_dr_pr = zeros(1,length(bb));
for bin = 1:numel(bb)-1
    pointsInBin = find(binIndices == (bin));
    pointsInBin = pointsInBin(~isnan(pointsInBin));
    prob_dr_pr(bin) = mean(gampdf(zeros(1,length(pointsInBin)-1)+bb(bin), gam_shapes_(1)+gam_phis_(1)*dv_k(pointsInBin(1:end-1)), gam_scales_(1)).*p_turn_k(pointsInBin(2:end)));
    prob_dr_wv(bin) = mean(gampdf(zeros(1,length(pointsInBin)-1)+bb(bin), gam_shapes_(2)+gam_phis_(2)*dv_k(pointsInBin(1:end-1)), gam_scales_(2)).*(1-p_turn_k(pointsInBin(2:end))));
end
plot(bb, prob_dr_pr)
plot(bb, prob_dr_wv)
% gam_pr = gampdf(bb, gam_shapes_(1), gam_scales_(1));
% gam_wv = gampdf(bb, gam_shapes_(2), gam_scales_(2));
% scal = sum(gam_pr*p_brw + gam_wv*p_wv);
% plot(bb, p_brw/scal*gam_pr)
% plot(bb, p_wv/scal*gam_wv)
% plot(bb, p_wv/scal*gam_wv + p_brw/scal*gam_pr)

% gam_pr = gamrnd(gam_shapes_(1)+gam_phis_(1)*dv(1:end-1), gam_scales_(1)).*Pturns(1:end-1);
% gam_wv = gamrnd(gam_shapes_(2)+gam_phis_(2)*dv(1:end-1), gam_scales_(2)).*(1-Pturns(1:end-1));
% [counts_pr, edges_pr] = histcounts(gam_pr, bb);
% [counts_wv, edges_wv] = histcounts(gam_wv, bb);
% scal = sum(counts_pr*p_brw + counts_wv*p_wv);
% plot(bb(2:end), 1/scal*counts_pr)
% plot(bb(2:end), 1/scal*counts_wv)
% plot(bb(2:end), 1/scal*(counts_wv+counts_pr))

%% transitional kernels
figure
alpha_tran_ = squeeze(mmhat.wts_state(1,2,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(fliplr(K_trans_)); hold on
alpha_tran_ = squeeze(mmhat.wts_state(2,1,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(fliplr(K_trans_)); hold on
plot(squeeze(K_ij(1,2,:)),'--')
plot(squeeze(K_ij(2,1,:)),'--'); hold off

%% debugging
% x = [reshape(mmhat.wts_state,1,[])    reshape(mmhat.A,1,[])];
% nll_MLR(x, dc, alphas_, Basis, mask)

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
function [logli] = logli_trans(mm,xx,yy,mask)
    dc = xx(1,:);  % input driving the state transitions
    nStates = size(mm.A,2);
    T = length(xx);
    logli = zeros(nStates, nStates, T);
    alpha_ij = mm.wts_state;
    Basis = mm.basis;
    Kij = zeros(nStates, nStates, size(Basis,1));  % transition filters
    Pij = mm.A;
%     Pij = Pij./sum(Pij, 2);  % transition matrix
    for ii = 1:nStates
        for jj = 1:nStates
            if ii ~=jj
                K_ij_dc = (squeeze(alpha_ij(ii,jj,:))'*Basis');  % reconstruct kernel
                Kij(ii,jj,:) = K_ij_dc;
                logli(ii,jj,:) = log(exp(conv_kernel([dc(2:end) 0], K_ij_dc)) + (Pij(ii,jj)));  % taking the log of numerater, apparently the time shift is important!!
            else
                logli(ii,jj,:) = 0 + log(0+Pij(ii,jj));  % remove diagonal components, important to remove zero-kernels!
            end
        end
        logli(ii,:,:) = logli(ii,:,:) - logsumexp(logli(ii,:,:));
    end
end

function mm = runMstep_trans(mm, xx, yy, xis, mask)
    dc = xx(1,:);  % input driving the state transitions
    nStates = size(mm.wts_state,1);
    Basis = mm.basis;
    Pij_ = mm.A;  % baseline transitions
    lfun = @(x)nll_MLR(x, dc, xis, Basis, mask);  % multi-nomial logistic regression loss function
    prs0 = [reshape(mm.wts_state,1,[])    reshape(Pij_',1,[])];  % concatenate flatten parameters
    LB = [zeros(1,length(prs0))-100    zeros(1,nStates^2)-100];
    UB = [zeros(1,length(prs0))+100    zeros(1,nStates^2)+100];
    Aeq = zeros(nStates, length(prs0));   % transition probability sums to one!
    lw = length(reshape(mm.wts_state,1,[]));
    for ii = 1:nStates
        Aeq(ii,lw+(ii-1)*nStates+1:lw+ii*nStates) = ones(1,nStates);  % "smart" parameterizing with the last nStates^2 being state transitions
    end
    beq = ones(nStates,1);  % this all makes the constrain: Ax = b
%     [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],Aeq,beq,LB,UB,[]);  % constraint needed only    if A is the probability matrix itself!
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[]);
    mm.wts_state = reshape(x(1:end-nStates^2), nStates, nStates, []);   % update weights
    mm.A = reshape(x(end-nStates^2+1:end), nStates, nStates)';  % update baseline, note the transpose here, at flatten, and corresponding to the constaints, important!!
end

function nll = nll_MLR(x, dc, xis, Basis, mask)
    
%     logP = log(alphas_);  % target transitional probability
    logP_ = xis*0;  % reconstructed transitional probability (with driven and baseline)
    nStates = size(xis,1);
    Kij = zeros(nStates, nStates, size(Basis,1));
    alpha_ij = reshape(x(1:end-nStates^2), nStates, nStates,[]);  % weights on basis for kernels
    Pij = reshape(x(end-nStates^2+1:end), nStates, nStates)';  % basline probability, important to transpose given the way its flattened!!
%     Pij = Pij./sum(Pij, 2);  % transition matrix
    for ii = 1:nStates
        for jj = 1:nStates
            if ii ~= jj
                K_ij_dc = (squeeze(alpha_ij(ii,jj,:))'*Basis');  % reconstruct kernel
                Kij(ii,jj,:) = K_ij_dc;
                logP_(ii,jj,:) = log(exp(conv_kernel([dc(2:end) 0], K_ij_dc)) + (Pij(ii,jj)));
            else
                logP_(ii,jj,:) = 0 + log(0+Pij(ii,jj));
            end
        end
        logP_(ii,:,:) = logP_(ii,:,:) - logsumexp(logP_(ii,:,:));
    end
    
    nll = -squeeze(sum(sum(xis(:,:,1:end) .* logP_(:,:,1:end))))'*mask(1:end)';  %%% currently debugging if this is mismatched!!! %%%
end

function [logli] = logli_GTmVMG(mm, xx, yy, mask)

% Compute log-likelihood term under a mixture of von Mesis and gamma
% distribution, with input-dependent mixing weights
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%      .wts   [1 len(param) K] - per-state parameters for the model
%      .basis [len(alpha) len(kernel)] - basis functions used for the kernel
%      .lambda -scalar for regularization of the logistic fit
%    xx [2 T] - inputs (time series of dc,dcp, and dth)
%    yy [2 T] - outputs (dth angles and dis speed)
%
% Output
% ------
%  logpy [T K] - loglikelihood for each observation for each state
    
    %%% loading parameters
    THETA = mm.wts;
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
        gamm_shape_turn = THETA(1,15,k);     % shape of gamma for turn
        gamm_shape_wv = THETA(1,16,k);     % shape of gamma for wv
        gamm_scale_turn = THETA(1,17,k);     % scale of gamma for turn
        gamm_scale_wv = THETA(1,18,k);     % scale of gamma for wv
        gamm_AR_turn = THETA(1,19,k);      % weight for turn
        gamm_AR_wv = THETA(1,20,k);        % weight for wv

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
        VM_wv = C * exp(kappa_wv^2*cos(( filt_dcp - dth(2:end) )*d2r));  %von Mises distribution
%         gamm_wv = 1/(gamma(gamm_shape_wv)*gamm_scale_wv^gamm_shape_wv) .*dv.^(gamm_shape_wv-1) .* exp(-dv./gamm_scale_wv);
        gamm_wv = gampdf(dv(2:end), gamm_shape_wv + gamm_AR_wv*dv(1:end-1), gamm_scale_wv);
        
        %%% turning analge model
        VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
        gamm_turn = gampdf(dv(2:end), gamm_shape_turn + gamm_AR_turn*dv(1:end-1), gamm_scale_turn);
        
        %%% marginal probability
        marginalP = (1-P).*VM_wv.*gamm_wv + P.*VM_turn.*gamm_turn;
%         pos = (~isinf(log(dv)));
%         ll_gamm = (gamm_shape-1)*log(dv.*pos) - dv/gamm_scale - gamm_shape*log(gamm_scale) - gammaln(gamm_shape);
        logli(2:end,k) = ( mask(2:end).* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));  
        %+ mask(2:end).*ll_gamm(2:end);
        %
    end
end

function mm = runMstep_GTmVMG(mm, xx, yy, gams, mask)
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
LB = [ones(1,12)*-10, 0, 0.,0, 0.,  0.,  0. 0 0];
UB = [ones(1,12)*10, 20, 20.,20, 20., 100, 100 1 1];
% prs0 = rand(1,10);
% prs0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.] ;
% prs0 = [9.9763  -0.5343  -0.0776   0.1238  -0.0529   0.5335   7.7254  367.3817  0.1990  1.0000  0.1000]; %single mGLM
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

for jj = 1:nStates
    lfun = @(x)nll_mVMG(x, dth, dcp, dc, dv, gamnrm(jj,:), Basis, lambda, mask);
%     prs0 = prs0 + prs0.*randn(1,length(UB))*0.5;
    prs0 = mm.wts(:,:,jj); %+ mm.wts(:,:,jj).*(2*(rand(1,length(UB))-0.5))*0.5;  %from last time!
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
%     x = fminunc(lfun, prs0);
    mm.wts(:,:,jj) = x; % weighted optimization
end


end

function [nll] = nll_mVMG(THETA, dth, dcp, dc, dv, gams, Basis, lambda, mask)
    
    %%% Assume we parameterize in such way first
    alpha_h = THETA(1:4);       % kernel for dth angle history kernel (weights on kerenl basis)
    alpha_dc = THETA(5:8);      % kernel for dC transitional concentration difference (weights on kerenl basis)
    alpha_dcp = THETA(9:12);    % kernel for dCp perpendicular concentration difference (weights on kerenl basis)
    kappa_turn = THETA(13)^0.5;     % vairance of the sharp turn von Mises
    kappa_wv = THETA(14)^0.5;       % variance of von Mises
    ks_gamm = THETA(15:16);  % shapes for gammas
    thetas_gamm = THETA(17:18);  % scales for gammas
    phi_gamm = THETA(19:20);  % AR-weights for gamma
    
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
    gamm_wv = gampdf(dv(2:end), ks_gamm(2) + phi_gamm(2)*dv(1:end-1), thetas_gamm(2));
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
    gamm_turn = gampdf(dv(2:end), ks_gamm(1) + phi_gamm(1)*dv(1:end-1), thetas_gamm(1));
    
    %%% marginal probability
    marginalP = (1-P).*VM.*gamm_wv + P.*VM_turn.*gamm_turn;
    nll_dth = -nansum( mask(2:end) .* gams(2:end) .* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));
%     pos = (~isinf(log(dv)));
%     ll_gam = (theta_gamm(1)-1)*nansum(gams(2:end).*log(dv(2:end)).*pos(2:end)) - nansum(gams(2:end).*dv(2:end))/theta_gamm(2) - nansum(gams(2:end).*theta_gamm(1)*log(theta_gamm(2)) + gams(2:end).*gammaln(theta_gamm(1)));
    nll = nll_dth;% - ll_gam;
end

function [NLL] = nLL_gamma(THETA, dis, mask)

    %%% loading parameters
    k1 = THETA(1);
    theta1 = THETA(2);
    
    %%% log-likelihood
    N = length(dis);
    dis = mask(1:end).*dis;  % remove masked elements
    pos = (~isinf(log(dis)));
    L1 = (k1-1)*nansum(log(dis).*pos) - nansum(dis)/theta1 - N*(k1*log(theta1) + gammaln(k1));
%     L2 = (k2-1)*nansum(alpha(:,2)'.*log(dis)) - nansum(alpha(:,2)'.*dis)/theta2 - 1*(N*k2*log(theta2) + N*log(gamma(k2)));
    %%% log_likelihood = sum((k-1)*log(data) - data/theta - k*log(theta) -gammaln(k));  %nansum(alpha(:,2))
    
    marginalP = L1;% + L2;
    NLL = -marginalP;
end
