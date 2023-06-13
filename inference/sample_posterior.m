% sample_posterior
%%%
% a test script to sample posterior parameter distribution for the mVM
% project, in order to characeterize parameter spread
%%%

% for now, this continues from the Chemotaxis_inference script with those
% _fit data vector ready.
%% likelihood function
% for concatenated tracks
% logLike = @(x)nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);

% for Data structure
Data_fit = Data(drange(1:100));  %Data(1:100); %
logLike = @(x)pop_nLL(x, Data_fit);

%% Prior information
logprior = @(x)logPrior_mVM(x, cosBasis, 0.1);

%% Find the posterior distribution using GWMCMC
% Now we apply the MCMC hammer to draw samples from the posterior.
% first we initialize the ensemble of walkers in a small gaussian ball 
% around the m0 estimate. 

%%% information for maximum likelihood
% LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0 -inf, 1e-0, -inf, 1e-1, 1e-0*2, 0.1    -inf -180];
% UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1    inf 180];
% prs0 = [5, 0.01, randn(1,nB)*10, 0.01, 10, 25, 10, 25, 5, 1.    0 10]; 
% K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7)/1; Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); sb = x(12); K2_ = x(13);  gamma = x(14);
%%%

x0 = [80, 0.6, randn(1,nB)*10, 0.01, -1, 2, 1, 20, 5, 1.  0 10];
x0 = x_MLE;  % from MLE!!
n_walkers = 70;
% ball = randn(length(x0),n_walkers) .* repmat(x0',1,n_walkers)*0.1;  % scaled noise
ball = randn(length(x0),n_walkers)*0.1;  % uniform noise
mball=bsxfun(@plus,x0',ball);
% mball(2,:) = min(mball(2,:),1);  % specific for parameter bound

%% Apply the hammer:
%
% Draw samples from the posterior. 
%
tic
% m=gwmcmc(mball,{logprior logLike},10000,'burnin',.3,'stepsize',2);
% [m, logps]=gwmcmc(mball,{logprior logLike},10000,'burnin',.3,'stepsize',2);
[m, logps] = gwmcmc(mball,{logprior logLike},10000,'ThinChain',5,'burnin',.2);
toc

%% post-analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the finite runs
check_inf = sum(squeeze(logps(1,:,:)),2);  % check if logp for prior was extremely small (out of bound)
pos = find(check_inf<-10^10);
m_bound = m;
logps_clean = logps;
logps_clean(:,pos,:) = [];
m_bound(:,pos,:) = [];

%% Auto-correlation function
figure
[C,lags,ESS]=eacorr(m_bound);
plot(lags,C,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
title('Markov Chain Auto Correlation')

%% Corner plot of parameters
figure
ecornerplot(m_bound,'ks',true,'color',[.6 .35 .3])

%%% save('/projects/LEIFER/Kevin/Data_learn/N2/..._mcmc.mat','m'

