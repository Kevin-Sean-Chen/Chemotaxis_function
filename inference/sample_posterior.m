% sample_posterior
%%%
% a test script to sample posterior parameter distribution for the mVM
% project, in order to characeterize parameter spread
%%%

% for now, this continues from the Chemotaxis_inference script with those
% _fit data vector ready.
%% likelihood function
logLike = @(x)nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);

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

x0 = [50, 0.5, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.  0 10]; 
ball = randn(length(x0),60)*0.1;
mball=bsxfun(@plus,x0',ball);

%% Apply the hammer:
%
% Draw samples from the posterior. 
%
tic
% m=gwmcmc(mball,{logprior logLike},10000,'burnin',.3,'stepsize',2);
m=gwmcmc(mball,{logprior logLike},20000,'burnin',.5,'stepsize',5);
toc

%% post-analysis
%% Auto-correlation function
figure
[C,lags,ESS]=eacorr(m);
plot(lags,C,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
title('Markov Chain Auto Correlation')

%% Corner plot of parameters
figure
ecornerplot(m,'ks',true,'color',[.6 .35 .3])

