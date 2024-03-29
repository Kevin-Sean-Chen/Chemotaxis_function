% CV_mVMHMM
% high level script that calls fitting fucnctions, used to scan hyper
% parameters and used for cross-validation given group data.

%% load data
% loading track-based Data structure, with dc, dcp, dth, mask information

%% Model fitting options
% states, prior, and regularization
nStates = 2; % number of latent states
alpha = 1.;  % Dirichlet shape parameter as a prior
kappa = 0.2;  % upweighting self-transition for stickiness
lambs = zeros(nStates);  % constraints on turning amplitude in for each state

% iteration variables
maxiter = 50;  % EM iterations
dlogp_tol = 1e-3;  % tolerance for delta logP

% construct an option structure for cleaness
fit_option.nStates = nStates;
fit_option.alpha = alpha;
fit_option.kappa = kappa;
fit_option.lambs = lambs;
fit_option.maxiter = maxiter;
fit_option.dlogp_tol = dlogp_tol;

%% loop for CV
K = 7;  % K-fold cross-validation
rep = 1;  % repetition might be required
testLL = zeros(rep,K);
data_id = 1:length(Data);  % original track ID
indices = crossvalind('Kfold',data_id, K);  % indices for CV
logTraces = cell(rep,K);  % record all log-traces

for ri = 1:rep
for ki = 1:K
    % define train and test sets
    test_set = (indices==ki);
    train_set = ~test_set;
    
    % load individual tracks and concatenation
    Data_train = Data(train_set);
    Data_test = Data(test_set);
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    [xx_test, yy_test, mask_test] = data2xy(Data_test);
    
    % train the model
    [mmhat, logpTrace] = runEM_mVMHMM(Data_train, xx_train, yy_train, mask_train, fit_option);
    logTraces{ri,ki} = logpTrace;
    
    % run model test on held-out data with fitted model parameters
    %%% for null model %%%
    mmhat.wts(1,:,1) = [mmhat.wts(1,1,1), zeros(1,nB), 0, 1, 0, 1, mmhat.wts(1,end-1,1), mmhat.wts(1,end,1)];
    %%%%%%%
    [logp_test, gams_test, xisum_test] = runFB_GLMHMM(mmhat,xx_test,yy_test,mask_test);
    testLL(ri, ki) = logp_test/length(mask_test);

end
end

%% plotting test LL
figure;
plot(testLL','o')
ylabel('test LL')

%% plotting logP traces
figure;
for ki = 1:K
    temp = logTraces{1,ki};
    temp(temp==0) = NaN;
    plot(temp,'k-','LineWidth',2)
    hold on
end
xlabel('iterations')
ylabel('$\mathcal{L}$','Interpreter','latex')
