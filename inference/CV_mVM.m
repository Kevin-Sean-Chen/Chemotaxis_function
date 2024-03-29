% CV_mVM
% this is used for the mixture of von Mesis model, without latent states
% high level script that calls fitting fucnctions, used to scan hyper
% parameters and used for cross-validation given group data.

%% load data
% loading track-based Data structure, with dc, dcp, dth, mask information
% load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app.mat')
% load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai2.mat')
load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test.mat')

%% Model fitting options
% loading data structure
drange = randperm(length(Data));%[1:length(Data)]; %
Data_fit = Data(drange(:));  %Data(1:100); %

% fitting options, bounds, and initial conditions
opts = optimset('display','iter');
LB = [1e-0, 1e-1, ones(1,nB)*-inf, 0 -inf, 1e-0, -inf, 1e-1, 0., 0.1, -inf, -180];
UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1, inf, 180];
prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1., 0, 10];
prs0 = prs0 + prs0.*randn(1,length(UB))*0.2;
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

%% loop for CV
K = 10;  % K-fold cross-validation
rep = 1;  % repetition might be required
testLL = zeros(rep,K);  % test log-likelihood
data_id = 1:length(Data);  % original track ID
indices = crossvalind('Kfold',data_id, K);  % indices for CV
logVal = zeros(rep,K);  % record fitted log-values
fitted_params = zeros(rep,K, length(prs0));  % record the fitted parameters

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
    lfun = @(x)pop_nLL(x, Data_train);
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    logVal(ri,ki) = fval / length(yy_train);
    fitted_params(ri,ki,:) = x;
    
    % run model test on held-out data with fitted model parameters
    %%% for null model %%%
%     mmhat.wts(1,:,1) = [mmhat.wts(1,1,1), zeros(1,nB), 0, 1, 0, 1, mmhat.wts(1,end-1,1), mmhat.wts(1,end,1)];
    null_params = [x(1:2), zeros(1,nB),x(7),0,1,0,1,x(12),x(13),0,10];
    null_ll = lfun(null_params);
    %%%%%%%
    testLL(ri,ki) = -(pop_nLL(x, Data_test)-null_ll) / length(yy_test) / (5/14);
    %%%% test for different model ablation here~~~
    %%%%%%%%%%%
    %%%%%%%%%%%
    %%%%%%%%%%%

end
end

%% plotting test LL
figure;
plot(testLL','o')
ylabel('test LL')
