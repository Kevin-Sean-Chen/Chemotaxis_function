% dPAW_gof_strains
%%% this script is used for revision resopnse to the goodness-of-fit across
%%% strains for dPAW.
%%% the idea is to scan through fitted parameters and measure
%%% log-likelihood, then showcase the good and bad fits...

%% load data
%%% mutant fits
% temp = load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/mle_mut_params7.mat');
temp = load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/mutant_fit_vars.mat');
mle_mut_params = temp.mle_mut_params;

%%% mutant Data
data_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/Data_files';
fileList = dir(fullfile(data_path, '*.mat'));

%%% N2 Data
data_n2 = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat'};

%%% N2 fits
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/data4plots/Kfold_mle_param7.mat');
temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/data4plots/classify_learn3_vars7.mat');
mle_params = temp.mle_params;

%% setup numbers
n_strains = 5;
n_params = 13+0;
L = length(fileList); 
all_ids = cell(1,L+3); % all mutants and three N2s, for ids

rng(37)

%% looping -- rerun mutant fits
rep = 3;
mle_mut_params = zeros(L, n_params); % L x N
trainLL = zeros(L,1);
testLL = zeros(L,1);
%%
for li = 1:numel(fileList)
    
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    
    %%% training, with small repetition
    x_temp = zeros(rep,n_params);
    fv_temp = zeros(1,rep);
    data_id = 1:length(Data);
    ids = crossvalind('Kfold',data_id, 2);
    for ri = 1:rep   % test with repreats to find the max MLE
        [x_MLE, fval] = MLE_mGLM(Data(ids==1));  %(randperm(length(Data))));
        fv_temp(ri) = fval;
        x_temp(ri,:) = x_MLE;
    end
        
    % record MLE
    [minv,pos] = min(fv_temp);
    x_MLE = x_temp(pos,:);
    mle_mut_params(li, :) = x_MLE;
    trainLL(li) = minv;
    
    % test LL
    testLL(li) = dPAW_gof(Data(ids==2), x_MLE);
    all_ids{li+3} = ids;  % record cv index
end

%%
GOF = zeros(3, n_strains + 1);
temp = squeeze(reshape(testLL,[3, n_strains, 1]));
GOF(:,2:end) = temp*5/14*log(2);
GOF = GOF([1,3,2],:);

%% directly compute for mutants
% mut_param = reshape(mle_mut_params,[3,n_strains,n_params]);
% L = length(fileList);  % K-fold cross-validation
% GOF = zeros(3, n_strains + 1);  % +1 for N2 as the first
% ii = 1;
% for ss = 2:n_strains+1  % strain
%     for cc = 1:3  % condition
%         %%% load data
%         fileName = fileList(ii).name;
%         filePath = fullfile(data_path, fileName);
%         load(filePath);
%         ii = ii+1;
%         
%         %%% fit params
%         GOF(cc,ss) = dPAW_gof(Data(floor(length(Data)/1):end), squeeze(mut_param(cc,ss-1,:)));
%     end
% end

%% for N2 fitting!
rep = 3;
mle_params = zeros(3, n_params);
trainLL_N2 = zeros(3,1);
testLL_N2 = zeros(3,1);
for cc = 1:3  % condition
    %%% load data
    load(data_n2{cc});
    data_sub = Data(1:150);  % sub-sample to compare to mutants
    %%% training, with small repetition
    x_temp = zeros(rep,n_params);
    fv_temp = zeros(1,rep);
    data_id = 1:length(data_sub);
    ids = crossvalind('Kfold',data_id, 2);
    for ri = 1:rep   % test with repreats to find the max MLE
        [x_MLE, fval] = MLE_mGLM(data_sub(ids==1));  %(randperm(length(Data))));
        fv_temp(ri) = fval;
        x_temp(ri,:) = x_MLE;
    end
        
    % record MLE
    [minv,pos] = min(fv_temp);
    x_MLE = x_temp(pos,:);
    mle_params(cc, :) = x_MLE;
    trainLL_N2(cc) = minv;
    
    % test LL
    testLL_N2(cc) = dPAW_gof(Data(ids==2), x_MLE);
    all_ids{cc} = ids;
end

%% loading past CV
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/data4plots/classify_learn3_vars7.mat');
% testLL_N2 = temp.testLL;
% temp_m = squeeze(mean(testLL_N2(:,:,1),1));

%% put in GOF array
% GOF(:,1) = temp_m; %
GOF(:,1) = testLL_N2([1,3,2]);

%% directly compute for N2
% for cc = 1:3  % condition
%     %%% load data
%     load(data_n2{cc});
%     %%% fit params
% %     GOF(cc,1) = dPAW_gof(Data(selected(cc,:)), squeeze(mle_params(cc,:)));
%     GOF(cc,1) = dPAW_gof(Data(:), squeeze(mle_params(7,cc,:)));
% end

%% plot results
xtickLabels = {'N2','AIA','AIB-','AIY-','AIZ-','RIA-'};
figure
bar(GOF')
xticklabels(xtickLabels); set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('gof (bits/s)')

%% check the bad fits...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(data_n2{1});
% fileName = fileList(7).name; filePath = fullfile(data_path, fileName); load(filePath);
[xxf, yyf, alltrials, time] = data2xy(Data(:));
figure;
hist(xxf(2,:),100)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function
function [LL_fit] = dPAW_gof(Data, x)
    [xx_train, yy_train, mask_train] = data2xy(Data);
%     nB = 4;
%     [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
%     ang_fit = yy_train;
%     dcp_fit = xx_train(2,:);
%     ddc_fit = xx_train(1,:);
%     trials_fit = mask_train;
% %     fval = -nLL_kernel_hist5(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
%     fval = (-nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit) - ...
%             -nLL_randomwalk(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit) );
%     LL_fit = fval / sum(trials_fit) / (5/14) / log(2);  % LL/s

    %%% testing...
    LL_fit = (-pop_nLL(x, Data) - -pop_nLL_randwalk(x, Data)) / sum(mask_train) / (5/14) / log(2);
end

function [x_MLE, fval] = MLE_mGLM(Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
%     lfun = @(x)nLL_kernel_hist5(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    lfun = @(x)nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    opts = optimset('display','iter');
%     LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1  -inf, -inf];
%     UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 100  inf, inf];
%     prs0 = [50, 0.2, randn(1,nB)*10, 0.01, -1, 2, 1, 25, 10 -10, -1];
    LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1  0.2];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 100, 1];
%     prs0 = [50, 0.2, randn(1,nB)*10, 0.01, -1, 2, 1, 25, 10 0.5];
    prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.01;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end
