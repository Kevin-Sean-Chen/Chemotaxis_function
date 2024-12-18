% dPAW_gof_strains
%%% this script is used for revision resopnse to the goodness-of-fit across
%%% strains for dPAW.
%%% the idea is to scan through fitted parameters and measure
%%% log-likelihood, then showcase the good and bad fits...

%%
% load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/mutant_gof_vars.mat')
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
n_params = 13+1+0;
L = length(fileList); 
all_ids = cell(1,L+3); % all mutants and three N2s, for ids

rng(13)

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

%% build GOF matrix
GOF = zeros(3, n_strains + 1);
temp = squeeze(reshape(testLL,[3, n_strains, 1]));
GOF(:,2:end) = temp*5/14*log(2);
GOF = GOF([1,3,2],:);

%% for N2 fitting!
rng(37) % 42; 37!

rep = 3;
mle_params = zeros(3, n_params);
trainLL_N2 = zeros(3,1);
testLL_N2 = zeros(3,1);
for cc = 1:3  % condition
    %%% load data
    load(data_n2{cc});
    temp = randperm(length(Data));
    data_sub = Data(temp(1:150));  %150 sub-sample to compare to mutants

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
    testLL_N2(cc) = dPAW_gof(data_sub(ids==2), x_MLE);
    all_ids{cc} = ids;
end

%% loading past CV
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/data4plots/classify_learn3_vars7.mat');
% testLL_N2 = temp.testLL;
% temp_m = squeeze(mean(testLL_N2(:,:,1),1));

%% put in GOF array
% GOF(:,1) = temp_m; %
GOF(:,1) = testLL_N2([1,3,2])*5/14*log(2);

%% directly compute for N2
% for cc = 1:3  % condition
%     %%% load data
%     load(data_n2{cc});
%     %%% fit params
% %     GOF(cc,1) = dPAW_gof(Data(selected(cc,:)), squeeze(mle_params(cc,:)));
%     GOF(cc,1) = dPAW_gof(Data(:), squeeze(mle_params(7,cc,:)));
% end

%% plot results
xtickLabels = {'N2','AIA-','AIB-','AIY-','AIZ-','RIA-'};
figure
bar(GOF')
xticklabels(xtickLabels); set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('gof (bits/s)')

%% check the bad fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading example
% remember to match the test Data and the parameters!!!
%%% for N2
% n2_id = 3;
% load(data_n2{n2_id});
% [xxf, yyf, alltrials, time] = data2xy(Data(all_ids{n2_id}==2));
% x = squeeze(mle_params(n2_id,:));
%%% for mutant
mut_id = 7;
fileName = fileList(mut_id).name; filePath = fullfile(data_path, fileName); load(filePath);
[xxf, yyf, alltrials, time] = data2xy(Data(all_ids{3+mut_id}==1));
x = squeeze(mle_mut_params(mut_id,:));
fileName

%% model fits
%%% unpack parameters
K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7)+0.; Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); K2_ = x(12);  gamma=0.2 %gamma = x(15); %
base_dc = x(13); base_dcp = x(14);

%%% process data
dcp_fit = xxf(2,:);
ang_fit = yyf(1,:);
ddc_fit = xxf(1,:);

%%% reconstructions
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
xx = 0:length(cosBasis)-1;
K_dcp_rec = Amp*exp(-xx/tau);
filt_dcp = conv_kernel(dcp_fit, K_dcp_rec);
K_dc_rec = B_*cosBasis';
filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
xx_h = 1:length(xx)*1;
K_h_rec = Amp_h*exp(-xx_h/tau_h);
filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
dc_dth = filt_ddc + 1*filt_dth;
Pturns = (A_-C_) ./ (1 + exp( -(dc_dth + base_dc))) + C_; 
Pturn_fac = sum(Pturns)/length(Pturns);

%%% plotting
figure
subplot(121)
nbins = 100;
[counts, edges] = histcounts((ang_fit - filt_dcp - base_dcp), nbins);
nCounts = (counts)/sum(counts);
b = bar((edges(2:end)+edges(1:end-1))/2*pi/180, nCounts, 'FaceAlpha',0.8); hold on
b.FaceColor = [0.5 0.5 0.5];   b.EdgeColor = 'none';
bb = (edges(2:end)+edges(1:end-1))/2*pi/180;  %edges*pi/180;
wv_dense = 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb ))*(1-Pturn_fac);
pr_dense = (1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi))*Pturn_fac;
scal_fac = sum(wv_dense+pr_dense);;
plot(bb, ( wv_dense * 1/scal_fac), 'g','LineWidth',2); hold on
plot(bb, ( pr_dense * 1/scal_fac), 'r--', 'LineWidth',2)
set(gcf,'color','w'); set(gca,'Fontsize',20); xlabel('d\theta'); ylabel('pdf')
xticks([-pi, pi]); xticklabels({'-\pi', '\pi'});
subplot(122)
b = bar((edges(2:end)+edges(1:end-1))/2*pi/180, nCounts, 'FaceAlpha',0.8); hold on
b.FaceColor = [0.5 0.5 0.5];   b.EdgeColor = 'none';
plot(bb, ( wv_dense * 1/scal_fac), 'g','LineWidth',2); hold on
plot(bb, ( pr_dense * 1/scal_fac), 'r--', 'LineWidth',2)
ylim([0, 0.005]); set(gcf,'color','w'); set(gca,'Fontsize',20);xlabel('d\theta');
xticks([-pi, pi]); xticklabels({'-\pi', '\pi'});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function
function [LL_fit] = dPAW_gof(Data, x)
    [xx_train, yy_train, mask_train] = data2xy(Data);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
%     fval = -nLL_kernel_hist5(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    fval = (-nLL_kernel_hist5(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit) - ...
            -nLL_randomwalk(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit) );
    LL_fit = fval / sum(trials_fit) / (5/14) / log(2);  % LL/s

    %%% testing...
%     LL_fit = (-pop_nLL(x, Data) - -pop_nLL_randwalk(x, Data)) / sum(mask_train) / (5/14) / log(2);
end

function [x_MLE, fval] = MLE_mGLM(Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_kernel_hist5(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    opts = optimset('display','iter');
    LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.01 -inf, 1e-0, -inf, 1e-1, 1  -inf, -inf];% 0.05];
    UB = [200, 1., ones(1,nB)*inf, 0.5, inf, 50, inf, 100, 100  inf, inf];% 1];
    prs0 = [50, 0.2, randn(1,nB)*10, 0.1, -1, 2, 1, 25, 10 -10, -1 ];%0.2];

    prs0 = prs0 + prs0.*randn(1,length(UB))*0.01;
    %%% tesing in-equality constraints!
    Ai = zeros(1,length(LB));
    Ai([2,7]) = [-1,1];
    bi = 0;  %Ax<=b
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,Ai,bi,[],[],LB,UB,[],opts);
    x_MLE = x;
end
