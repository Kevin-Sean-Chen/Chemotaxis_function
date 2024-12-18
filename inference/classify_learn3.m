% classifiy_learn
%%% classify learning conditions given parameter samples
%%% with proper K-fold cross-validation in earch loop
%% load tracks from each condition
% load data files
datas = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test2.mat'};
     
%% loop through K-folds, conditions, and scaling
rep = 5;  % repetition per fold to make sure that it is the better fit
K = 7;  % K-fold cross-validation
cond = length(datas);  % all exp conditions
n_params = 13;  % number of parameters in our model for now
mle_params = zeros(K, cond, n_params); % K x c x N
null_params = zeros(K, cond, n_params);  % K x c x N, for random walk null model
bin_params = zeros(K, cond, 1);  % K x c x 1, for binomial distribution (chemotaxis index)
kc_control_params = zeros(K, cond, n_params);  % control anaylsis ablating kernels
kcp_control_params = zeros(K, cond, n_params);
Mdls = cell(1,K);  % models for naive bayes for dC clustering

%% %%%%%%%  TRAINING  %%%%%%% 
%%% seperate indices a head of time
ids = cell(1,3);
data_ids = zeros(3, 400);  %place holder for scambled ID for each condition
for ci = 1:cond
    load(datas{ci})  % load all Data tracks
    idvec = randperm(length(Data));
    data_ids(ci,:) = idvec(1:size(data_ids,2));
    Data = Data(data_ids(ci,:));
    data_id = 1:length(Data);  % original track ID
    ids{ci} = crossvalind('Kfold',data_id, K);  % indices for CV
end

%%% training MLEs
for ci = 1:cond  % C conditions
    %%% load a given condition
    load(datas{ci})  % load all Data tracks
    Data = Data(data_ids(ci,:)); % control the number of tracks 
    indices = ids{ci};  % indices for CV (pre-assigned!)
    
    for ki = 1:K  % K-fold
        % using only training set here (flipped the indices for more testing sets)
        test_set = (indices==ki);
        train_set = ~test_set;
        Data_train = Data(train_set); %(test_set);

        % training, with small repetition
        x_temp = zeros(rep,n_params);
        fv_temp = zeros(1,rep);
        for ri = 1:rep   % test with repreats to find the max MLE
            [x_MLE, fval] = MLE_mGLM(Data_train);
            fv_temp(ri) = fval;
            x_temp(ri,:) = x_MLE;
        end
        
        % record MLE
        [~,pos] = min(fv_temp);
        x_MLE = x_temp(pos,:);
        mle_params(ki, ci, :) = x_MLE;  % can replace this with repetition and choose best fit in the future..........
        [null_mle, null_fval] = MLE_randwalk(Data_train);  % random walk MLE
        null_params(ki,ci, :) = null_mle;
        [p_mle, p_fval] = MLE_binomial(Data_train);  % binomial chemotaxis MEL
        bin_params(ki,ci,1) = p_mle;
        
        kc_control_params(ki,ci,:) = MLE_wo_Kc(Data_train);
        kcp_control_params(ki,ci,:) = MLE_wo_Kcp(Data_train);
    
    end
end

% Naive Bayes for dC calculation
for ki = 1:K
    dc_train = [];
    datal_c = zeros(1,3);
    for ci = 1:cond  % C conditions
        load(datas{ci})  % load all Data tracks
        Data = Data(data_ids(ci,:)); % control the number of tracks 
        indices = ids{ci}; 
        test_set = (indices==ki);
        train_set = ~test_set;
        Data_train = Data(train_set);
        datal_c(ci) = sum(train_set);
        for dd = 1:datal_c(ci)
            dc_train = [dc_train   (Data_train(dd).dc(end) - Data_train(dd).dc(1))];
        end
    end
    labels = [make_label(datal_c(1),'app')', make_label(datal_c(2),'nai')', make_label(datal_c(3),'ave')']';
    mss = dc_train';
    Mdls{ki} = fitcnb(mss,labels,'ClassNames',{'app','nai','ave'});
end

%% %%%%%%%  TESTING  %%%%%%% 
%%% testing with MLEs
scal = 9;  % data length portions
min_scal = 0.35; %.35  % 0-1, for rescaling test data
rep_samp = 10;  %10% repeat sub-sampling for better average performance
cv_class = zeros(K, cond, scal, rep_samp); % K x c x T x s
data_len = zeros(K, cond, scal, rep_samp); % K x c x T x s
bin_class = zeros(K, cond, scal, rep_samp); % K x c x T x s
null_class = zeros(K, cond, scal, rep_samp); % K x c x T x s
cv_perf = zeros(cond, scal, rep_samp);  % c x T x s, this is averaged over K folds
cv_perf_bin = zeros(cond, scal, rep_samp);  % c x T x s,  for binomial comparison
cv_perf_null = zeros(cond, scal, rep_samp);  % c x T x s,  for null-model comparison
testLL = zeros(K,cond, 3);  % record the test LL (three conditions for analysis)

for ci = 1:cond  % C conditions
    %%% load a given condition
    load(datas{ci})  % load all Data tracks
    Data(data_ids(ci,:));
    indices = ids{ci};  % indices for CV (pre-assigned!)
    
    for ki = 1:K  % K-fold
        % using only training set here
        test_set = (indices==ki);
        train_set = ~test_set;
        Data_test = Data(test_set);%(train_set);
    
        % record the test LL
        x_MLE = squeeze(mle_params(ki, ci, :))';  % recall MLE
        [xx_test, yy_test, mask_test] = data2xy(Data_test);
        
        %%% testing null models
        info_norm_factor = (sum(mask_test)*(5/14)) * log(2);
        
        x_null = squeeze(null_params(ki, ci, :))';  % individually fitted null model parameters
        testLL(ki,ci,1) = -(pop_nLL(x_MLE, Data_test) - pop_nLL_randwalk(x_null, Data_test) ) / info_norm_factor;  % full model - randomwalk
%         testLL(ki,ci,1) = -(pop_nLL(x_MLE, Data_test) - pop_nLL(x_null, Data_test) ) / info_norm_factor;
        
        x_wo_Kc = [x_wo_Kc(1:2), zeros(1,4), x_wo_Kc(7:13)];
        testLL(ki,ci,2) = -(pop_nLL(x_MLE, Data_test) - pop_nLL(x_wo_Kc, Data_test) ) / info_norm_factor; % zero-out kernels
       
        x_wo_Kcp = [x_wo_Kcp(1:7),zeros(1,2), x_wo_Kcp(10:13)];
        testLL(ki,ci,3) = -(pop_nLL(x_MLE, Data_test) - pop_nLL(x_wo_Kcp, Data_test) ) / info_norm_factor;
        
        % testing scaled with time and bootstrap
        for ri = 1:rep_samp
            scal_vec = floor(linspace(1, min_scal*length(Data_test), scal));
            data_select_vec = [1:length(Data_test)];
            for si = 1:scal
                samps = randsample(data_select_vec, scal_vec(si));  % sample without replacement
                cv_class(ki, ci, si, ri) = argmaxLL(Data_test(samps), squeeze(mle_params(ki,:,:)));  % selection
                bin_class(ki, ci, si, ri) = argmaxLL_bin(Data_test(samps), squeeze(bin_params(ki,:,:)));  % compared to binomial model
                null_class(ki, ci, si, ri) = argmaxLL_null(Data_test(samps), squeeze(null_params(ki,:,:)));  % compared to null random-walk model
                [xx, yy, mm] = data2xy(Data_test(samps));  % concatenate data
                data_len(ki, ci, si, ri) = data_len(ki, ci, si, ri) + length(yy);  % average length(?)..........
            end
        end
    end
end

%%% lastly, quantify performance
for ci = 1:cond
    for si = 1:scal
        for ri = 1:rep_samp
            cv_perf(ci,si,ri) = length(find(cv_class(:, ci, si,ri)==ci))/K;
            cv_perf_bin(ci,si,ri) = length(find(bin_class(:, ci, si,ri)==ci))/K;
            cv_perf_null(ci,si,ri) = length(find(null_class(:, ci, si,ri)==ci))/K;
        end
    end
end

cv_nbc = zeros(K, scal, rep_samp);  % K x T x s
%%% Naiv Bayes testing
for ri = 1:rep_samp
for ki = 1:K
    dc_test = [];
    datal_c = zeros(1,3);
    for ci = 1:cond  % C conditions
        load(datas{ci})  % load all Data tracks
        Data = Data(data_ids(ci,:)); % control the number of tracks 
        test_set = (indices==ki);
        train_set = ~test_set;
        Data_test = Data(test_set);
        datal_c(ci) = sum(test_set);
        for dd = 1:datal_c(ci)
            dc_test = [dc_test   (Data_test(dd).dc(end) - Data_test(dd).dc(1))];
        end
    end
    labels = [make_label(datal_c(1),'app')', make_label(datal_c(2),'nai')', make_label(datal_c(3),'ave')']';
    scal_vec = floor(linspace(1, min_scal*length(Data_test), scal));
    data_select_vec = [1:length(dc_test)];
    for si = 1:scal
        samps = randsample(data_select_vec, scal_vec(si));  % sample without replacement
        predicted_labels = predict(Mdls{ki}, dc_test(samps)');
        cv_nbc(ki, si, ri) = sum(cellfun(@isequal, predicted_labels, labels(samps))) / (scal_vec(si));
    end
end
end
cv_nbc = squeeze(mean(cv_nbc,1)); 

%% plot cross-validated performance as a function of data length
figure;
rtime = squeeze(mean(mean(mean(data_len,1),2),4))*5/14/60/60;
figure;
plot(rtime, mean(mean(cv_perf,3)),'k-o')
hold on
errorbar(rtime, mean(mean(cv_perf,3)), std(mean(cv_perf,3))/sqrt(K), '.k')

plot(rtime, mean(mean(cv_perf_bin,3)),'g--')
errorbar(rtime, mean(mean(cv_perf_bin,3)), std(mean(cv_perf_bin,3))./sqrt(K.*scal_vec), '.g')

plot(rtime, mean((cv_nbc),2),'r--')
errorbar(rtime, mean((cv_nbc),2), std((cv_nbc),0,2)'./sqrt(K.*scal_vec), '.r')

plot(rtime, mean(mean(cv_perf_null,3)),'b-*')
errorbar(rtime, mean(mean(cv_perf_null,3)), std(mean(cv_perf_null,3))./sqrt(K.*scal_vec), '*b')

xlabel('mean data length (hr)')
ylabel('cross-validation (ratio)')
% legend({'appetitive','naive','aversive','mean prediction', 'binomial prediction', '\Delta C', 'behavior'})
legend({'mean prediction','','binomial prediction', '','\Delta C', '','behavior',''})
set(gcf,'color','w'); set(gca,'Fontsize',20);

figure
for ii = 1:3
    temp_cv = squeeze(cv_perf(ii,:,:));
    plot(rtime, mean(temp_cv,2),'-o')
    hold on
    errorbar(rtime, mean(temp_cv,2), std(temp_cv,[],2)/sqrt(K), '.k')
end
legend({'appetitive','','naive','','aversive',''})
xlabel('mean data length (hr)')
ylabel('cross-validation (ratio)')
set(gcf,'color','w'); set(gca,'Fontsize',20);

%% some post analysis for variability!
ttl = {'appetitive','naive','aversive'};
figure
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
for cc = 1:3
    subplot(1,3,cc);
for ii = 1:K
    temp = squeeze(mle_params(ii,cc,3:6))'*cosBasis';
    plot(tt,temp)
    hold on
end; title(ttl{cc})
end

figure; 
xx = 0:length(cosBasis)-1;
for cc = 1:3
    subplot(1,3,cc)
for ii = 1:K
    amp = mle_params(ii,cc,8); tau = mle_params(ii,cc,9); temp = (-amp*exp(-xx/tau));
    plot(tt,temp)
    hold on
end; title(ttl{cc})
end

%% information analysis
figure
col = {'b','k','r'};
tils = {'full','K_c','K_{c^\perp}'};
for nn = 1:3
    subplot(1,3,nn)
    temp_m = squeeze(mean(testLL(:,:,nn),1));
    temp_s = squeeze(std(testLL(:,:,nn),0,1))/sqrt(K);
    for cc = 1:3
        bar(cc, temp_m(cc), 'FaceColor',col{cc})
        hold on
        errorbar(cc, temp_m(cc), temp_s(cc), 'k.');
    end
    title(tils{nn})
    xticklabels([]);
end

%% functionals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_MLE, fval] = MLE_mGLM(Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);  %nLL_kernel_hist2
    opts = optimset('display','iter');
    LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.2];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 50, 1];
    prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end

function [x_MLE, fval] = MLE_randwalk(Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_randomwalk(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    opts = optimset('display','iter');
    LB = [1e-0, 0, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.1];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 50, 1];
    prs0 = [50, 0.05, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end

function [x_MLE, fval] = MLE_binomial(Data_train)
    lfun = @(x)pop_nLL_binomial(x, Data_train);
    opts = optimset('display','iter');
    LB = 0;
    UB = 1;
    prs0 = 0.5;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end

function [x_MLE] = MLE_wo_Kc(Data_train);
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_wo_kc(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    opts = optimset('display','iter');
    LB = [1e-0, 0, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.2];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 50, 1];
    prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end

function [x_MLE] = MLE_wo_Kcp(Data_train);
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_wo_kcp(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    opts = optimset('display','iter');
    LB = [1e-0, 0, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.2];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 50, 1];
    prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end

function [NLL] = pop_nLL_kc(THETA, D)
    [xx_train, yy_train, mask_train] = data2xy(D);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    NLL = nLL_wo_kc(THETA, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
end

function [NLL] = pop_nLL_kcp(THETA, D)
    [xx_train, yy_train, mask_train] = data2xy(D);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    NLL = nLL_wo_kcp(THETA, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
end

function [NLL] = pop_nLL_rand2(THETA, D)
    [xx_train, yy_train, mask_train] = data2xy(D);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    NLL = nLL_randomwalk(THETA, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
end


function predict_lambda = argmaxLL(data, mle_params)
%%%
% max-likelihood classification with the full mixture GLM model
%%%
    [Basis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    [xx, ang_fit, trials_fit] = data2xy(data);  % turn data sturcute into vectors for plug-in evaluation
    dcp_fit = xx(2,:);  % dcp concatented
    ddc_fit = xx(1,:);  % dc concatentated
    lls = zeros(1,3);  % three conditions
    for c = 1:3
        lls(c) = -nLL_kernel_hist2(mle_params(c,:), ang_fit, dcp_fit, ddc_fit, Basis, .1, trials_fit) / length(xx);  % normalize by length?
        
    end
%     lls
    predict_lambda = argmax(lls);
end

function predict_lambda = argmaxLL_bin(data, bin_params)
%%%
% max-likelihood classification with the binomial chemotaxis model
%%%
    lls = zeros(1,3);  % three conditions
    for c = 1:3
        lls(c) = -pop_nLL_binomial(bin_params(c), data);
    end
    predict_lambda = argmax(lls);
end

function predict_lambda = argmaxLL_null(data, bin_params)
%%%
% max-likelihood classification with the null behavioral model
%%%
    lls = zeros(1,3);  % three conditions
    for c = 1:3
        lls(c) = -pop_nLL_randwalk(bin_params(c,:), data);
    end
    predict_lambda = argmax(lls);
end

function labels= make_label(ns, label)
%%%
% make ns labels for naive Bayes classification
%%%
    labels = cell(ns,1);
    for ii = 1:ns
        labels{ii} = label;
    end
end