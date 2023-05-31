% classifiy_learn
%%% classify learning conditions given parameter samples
%%% with proper K-fold cross-validation in earch loop

%% load tracks from each condition
% load data files
datas = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave.mat'};
     
%% loop through K-folds, conditions, and scaling
rep = 1;  % repetition per fold to make sure that it is the better fit
K = 10;  % K-fold cross-validation
scal = 5;  % data length portions
min_scal = 0.2;  % 0-1, for rescaling test data
cond = length(datas);  % all exp conditions
n_params = 13;  % number of parameters in our model for now
mle_params = zeros(K, cond, n_params); % K x c x N
cv_class = zeros(K, cond, scal); % K x c x T
data_len = zeros(K, cond, scal); % K x c x T
cv_perf = zeros(cond, scal);  % c x T, this is averaged over K folds

%%% K x c x T times optimization  %%%
for ki = 1:K  % K-fold loop
    
    for ci = 1:cond  % experimental conditions
        load(datas{ci})  % load all Data tracks
        data_id = 1:length(Data);  % original track ID
        indices = crossvalind('Kfold',data_id, K);  % indices for CV
        
        % define train and test sets
        test_set = (indices==ki);
        train_set = ~test_set;
        
        % load individual tracks and concatenation... changing this so
        % there is more testing to scan through...
        Data_train = Data(test_set);
        Data_test = Data(train_set);
        
        % training
        [x_MLE, fval] = MLE_mGLM(Data_train);
        rep;
        mle_params(ki, ci, :) = x_MLE;  % can replace this with repetition and choose best fit in the future..........
        
        % testing scaled with time
        scal_vec = fliplr(floor(min_scal./[1:scal].*length(Data_test)));  % scaled data length, can do proper tiled timing sampling in the future..........
        data_select_vec = [1:length(Data_test)];
        for si = 1:scal
            samps = randsample(data_select_vec, scal_vec(si));  % sample without replacement
            cv_class(ki, ci, si) = argmaxLL(Data_test(samps), squeeze(mle_params(ki,:,:)));  % selection
            [xx, yy, mm] = data2xy(Data_test(samps));  % concatenate data
            data_len(ki, ci, si) = data_len(ki, ci, si) + length(yy);  % average length(?)..........
        end
    end    
    
        
end

for ci = 1:cond
    for si = 1:scal
        cv_perf(ci,si) = length(find(cv_class(:, ci, si)==ci))/K;
    end
end

%%
figure;
tlt = data_len; %squeeze(mean(data_len,2));
plot(tlt(:)*5/14/60/60,cv_class(:),'o')
xlabel('data length (hour)')
ylabel('class')

rtime = squeeze(mean(mean(data_len,1),2))*5/14/60/60;
figure;
plot(rtime, cv_perf','-o')
hold on
plot(rtime, mean(cv_perf),'k-o')
xlabel('mean data length (hr)')
ylabel('cross-validation (ratio)')
legend({'appetitive','naive','aversive','mean prediction'})
set(gcf,'color','w'); set(gca,'Fontsize',20);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classify with log-likelihood!
% load mle parameters
% mle_params = zeros(3,13);
% test = load('/projects/LEIFER/Kevin/Data_learn/app_param_new.mat');
% mle_params(1,:) = test.x;
% test = load('/projects/LEIFER/Kevin/Data_learn/nai_param_new.mat');
% mle_params(2,:) = test.x;
% test = load('/projects/LEIFER/Kevin/Data_learn/ave_param_new.mat');
% mle_params(3,:) = test.x;
%      
% % CV settings
% rep = 20;  % repeat cross-validations
% scal = 5;  % data length portions
% cv_class = zeros(rep,scal);  % repeats x data length
% datals = cv_class*1;  % record actual data length
% cv_perf = zeros(3,scal);  % record performance (%) across repeats
% 
% for dd = 1:3
%     load(datas{dd})
%     Datai = Data(1:end);  % load Data structure
%     scal_vec = fliplr(floor(.5./[1:scal].*length(Datai)));  % scaled data length
%     data_select_vec = [1:length(Datai)];  % select from track ID
% for rr = 1:rep
%     for sc = 1:scal
%         samps = randsample(data_select_vec, scal_vec(sc));  % sample without replacement
%         cv_class(rr, sc) = argmaxLL(Datai(samps), mle_params, cosBasis);  % selection
%         [xx, yy, mm] = data2xy(Datai(samps));  % concatenate data
%         datals(rr, sc) = length(yy);
%     end    
% end
% 
% for sc = 1:scal
%     cv_perf(dd,sc) = length(find(cv_class(:,sc)==dd))/rep;
% end
%     
% end
% 
% cv_class

%%
function [x_MLE, fval] = MLE_mGLM(Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
    opts = optimset('display','iter');
    LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.1];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1];
    prs0 = [50, 0.5, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end
function predict_lambda = argmaxLL(data, mle_params)
    [Basis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    [xx, ang_fit, trials_fit] = data2xy(data);  % turn data sturcute into vectors for plug-in evaluation
    dcp_fit = xx(2,:);  % dcp concatented
    ddc_fit = xx(1,:);  % dc concatentated
    lls = zeros(1,3);  % three conditions
    for c = 1:3
        lls(c) = -nLL_kernel_hist2(mle_params(c,:), ang_fit, dcp_fit, ddc_fit, Basis, .0, trials_fit) / length(xx);  % normalize by length?
    end
%     lls
    predict_lambda = argmax(lls);
end

function labels= make_label(ns, label)
    labels = cell(ns,1);
    for ii = 1:ns
        labels{ii} = label;
    end
end