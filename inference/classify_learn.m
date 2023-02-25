% classifiy_learn
%%% classify learning conditions given parameter samples

%% load mcmc sample parameters
test = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/app_mcmc.mat');
m_app = test.m_bound;
test = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/nai_mcmc.mat');
m_nai = test.m_bound;
test = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/ave_mcmc.mat');
m_ave = test.m_bound;

%% make labels
n_ap = size(m_app,2);
n_na = size(m_nai,2);
n_av = size(m_ave,2);
labels = [make_label(n_ap,'app')', make_label(n_na,'nai')', make_label(n_av,'ave')']';

mss = [m_app(8,:,1)'; m_nai(8,:,1)'; m_ave(8,:,1)'];  %Kdc^p
% mss = [vecnorm(m_app(3:6,:,1),2,1)'; vecnorm(m_nai(3:6,:,1),2,1)'; vecnorm(m_ave(3:6,:,1),2,1)']; %K_dc

%%
Mdl = fitcnb(mss,labels,'ClassNames',{'app','nai','ave'})
CVMdl = crossval(Mdl);
Loss = kfoldLoss(CVMdl)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classify with dC
datal = 200; %200 or 100000
test = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/app_ddc.mat');  % _dc or _ddc
dc_app = test.alldeltaC(1:datal);  %alldC or alldeltaC
test = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/nai_ddc.mat');
dc_nai = test.alldeltaC(1:datal);
test = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/ave_ddc.mat');
dc_ave = test.alldeltaC(1:datal);

labels = [make_label(datal,'app')', make_label(datal,'nai')', make_label(datal,'ave')']';
mss = [dc_app  dc_nai  dc_ave]';

%%
Mdl = fitcnb(mss,labels,'ClassNames',{'app','nai','ave'})
CVMdl = crossval(Mdl);
Loss = kfoldLoss(CVMdl)

%%
figure;
fname = {'K_{dC^{^\perp}}','K_{dC}','\Delta C','C','random'}; 
bar([1:length(fname)],(1-[0.09,0.2984,0.50,0.57,2/3])*100)
set(gca, 'XTick', 1:length(fname),'XTickLabel',fname);
ylabel('cross-validation performance (%)')
xlim([.5,length(fname)+.5])
% ylim([0, 2/3+0.05])

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classify with log-likelihood!
% load mle parameters
mle_params = zeros(3,15);
test = load('/projects/LEIFER/Kevin/Data_learn/app2_param.mat');
mle_params(1,:) = test.x;
test = load('/projects/LEIFER/Kevin/Data_learn/nai4_param.mat');
mle_params(2,:) = test.x;
test = load('/projects/LEIFER/Kevin/Data_learn/ave2_param.mat');
mle_params(3,:) = test.x;

datas = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave.mat'};
rep = 10;  % repeat cross-validations
scal = 5;  % data length portions
data_select_vec = [1:length(Data)];  % select from track ID
scal_vec = fliplr(floor(0.1./[1:scal].*length(Data)));
cv_class = zeros(rep,scal);  % repeats x data length
datals = cv_class*1;  % record actual data length
for dd = 1:3
    load(datas{dd})
for rr = 1:rep
    for sc = 1:scal
        samps = randsample(data_select_vec, scal_vec(sc));
        cv_class(rr, sc) = argmaxLL(Data(samps), mle_params, cosBasis);  % selection
        [xx, yy, mm] = data2xy(Data(samps));
        datals(rr, sc) = length(yy);
    end
end
end

cv_class
%%
figure;
plot(datals(:)*5/14/60/60,cv_class(:),'o')
xlabel('data length (hour)')
ylabel('class')

%%
function predict_lambda = argmaxLL(data, mle_params,Basis)
    [xx, ang_fit, trials_fit] = data2xy(data);  % turn data sturcute into vectors for plug-in evaluation
    dcp_fit = xx(2,:);  % dcp concatented
    ddc_fit = xx(1,:);  % dc concatentated
    lls = zeros(1,3);  % three conditions
    for c = 1:3
        lls(c) = -nLL_kernel_hist2(mle_params(c,:), ang_fit, dcp_fit, ddc_fit, Basis, .0, trials_fit);
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