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
bar([1:length(fname)],[0.09,0.2984,0.50,0.57,2/3])
set(gca, 'XTick', 1:length(fname),'XTickLabel',fname);
ylabel('cross-validation loss')
xlim([.5,length(fname)+.5])
ylim([0, 2/3+0.05])
%%
function labels= make_label(ns, label)
    labels = cell(ns,1);
    for ii = 1:ns
        labels{ii} = label;
    end
end