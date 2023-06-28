% simulation_MLE
% script to simulate chemotaxis tracks and calculate chemotaxis index with
% the fitted MLE

%% load array of fitted MLE
% mle = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param3.mat');
% mle = mle.mle_params;

%% specify simulation parameters
clear specs
specs = struct();

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
M = Cmap.vq1;
M = fliplr(flipud(M));
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);

specs.M = M;
specs.fr = 14/5;
specs.cosBasis = cosBasis;
specs.T = floor(20*60*14/5);
specs.dt = 1;
specs.REP = 50;

%% simulation loop
kfold = size(mle,1);
CIs = zeros(3, kfold);
for ci = 1:3
    for ki = 1:kfold
        x = squeeze(mle(ki,ci,:))';
        [tracks, CI] = param2tracks(x, specs, alldis);  % check if alldis is given from data
        
        CIs(ci,ki) = CI;
        ci
        ki
    end
end

%% analysis
%%%%%% maybe include index for two strategies later?
figure;
bar(CIs)

figure;
ctr = [1,2,3];
hBar = bar(mean(CIs'));
for kl = 1:3
    bar(ctr(kl), mean(CIs(kl,:)))
    hold on
    errorbar(ctr(kl), mean(CIs(kl,:)), std(CIs(kl,:))/sqrt(kfold), '.k')               
end
names = {'appetitive'; 'naive'; 'aversive'};
ylabel('CI')
set(gca,'xtick',[1:3],'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
title('simulation with MLE parameters')

%%
figure; imagesc(M); hold on
for ii=1:specs.REP
plot(tracks(ii).xy(1,:),tracks(ii).xy(2,:))
hold on
end