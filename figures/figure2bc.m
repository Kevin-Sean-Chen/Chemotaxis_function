% Figure 2b,c
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2b
% extracted from script 'simulation_MLE.m'

%% load data
load(fullfile(datadir,'data4plots', 'Kfold_mle_param7.mat'));
Cmap = load(fullfile(datadir,'data4plots', 'Landscape_low_0623_2.mat'));
M = Cmap.vq1;
M = fliplr(flipud(M));
mle = mle_params;

%% set configuration
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
specs.M = M;
specs.fr = 14/5;
specs.cosBasis = cosBasis;
specs.T = floor(30*60*14/5);
specs.dt = 1;
specs.REP = 100;

%% simulation loop
kfold = size(mle,1);
CIs = zeros(3, kfold);
BWs = zeros(3,kfold, 2);
for ci = 1:3
    for ki = 1:kfold
        x = squeeze(mle(ki,ci,:))';  %squeeze(median(mle(:,ci,:),1))';  % 
        [tracks, CI] = param2tracks(x, specs, []);  % check if alldis is given from data
        [brw_index, wv_index] = track2strat(tracks, M);
        CIs(ci,ki) = CI;
        BWs(ci,ki,:) = [brw_index, wv_index];
        ci
        ki
    end
end

%% plotting
% raw values
figure;
bar(CIs)
names = {'appetitive'; 'naive'; 'aversive'};
set(gca,'xticklabel',names,'FontSize',20); set(gcf,'color','w');  % compare this with experimental results in fig.1d

% stats
figure;
ctr = [1,2,3];
hBar = bar(median(CIs'));
for kl = 1:3
    bar(ctr(kl), median(CIs(kl,:)))
    hold on
    errorbar(ctr(kl), median(CIs(kl,:)), std(CIs(kl,:))/sqrt(kfold), '.k')               
end
names = {'appetitive'; 'naive'; 'aversive'};
ylabel('CI')
set(gca,'xtick',[1:3],'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
title('simulation with MLE parameters')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2c
% extracted from script 'simulation_MLE.m'

%% load data
load(fullfile(datadir,'data4plots', 'Kfold_mle_param3.mat'));
mle = mle_params;

%% example tracks
rng(15)
specs.REP = 1;
figure; imagesc(M); hold on
ki = 4;
for ci = 1:3
    x = squeeze(mle(ki,ci,:))';
    [tracks, CI] = param2tracks(x, specs, []);  % comment out random initial points for this simulation
    plot(tracks(1).xy(1,:),tracks(1).xy(2,:))
    hold on
    plot(tracks(1).xy(1,1),tracks(1).xy(2,1), 'go')
    plot(tracks(1).xy(1,end),tracks(1).xy(2,end), 'ro')
end
