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
specs.T = floor(30*60*14/5);
specs.dt = 1;
specs.REP = 50;

%% simulation loop
kfold = size(mle,1);
CIs = zeros(3, kfold);
BWs = zeros(3,kfold, 2);
for ci = 1:3
    for ki = 1:kfold
        x = squeeze(median(mle(:,ci,:),1))';  %squeeze(mle(ki,ci,:))';
        [tracks, CI] = param2tracks(x, specs, []);  % check if alldis is given from data
        
        [brw_index, wv_index] = track2strat(tracks, M);
        
        CIs(ci,ki) = CI;
        BWs(ci,ki,:) = [brw_index, wv_index];
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

%% visualize tracks
figure; imagesc(M); hold on
for ii=1:specs.REP
plot(tracks(ii).xy(1,:),tracks(ii).xy(2,:))
hold on
end

%% example tracks
specs.REP = 1;
figure; imagesc(M); hold on
ki = 4;
for ci = 1:3
    x = squeeze(mle(ki,ci,:))';
    [tracks, CI] = param2tracks(x, specs, []);
    plot(tracks(1).xy(1,:),tracks(1).xy(2,:))
    hold on
end

%% check run distribution form simulated tracks
load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat');
cond_id = 1;
specs.REP = 100;
x = squeeze(mle(1,cond_id,:))';
x_temp= x;
% x_temp(10) = -x(10)*0;
[tracks, CI] = param2tracks(x_temp, specs, []);
%% plots
thre = 50;
run_dur_model = [];
run_dur_data = [];
for ii = 1:100
    dth_i = tracks(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_model = [run_dur_model   diff(pos)];
    dth_i = Data(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_data = [run_dur_data   diff(pos)];
end

figure
subplot(121)
[xx_train, yy_train, mask_train] = data2xy(tracks); %tracks or Data
hist(yy_train,100)
subplot(122)
hist(run_dur_model,50);
figure
subplot(121)
[xx_train, yy_train, mask_train] = data2xy(Data);
hist(yy_train,100)
subplot(122)
hist(run_dur_data*14/5,50);

%%
figure
subplot(121)
[xx_train, yy_train, mask_train] = data2xy(tracks); %tracks or Data
histogram(yy_train,100, 'Normalization', 'probability',  'FaceAlpha', 0.7); hold on
[xx_train, yy_train, mask_train] = data2xy(Data);
histogram(yy_train,100, 'Normalization', 'probability',  'FaceAlpha', 0.7)
subplot(122)
numBins = 100;
binEdges = linspace(min([run_dur_model run_dur_data*14/5]), max([run_dur_model run_dur_data*14/5]), numBins + 1);
histogram(run_dur_model, binEdges, 'Normalization', 'probability', 'FaceAlpha', 0.5); hold on
histogram(run_dur_data*14/5, binEdges, 'Normalization', 'probability',  'FaceAlpha', 0.5);

%%
[dur, ydata] = hist(run_dur_data,150);
[dur, ydata] = hist(run_dur_model,150);
% dur = max([run_dur_model run_dur_data*14/5]);
y = ydata;
F = @(x,dur)x(1)*exp(-dur/x(2)) + x(3)*exp(-dur/x(4));
x0 = [100 10 1 50] ;
xunc = lsqcurvefit(F, x0, dur, y);
tlist = linspace(min(dur), max(dur));   % Plot Finer Resolution
tc = 1/(1/xunc(2)-1/xunc(4)) * log(xunc(1)/xunc(3)) / (5/14)
