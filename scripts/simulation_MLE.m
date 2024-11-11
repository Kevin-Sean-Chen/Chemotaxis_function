% simulation_MLE
% script to simulate chemotaxis tracks and calculate chemotaxis index with
% the fitted MLE

%% load array of fitted MLE
mle = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param3.mat');
% mle = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param7.mat');
mle = mle.mle_params;

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

%% analysis
%%%%%% maybe include index for two strategies later?
figure;
bar(CIs)

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

%%
xx = [1 2 3; 5 6 7];
figure
% bar(squeeze(mean(BWs,2))'*1)
bar(xx,squeeze(mean(BWs,2))'*1); hold on
errorbar(xx, squeeze(mean(BWs,2))', squeeze(std(BWs,[],2))'/sqrt(kfold),'o')

%% visualize tracks
figure; imagesc(M); hold on
for ii=1:specs.REP
plot(tracks(ii).xy(1,:),tracks(ii).xy(2,:))
hold on
end

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

%% check run distribution form simulated tracks
load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat');
cond_id = 1;
specs.REP = 100;
x = squeeze(mle(2,cond_id,:))';
x_temp= x;
% x_temp(10) = -x(10)*0;
[tracks, CI] = param2tracks(x_temp, specs, []);
x_temp(10) = 0;  % without history
[tracks_woh, CI] = param2tracks(x_temp, specs, []);

%% plots
thre = 50;
run_dur_model = [];
run_dur_data = [];
run_dur_woh = [];
for ii = 1:100
    dth_i = tracks(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_model = [run_dur_model   diff(pos)];
    dth_i = tracks_woh(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_woh = [run_dur_woh   diff(pos)];
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
% subplot(121)
% [xx_train, yy_train, mask_train] = data2xy(Data);
% histogram(yy_train,100, 'Normalization', 'probability',  'FaceAlpha', 0.7); hold on
% [xx_train, yy_train, mask_train] = data2xy(tracks); %tracks or Data
% histogram(yy_train,100, 'Normalization', 'probability',  'FaceAlpha', 0.7);
% [xx_train, yy_train, mask_train] = data2xy(tracks_woh);
% histogram(yy_train,100, 'Normalization', 'probability',  'FaceAlpha', 0.7)
% subplot(122)
numBins = 40;
% binEdges = linspace(min([run_dur_model*0 run_dur_data*14/5]), max([run_dur_model*0 run_dur_data*14/5]), numBins + 1);
binEdges = linspace(0, 120, numBins + 1);
histogram(run_dur_data*14/5, binEdges, 'Normalization', 'probability',  'FaceAlpha', 0.5); hold on
histogram(run_dur_model, binEdges, 'Normalization', 'probability', 'FaceAlpha', 0.5);
histogram(run_dur_woh, binEdges, 'Normalization', 'probability',  'FaceAlpha', 0.5);

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

%% check concentration distributions
rng(1)
% load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat');
temp = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param7.mat');
mle = temp.mle_params;
cond_id = 2;
specs.REP = 50;
x = squeeze(mle(4,cond_id,:))';
x_temp= x;
% x_temp(10) = -x(10)*0;
specs.T = floor(30*60*14/5);
[tracks, CI] = param2tracks(x_temp, specs, []);
c_sim = extractfield(tracks, 'dc');
dcp_sim = extractfield(tracks, 'dcp');
load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat');
c_data = extractfield(Data, 'dc');
dcp_data = extractfield(Data, 'dcp');

%%
nb = 20;
figure;
subplot(121)
mm = min([c_sim, c_data]); MM = max([c_sim, c_data]);
beds = linspace(mm,MM, nb);
histogram(c_data, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'black','FaceAlpha', 0.5);
hold on;
histogram(c_sim, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 0.5);
xlabel('ppm'); ylabel('probability'); set(gcf,'color','w'); set(gca,'Fontsize',20);
subplot(122)
mm = min([dcp_sim, dcp_data]); MM = max([dcp_sim, dcp_data]);
beds = linspace(mm,MM, nb);
histogram(dcp_data, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'black','FaceAlpha', 0.5);
hold on;
histogram(dcp_sim, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 0.5);
xlabel('ppm/mm'); ylabel('probability'); set(gcf,'color','w'); set(gca,'Fontsize',20);
