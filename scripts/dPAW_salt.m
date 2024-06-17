% dPAW_salt
%%% for revision, we apply dPAW to salt data
%%% here we load in Data for 0-50 and 100-50 mM salt chemotaxis
%%% load in the fitted dPAW parameters runned with Chemotaxis_model_pop.m
%%% script and saved in the folder
%%% We compare kernels, tracks, and simulated CI

%% load fits and data
temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/dPAW_salt50_100_0511.mat');
x_100 = temp.x;
data100_fit = temp.Data_fit;
temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/dPAW_salt0_50_0513_2.mat');
x_0 = temp.x;
data0_fit = temp.Data_fit;

%% process kernels
xxs = {x_0, x_100};
figure;
for ii = 1:2
    [K_dc_rec, K_dcp_rec] = dPAW_x2k(xxs{ii});
    tt = [1:length(K_dc_rec)].*5/14;
    subplot(121)
    plot(tt, K_dc_rec); hold on
    xlim([0,12]); xlabel('time (s)'); ylabel('weights'); set(gcf,'color','w'); set(gca,'Fontsize',20);
    subplot(122)
    plot(tt, K_dcp_rec); hold on
    xlim([0,12]); xlabel('time (s)'); ylabel('weights'); set(gcf,'color','w'); set(gca,'Fontsize',20);
end

%% compute CI for data
temp = load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt100_50_0511.mat');
data_100_50 = temp.Data;
temp = load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50.mat');
data_0_50 = temp.Data;

CI_data = zeros(1,2);
CI_data(1) = Data2CI(data0_fit);
CI_data(2) = Data2CI(data100_fit);
CI_data

%% sim tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% landscape
rows=2500; cols=3000;
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
gradient_x = x_ * 1;
M0 = (y_*0+1) .* gradient_x;
[x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M100 = (y_*0+1) .* gradient_x;

%% specs for sim
rng(42)
clear specs
specs = struct();
specs.M = M0;
specs.fr = 14/5;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
specs.cosBasis = cosBasis;
specs.T = floor(30*60*14/5);
specs.dt = 1;
specs.REP = 100;
reps = 5;

CI_sim = zeros(reps, 2);

for rr = 1:reps
    rr
    [tracks, CI] = param2tracks(x_0, specs, []);
    CI_sim(rr, 1) = CI;
    specs.M = M100;
    [tracks, CI] = param2tracks(x_100, specs, []);
    CI_sim(rr, 2) = CI;
end

%% plot results
figure;
bar([1,2], mean(CI_sim,1)); hold on
errorbar([1,2], mean(CI_sim,1), std(CI_sim,1,1),'.k')
plot([1,2], CI_data, 'ro')
names = {'0-50'; '50-100'};
ylabel('CI')
set(gca,'xtick',[1:2],'xticklabel',names,'FontSize',20); set(gcf,'color','w');

%% showcase tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data
rng(18) %37
nexp = 15;
temp = randperm(100);
selected = temp(1:nexp);
figure;
subplot(121)
imagesc(M0); colorbar(); hold on
for ii = 1:nexp
    temp = data_0_50(selected(ii)).xy;
    plot(temp(1,:), temp(2,:),'k', 'LineWidth',1.5); hold on
    plot(temp(1,1), temp(2,1), 'g.','MarkerSize',13); plot(temp(1,end), temp(2,end), 'r.','MarkerSize',13);
end
set(gcf,'color','w'); set(gca,'Fontsize',20);
rng(5)
temp = randperm(100);
selected = temp(1:nexp);
subplot(122)
imagesc(M100); colorbar(); hold on
for ii = 1:nexp
    temp = data_100_50(selected(ii)).xy;
    plot(temp(1,:), temp(2,:),'k', 'LineWidth',1.5); hold on
    plot(temp(1,1), temp(2,1), 'g.','MarkerSize',13); plot(temp(1,end), temp(2,end), 'r.','MarkerSize',13);
end
set(gcf,'color','w'); set(gca,'Fontsize',20);

%% model
rng(42)
nexp = 15;
figure;
subplot(121)
imagesc(M0); colorbar(); hold on
specs.M = M0;
specs.REP = nexp;
[tracks, CI] = param2tracks(x_0, specs, []);
for ii = 1:nexp
    temp = tracks(ii).xy;
    plot(temp(1,:), temp(2,:),'k', 'LineWidth',1.5); hold on
    plot(temp(1,1), temp(2,1), 'g.','MarkerSize',13); plot(temp(1,end), temp(2,end), 'r.','MarkerSize',13);
end
set(gcf,'color','w'); set(gca,'Fontsize',20);
subplot(122)
imagesc(M100); colorbar(); hold on
specs.M = M100;
specs.REP = nexp;
[tracks, CI] = param2tracks(x_100, specs, []);
for ii = 1:nexp
    temp = tracks(ii).xy;
    plot(temp(1,:), temp(2,:),'k', 'LineWidth',1.5); hold on
    plot(temp(1,1), temp(2,1), 'g.','MarkerSize',13); plot(temp(1,end), temp(2,end), 'r.','MarkerSize',13);
end
set(gcf,'color','w'); set(gca,'Fontsize',20);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function
function [K_dc_rec, K_dcp_rec] = dPAW_x2k(x);
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7); Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); K2_ = x(12);  gamma = x(13);
    base_dc = 0;  base_dcp = 0;
    xx = 0:length(cosBasis)-1;
    K_dcp_rec = Amp*exp(-xx/tau);
    K_dc_rec = B_*cosBasis';
end

function [ci] = Data2CI(Data)
    cii = 0;
    for ii = 1:length(Data)
        temp = Data(ii).dc;
        if temp(end)>temp(1)
            cii = cii + 1;
        end
    end
    ci = 2*cii/length(Data) - 1;
end
