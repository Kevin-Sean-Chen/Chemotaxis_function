% Figure 3
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 c
% extracted from script 'exp_track_NaCl.m'
% example tracks for salt chemotaxis

%% load data
load(fullfile(datadir,'data4plots', 'Data_salt100_50.mat'))
% load(fullfile(datadir,'data4plots', 'Data_salt0_50_0513.mat'))

%% for salt linear gradient
Cmap = load(fullfile(datadir,'data4plots', 'Landscape_low_0623_2.mat'));
M = Cmap.vq1;
M = fliplr(flipud(M));
[rows, cols] = size(M);
[x_, y_] = meshgrid(linspace(50, 0, cols), 1:rows);  % for   0 to 50 mM
[x_, y_] = meshgrid(linspace(50, 100, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M = (y_*0+1) .* gradient_x;  
%figure; imagesc(M)
pix2mm = 1/30;

%% example tracks
data_i = Data([57:77]);
% data_i = Data(10:30);

%%% plotting
figure()
ax1 = axes;
imagesc(ax1,M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
set(gcf,'color','w'); set(gca,'Fontsize',20); xticks([]); yticks([]);
colormap()
caxis([0 100]);
hold on
ax2 = axes;
for ii = 1:length(data_i)
    xx = data_i(ii).xy(1,:);
    yy = data_i(ii).xy(2,:);
    ll = length(data_i(ii).xy);
    gg = linspace(0,1,ll);
%     plot(xx, yy, 'Color', [grayLevel grayLevel grayLevel]);
    patch(ax2, [xx nan]*pix2mm,[yy nan]*pix2mm,[gg nan],[gg nan], 'edgecolor', 'interp','LineWidth',2); 
    hold on
    plot(ax2,xx(1)*pix2mm, yy(1)*pix2mm,'g.', 'MarkerSize',25)
    plot(ax2,xx(end)*pix2mm, yy(end)*pix2mm,'r.', 'MarkerSize',25)
end
set(gca, 'YDir','reverse')
set(gca, 'XDir', 'reverse')  %%% for 100-50
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
c = gray;
colormap(ax2,c)
colormap(ax1)
xlim([1,3000].*pix2mm); ylim([1,2500].*pix2mm)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 d
% extracted from script 'index_salt.m'
% compute indeces for salt chemotaxis

%% linear salt environment
rows=2500; cols=3000;
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
gradient_x = x_ * 1;
M0 = (y_*0+1) .* gradient_x;
[x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M100 = (y_*0+1) .* gradient_x;

[ci_, brw_index, wv_index] = compute_index(Tracks(1:end), M0, 15);

%% looping files and conditions
rng(123)
fields_to_load = {'Path','Time','Runs','Pirouettes'};
temp = load(fullfile(datadir,'data4plots', 'salt_data_folder.mat'));
folder_all = {temp.salt_data_50_0, temp.salt_data_50_100};
track_learn_all = cell(1,2);
BWCs = cell(1,2);
cols = ['b','k','r'];
rep_samp = 10;
n_samps = 500;
figure
for cc = 1:2
    folder_names = folder_all{cc};
    BWC = zeros(length(folder_names), 3);
    sub_learn_tracks = cell(1,length(folder_names));
    Tracks = loadtracks(folder_names,fields_to_load);
    for ff = 1:rep_samp
        samp_id = randperm(length(Tracks));
        Tracks_i = Tracks(samp_id(1:n_samps));
        if cc == 1
            [ci_, brw_index, wv_index] = compute_index(Tracks_i, M0, 20);
            BWC(ff,:) = [ci_, brw_index, wv_index];
        else
            [ci_, brw_index, wv_index] = compute_index(Tracks_i, M100, 15);
            BWC(ff,:) = [ci_, brw_index, -wv_index];
        end
        
        sub_learn_tracks{ff} = Tracks;
    end
    
    plot(BWC',cols(cc)); hold on
    BWCs{cc} = BWC;
    track_learn_all{cc} = sub_learn_tracks;
end

%% bar plot
figure
m_0 = mean(BWCs{1},1);
m_100 = mean(BWCs{2},1);
s_0 = std(BWCs{1},[],1)/sqrt(rep_samp);
s_100 = std(BWCs{2},[],1)/sqrt(rep_samp);

hBar = bar([m_0; m_100]');

for k1 = 1:2
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');     
    ydt(k1,:) = hBar(k1).YData;  hold on
end
hold on 
errorbar(ctr, ydt, [m_0;m_100]*0,[s_0;s_100], '.k')         
hold off
ylabel('Index')
names = {'CI'; 'BRW'; 'WV'};
set(gca,'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
legend([hBar(1), hBar(2)], '0-50','50-100')
ylim([-1,1])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 g,h
% extracted from script 'dPAW_salt.m'
% here we load in Data for 0-50 and 100-50 mM salt chemotaxis
% load in the fitted dPAW parameters runned with Chemotaxis_model_pop.m
% script and saved in the folder
% We compare kernels, tracks, and simulated CI

%% load fits and data
temp = load(fullfile(datadir,'data4plots', 'dPAW_salt50_100_0511.mat'));
x_100 = temp.x;
data100_fit = temp.Data_fit;
temp = load(fullfile(datadir,'data4plots', 'dPAW_salt0_50_0513_2.mat'));
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
temp = load(fullfile(datadir,'data4plots', 'Data_salt100_50_0511.mat'));
data_100_50 = temp.Data;
temp = load(fullfile(datadir,'data4plots', 'Data_salt0_50.mat'));
data_0_50 = temp.Data;

CI_data = zeros(1,2);
CI_data(1) = Data2CI(data0_fit);
CI_data(2) = Data2CI(data100_fit);
CI_data

%% sim tracks (fig. 3h)
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
reps = 10;

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

%% showcase tracks (fig. 3g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data
rng(18)
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

%% time plot
%%% plotting
rng(18)
nexp = 15;
temp = randperm(100);
selected = temp(1:nexp);

figure()
ax1 = axes;
imagesc(ax1,M0,'XData',[0 size(M0,2)*1],'YData',[0 size(M0,1)*1]);
set(gcf,'color','w'); set(gca,'Fontsize',20); %xticks([]); yticks([]);
colormap()
hold on
ax2 = axes;
for ii = 1:nexp
    temp = data_0_50(selected(ii)).xy;
%     temp = data_100_50(selected(ii)).xy;
    temp = temp(:,1:10:end);
    xx = temp(1,:);
    yy = temp(2,:);
    ll = length(temp);
    gg = linspace(0,1,ll);
%     plot(xx, yy, 'Color', [grayLevel grayLevel grayLevel]);
    patch(ax2, [xx nan]*1,[yy nan]*1,[gg nan],[gg nan], 'edgecolor', 'interp','LineWidth',2); 
    hold on
    plot(ax2,xx(1)*1, yy(1)*1,'g.', 'MarkerSize',25)
    plot(ax2,xx(end)*1, yy(end)*1,'r.', 'MarkerSize',25)
end
xlim(ax2, [0, size(M0,2)]);
ylim(ax2, [0, size(M0,1)]);
set(gca, 'YDir','reverse')
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
c = gray;
colormap(ax2,c)
colormap(ax1)

%% model
rng(4)
nexp = 15;
figure;
subplot(121)
imagesc(M0); colorbar(); caxis([0 100]); hold on
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
imagesc(M100); colorbar(); caxis([0 100]); hold on
specs.M = M100;
specs.REP = nexp;
[tracks, CI] = param2tracks(x_100, specs, []);
for ii = 1:nexp
    temp = tracks(ii).xy;
    plot(temp(1,:), temp(2,:),'k', 'LineWidth',1.5); hold on
    plot(temp(1,1), temp(2,1), 'g.','MarkerSize',13); plot(temp(1,end), temp(2,end), 'r.','MarkerSize',13);
end
set(gcf,'color','w'); set(gca,'Fontsize',20);

%% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
