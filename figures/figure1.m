% Figure 1c,d
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1c
% extracted from script "exp_tracks.m"

%% load data
% load example tracks (comment out others to plot the given condition)
load(fullfile(datadir,'data4plots', 'exp_track_app.mat'));
% load(fullfile(datadir,'data4plots', 'exp_track_nai.mat'))
% load(fullfile(datadir,'data4plots', 'exp_track_ave.mat'))

% load odor map
Cmap = load(fullfile(datadir,'data4plots', 'Landscape_low_0623_2.mat'));
M = Cmap.vq1;
M = fliplr(flipud(M));

%% plotting
pix2mm = 1/31.5;
figure()
ax1 = axes;
imagesc(ax1,M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
colormap()
hold on
ax2 = axes;
for ii = 1:length(exp_track)
    xx = exp_track(ii).x';
    yy = exp_track(ii).y';
    ll = length(exp_track(ii).x);
    gg = linspace(0,1,ll);
    patch(ax2, [xx nan]*pix2mm,[yy nan]*pix2mm,[gg nan],[gg nan], 'edgecolor', 'interp','LineWidth',2); 
    hold on
    plot(ax2,xx(1)*pix2mm, yy(1)*pix2mm,'g.', 'MarkerSize',25)
    plot(ax2,xx(end)*pix2mm, yy(end)*pix2mm,'r.', 'MarkerSize',25)
end
set(gca, 'YDir','reverse')
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
c = gray;
colormap(ax2,c)
colormap(ax1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure d
% extracted with script "BRW_WV_CI.m"

%% load analyzed data
load(fullfile(datadir,'data4plots', 'BWCs2.mat'))

%% plotting bars
app_ = BWCs{1}'; nai_ = BWCs{2}'; ave_ = BWCs{3}';  % loading the summary statistics
figure
m_ap = mean(app_');
m_na = mean(nai_');
m_av = mean(ave_');
s_ap = std(app_')/sqrt(size(app_,2));
s_na = std(nai_')/sqrt(size(nai_,2));
s_av = std(ave_')/sqrt(size(ave_,2));

hBar = bar([m_ap;m_na;m_av]');
data_cell = {app_, nai_, ave_};
for k1 = 1:3
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');     
    ydt(k1,:) = hBar(k1).YData;  hold on
    plot(ctr(k1,:), data_cell{k1},'ko')
end
hold on 
errorbar(ctr, ydt, [m_ap;m_na;m_av]*0,[s_ap;s_na;s_av], '.k')         
hold off
ylabel('Index')
names = {'CI'; 'BRW'; 'WV'};
set(gca,'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
legend([hBar(1), hBar(2),hBar(3)], 'Appetitive','Naive','Aversive')

%% alternative way to compute from the large cell with all tracks
%% load analyzed data
load(fullfile(datadir,'data4plots', 'track_learn_all.mat'))

%%
cols = ['b','k','r'];
figure
for cc = 1:3  % learning conditions
    sub_learn_tracks = track_learn_all{cc};
    BWC = zeros(length(sub_learn_tracks), 3);
    for ff = 1:length(sub_learn_tracks)  % load data for each condition
        Tracks = sub_learn_tracks{ff};
        [ci_, brw_index, wv_index] = compute_index(Tracks, M, 30); %20
        BWC(ff,:) = [ci_, brw_index, wv_index];
    end
    plot(BWC',cols(cc)); hold on
    BWCs{cc} = BWC;  % can now use this cell to plot the same bar plot above
end
