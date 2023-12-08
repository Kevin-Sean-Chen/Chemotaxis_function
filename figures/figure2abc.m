% Figure 2a,b,c
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2a schematic
% extracted from script 'exp_tracks.m'

%% load data
% load example chemotaxis trajectory
load(fullfile(datadir,'data4plots', 'Data_exp_dcdthdxy.mat'));
Cmap = load(fullfile(datadir,'data4plots', 'Landscape_low_0623_2.mat'));
M = Cmap.vq1;
M = fliplr(flipud(M));

%% plotting the example
% track ID and time window
id = 5;
wind = [650:1780];

% normalize time with color
n = length(wind);
temp = [n:-1:1]*230/n;
ccd = uint8(zeros(4,n));
for i=1:n
    ccd(:,i) = uint8(temp(i));
end
time_x = [1:n]*5/14;
pix2mm = 1/31.5;

figure;
ax1 = axes;
imagesc(ax1,M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
colormap()
hold on
ax2 = axes; 
xx = Data(id).xy(1,wind);
yy = Data(id).xy(2,wind);
ll = length(xx);
gg = fliplr(linspace(0,1,ll));%linspace(0,1,ll);
patch(ax2, [xx nan]*pix2mm,[yy nan]*pix2mm,[gg nan],[gg nan], 'edgecolor', 'interp','LineWidth', 5); hold on
    
plot(ax2,Data(id).xy(1,wind(1))*pix2mm, Data(id).xy(2,wind(1))*pix2mm, 'g.','MarkerSize',40)
plot(ax2,Data(id).xy(1,wind(end))*pix2mm, Data(id).xy(2,wind(end))*pix2mm, 'r.','MarkerSize',40)
set(gca, 'YDir','reverse')
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
c = gray;
colormap(ax2,c)
colormap(ax1)

%% time series extracted
figure;
subplot(311)
p = plot(time_x,Data(id).dcp(wind), 'LineWidth',2)
ylabel('dC^{\perp}')
subplot(312)
p = plot(time_x,Data(id).dc(wind), 'LineWidth',2)
ylabel('dC')
subplot(313)
p = plot(time_x,Data(id).dth(wind), 'LineWidth',2)
hold on
ylabel('d\theta')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2b,c kernels
% computed with scripts "Hess4MLE.m" and "compareK.m"

%% load data
load(fullfile(datadir,'data4plots', 'kernels_mle_std.mat'));  % loads mle_params and MLE_std arrays

%% some post analysis for variability!
% for Kc kernel
ttl = {'appetitive','naive','aversive'};
col = {'b','k','r'};
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
K = size(mle_params,1);
figure
for cc = 1:3
%     subplot(1,3,cc);
    mlee = squeeze(median(mle_params(:,cc,:)));
    y = mlee(3:6)'*cosBasis';
    mle_hess = MLE_std(cc, 3:6);
    standardError = mle_hess*cosBasis';
    plot(tt,y,col{cc},'LineWidth',3)
    hold on
    
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y + standardError, fliplr(y - standardError)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    title(ttl{cc})
end

%% 
% for Kc_perp kernel
ttl = {'appetitive','naive','aversive'};
col = {'b','k','r'};
tt = [1:length(cosBasis)]*5/14;
figure
for cc = 1:3
%     subplot(1,3,cc);
    mlee = squeeze(median(mle_params(:,cc,:)));
    y = -mlee(8).*exp(-tt./mlee(9));
    mle_hess = MLE_std(cc,8:9)/sqrt(10);  % compute sem here
    standardError = -mle_hess(1).*exp(-tt./mle_hess(2));
    plot(tt,y,col{cc},'LineWidth',3)
    hold on
    
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y + standardError, fliplr(y - standardError)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    title(ttl{cc})
end
