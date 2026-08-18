%% make example staPAW demo from data
%% LOAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load folder with indivisual worm posture and time series
% inspect the behavior
% inference with staPAW
% save as vidieo...

%% loading data files
data_path = '/projects/LEIFER/Kevin/Data_salt/20240126_GWN_N2_salt_50_0/Data20240126_195905';

tracks_file = fullfile(data_path, 'analysis');
img_file = fullfile(data_path, 'individual_worm_imgs');

%% load with fields
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothX','SmoothY', 'LEDVoltages'};%,'Behaviors'};
% folder_names = getfoldersGUI();  %%% if using GUI
Tracks = loadtracks(data_path,fields_to_load);

%% load fitted model
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat')
rng(42)
rr = 1;
staPAW_fit = all_record(rr,2,1).params;  % two-state model  %%% num_cv, num_states, num_repeats

%% make kinemtatics structure
% make environment
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
M = Cmap.vq1;
M = fliplr(flipud(M));  %flipped camera
[rows, cols] = size(M);
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
% [x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M = (y_*0+1) .* gradient_x;

% load path
path_file = fullfile(tracks_file, 'Path.mat');
path = load(path_file);
paths = path.values;
n_tracks = length(paths);
cand = 1:n_tracks;

% make Data structure
[Data] = track2data(Tracks, cand, M);

%% inspection
% figure()
% for ii = 1:length(Data)
%     
%     subplot(311)
%     plot(Data(ii).dc)
%     ylabel('dc')
%     
%     subplot(312)
%     plot(Data(ii).dth)
%     ylabel('dth')
%     
%     subplot(313)
%     plot(Data(ii).dis)
%     ylabel('dis')
%     
%     sgtitle(['Track ' num2str(ii)])
%     
%     pause()
% end

%% MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% note good candiates for demo
demo_id = [1,8,9,11,13,17,35];
track_id = 11;

file_name = sprintf('worm_%d.mat', track_id);
worm_data = load(fullfile(img_file, file_name));
worm_images = worm_data.worm_images;

%%% down sample to match
bin = 5;
l_window = 1;
N = length(Data(track_id).dth);
img_idx = 1 + l_window*bin + (0:N-1)*bin;
worm_images_ds = worm_images(:,:,img_idx);

%% staPAW inference
[xxf, yyf, alltrials, alltime] = data2xy(Data(track_id));
alldis = extractfield(Data(track_id), 'dis');
yyf = [yyf; alldis];

wind_test = [1100:1600];%[3000:3500];
xx = xxf;%(:,wind_test);
yy = yyf;%(:,wind_test);
maskf = alltrials;
mask = maskf;%(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(staPAW_fit,xx,yy,mask);

p_z = gams_(1,:);  %%% check if it is first or second channel
p_z = conv(p_z, ones(1,5), 'same')/5;

%% visualization
%%
%%% static
time_vec = [1:length(xx)]*(5/14);
figure()
subplot(411)
plot(time_vec, xxf(1,:))
subplot(412)
plot(time_vec, yyf(1,:))
subplot(413)
plot(time_vec, yyf(2,:))
subplot(414)
imagesc(p_z)

figure()
subplot(121)
imagesc(worm_images(:,:,1))
subplot(122)
plot(Data(track_id).xy(1,:), Data(track_id).xy(2,:))

%% VIDEO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% dynamic

time_vec = (1:length(xx)) * (5/14);
% threshold for state
p_thresh = 0.2;

%%% size and vido settings
video = VideoWriter('worm_dynamics_good_exp.avi', 'Motion JPEG AVI');
speed_up = 5;
video.FrameRate = 14/5 * speed_up;
video.Quality = 95;
open(video);
set(groot, 'defaultAxesFontSize', 16)
set(groot, 'defaultTextFontSize', 16)

figure('Units','normalized', 'Position',[0.02 0.05 0.95 0.85])

% 3 columns for visualization, 2 columns for dynamics
t = tiledlayout(4, 5, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% =========================================================
%  Upper left: worm image
% =========================================================
ax_img = nexttile(1, [2 3]);

h_img = imagesc(flipud(worm_images_ds(:,:,1)));
axis image
axis off
title('Worm')

hold on

% red surrounding box for state
[nrow, ncol, ~] = size(worm_images_ds);

h_statebox = rectangle( ...
    'Position',[0.5 0.5 ncol-1 nrow-1], ...
    'EdgeColor','r', ...
    'LineWidth',5, ...
    'Visible','off');


%% =========================================================
%  Lower left: trajectory
% =========================================================
ax_xy = nexttile(11, [2 3]);

pix2mm = 1/31.5;

xtraj = Data(track_id).xy(1,:);
ytraj = Data(track_id).xy(2,:);

xlims = [min(xtraj) max(xtraj)];
ylims = [min(ytraj) max(ytraj)];

% concentration landscape
h_bg = imagesc(ax_xy, [1 cols], [1 rows], M);
hold on

set(ax_xy, 'YDir', 'reverse')

% concentration colormap
colormap(ax_xy, parula(256))

% colorbar
cb = colorbar(ax_xy);
cb.Label.String = 'Concentration (mM)';

% trajectory
plot(xtraj, ytraj, 'k-', 'LineWidth', 1.5);

h_pos = plot(xtraj(1), ytraj(1), ...
             'ro', ...
             'MarkerFaceColor','r', ...
             'MarkerSize',8);

axis equal
xlim(xlims)
ylim(ylims)
title('Trajectory')

% convert displayed ticks from pixels -> mm
xt = get(ax_xy, 'XTick');
yt = get(ax_xy, 'YTick');

set(ax_xy, ...
    'XTickLabel', compose('%.1f', xt * pix2mm), ...
    'YTickLabel', compose('%.1f', yt * pix2mm))

xlabel('x (mm)')
ylabel('y (mm)')

%% =========================================================
%  Right: dynamics
% =========================================================

ax1 = nexttile(4, [1 2]);
plot(time_vec, xxf(1,:), 'k')
ylabel('concentration')
hold on
line1 = xline(time_vec(1), 'r--', 'LineWidth',1.5);

ax2 = nexttile(9, [1 2]);
plot(time_vec, yyf(1,:), 'k')
ylabel('heading change')
hold on
line2 = xline(time_vec(1), 'r--', 'LineWidth',1.5);

ax3 = nexttile(14, [1 2]);
plot(time_vec, yyf(2,:), 'k')
ylabel('speed')
hold on
line3 = xline(time_vec(1), 'r--', 'LineWidth',1.5);

%% =========================================================
%  Bottom right: state probability
% =========================================================
ax4 = nexttile(19, [1 2]);

imagesc( ...
    [time_vec(1) time_vec(end)], ...
    [0 1], ...
    p_z);

% white -> red colormap
ncolor = 256;
redmap = [ ...
    ones(ncolor,1), ...
    linspace(1,0,ncolor)', ...
    linspace(1,0,ncolor)'];

colormap(ax4, redmap)
caxis(ax4, [0 1])

yticks([])
ylabel('P(z)')
xlabel('Time (s)')
hold on

line4 = xline(time_vec(1), 'k--', 'LineWidth',1.5);

% align all dynamics panels
linkaxes([ax1 ax2 ax3 ax4], 'x')

%% =========================================================
%  March through time
% =========================================================

for ii = 1:length(time_vec)

    % update worm image
    h_img.CData = (worm_images_ds(:,:,ii));

    % update position on trajectory
    h_pos.XData = Data(track_id).xy(1,ii);
    h_pos.YData = Data(track_id).xy(2,ii);

    % current time
    current_time = time_vec(ii);

    % update time markers
    line1.Value = current_time;
    line2.Value = current_time;
    line3.Value = current_time;
    line4.Value = current_time;

    % -----------------------------------------------------
    % state indicator
    % -----------------------------------------------------
    if p_z(ii) > p_thresh
        h_statebox.Visible = 'on';
    else
        h_statebox.Visible = 'off';
    end

    sgtitle(sprintf( ...
        'Frame %d   |   t = %.2f s   |   p_z = %.2f', ...
        ii, current_time, p_z(ii)));

    drawnow

    %%% manually move forward
%     pause()
    %%% video
    drawnow
    frame = getframe(gcf);
    writeVideo(video, frame);
    
end

close(video);
