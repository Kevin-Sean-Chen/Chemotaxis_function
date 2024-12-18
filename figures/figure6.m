% Figure 6
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6b,c
% extracted from script 'Path_triggered.m'

%% load data
load(fullfile(datadir,'data4plots', 'app_trigs.mat'));
load(fullfile(datadir,'data4plots', 'nai_trigs.mat'));
load(fullfile(datadir,'data4plots', 'ave_trigs.mat'));

%% analysis parameters
poly_degree = 3;  %polynomial fitting for moving window
filt = 14;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
bin = 7;  %down-sampling
windt = 14 * (1/(bin*fr)); %time window in seconds
acst = 4 * (1/(bin*fr));  % acausal window

%% bar plot
all_ang_trigs = {app_up_trigs, app_down_trigs, nai_up_trigs, nai_down_trigs, ave_up_trigs, ave_down_trigs};
dang_means = zeros(1,6);
dang_stds = zeros(1,6);

figure
for ii = 1:6
    temp_trigs = all_ang_trigs{ii};
    temp = mean(abs(temp_trigs(:,acst:acst+14)),2) - mean(abs(temp_trigs(:,1:acst)),2);  % post - pre absolute average
    temp = temp/(fr*bin);   % per time unit
    dang_means(ii) = mean(temp);
    dang_stds(ii) = std(temp)/sqrt(length(temp));
    
    bar(ii, dang_means(ii)); hold on
    errorbar(ii, dang_means(ii), dang_stds(ii), 'k.')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6d
% extracted from script 'Hess4MLE_opto.m'

%% load analyzed kernels
load(fullfile(datadir,'data4plots', 'opto_K_mle_std.mat'));

%% plotting Ks
ttl = {'appetitive','naive','aversive'};
col = {'b','k','r'};
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
K = size(mle_params_opto,1);
figure
for cc = 1:3
    subplot(1,3,cc);
    mlee = squeeze(nanmedian(mle_params_opto(cc,:,:),3));  %((mle_params_opto(cc,:,5))); %
    y_odor = mlee(3:6)*cosBasis';
    y_opto = mlee(7:10)*cosBasis';
    mle_hess_odor = MLE_std_opto(cc,3:6)*1;%4/sqrt(length(Data));   % odor
    mle_hess_opto = MLE_std_opto(cc,7:10)*1;  % opto
    
    standardError_odor = mle_hess_odor*cosBasis';
    yyaxis left; 
    plot(tt,y_odor,col{cc},'LineWidth',3); ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
    ylabel('K_{c}');
    hold on
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y_odor + standardError_odor, fliplr(y_odor - standardError_odor)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    hold off
    
    yyaxis right
    standardError_opto = mle_hess_opto*cosBasis';
    plot(tt,y_opto,'Color',col{cc},'LineWidth',1); yliml = get(gca,'Ylim');% ,'Alpha', 0.3)
    ylabel('K_{opto}');
    hold on
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y_opto + standardError_opto, fliplr(y_opto - standardError_opto)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold off
    
    if yliml(2)*ratio<yliml(1)
        set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
    else
        set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
    end

    xlabel('time (s)');
    set(gca,'FontSize',20); set(gcf,'color','w'); title(ttl{cc})
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5e
% extracted from script 'Path_triggered.m'

%% load data
% same triggered traces as 5b
all_ang_trigs = {app_up_trigs, app_down_trigs, nai_up_trigs, nai_down_trigs, ave_up_trigs, ave_down_trigs};
trig_plot = reshape(all_ang_trigs,2,3)';

titles = {'appetitive', 'naive', 'aversive'};
t_vec = [-acst:windt]*((bin*fr));
line_symb = {'-','--'};

%% plottomg
for ii = 1:3
    figure;
    title(titles{ii})
    patch([0 5 5 0], [20 20, 60 60], [0.7 0.7 0.9])
    hold on
    for jj = 1:2
        % load data and compute mean and sem
        trig_angs = trig_plot{ii,jj};
        mean_ang = mean(abs(trig_angs)) / (fr*bin);
        std_ang = std(abs(trig_angs) / (fr*bin)) / sqrt(size(trig_angs,1));
        % plotting
        plot(t_vec, mean_ang, 'k','LineStyle',line_symb{jj}, 'LineWidth',3)
        hold on
        xArea = [t_vec, fliplr(t_vec)];
        yArea = [mean_ang + std_ang, fliplr(mean_ang - std_ang)];
        fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xlabel('time (s)')
        ylabel('angular speed (|deg|/s)')
        set(gca,'FontSize',20); set(gcf,'color','w');
        ylim([20,60])
    end
end
