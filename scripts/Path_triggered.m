%%% Path_triggered
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
% fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothSpeed','AngSpeed','SmoothX','SmoothY','Behaviors', 'LEDVoltages'};
fields_to_load = {'Path','SmoothX','SmoothY', 'LEDVoltages'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% odor landscape
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');
Fcon = Fcon.F;
M = Cmap.vq1;
M = fliplr(flipud(M));

%% simple triggering analysis
poly_degree = 3;  %polynomial fitting for moving window
filt = 14;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
bin = 7;  %down-sampling
nn = length(Tracks); %number of worms selected
windt = 14 * (1/(bin*fr)); %time window in seconds
acst = 4 * (1/(bin*fr));  % acausal window
trig_angs = [];
kk = 0;
turn_thr = 100;  % threshold for turning
n_up = 0;  n_down = 0;  % record movement direction upon the pulse!

figure
for ii = 1:nn
    
    %%% calculate angles
    temp = Tracks(ii);  % for loading saved tracks
    temp1 = zeros(round(size(temp.Path,1)/1),2);
    temp1(:,1) = smooth(temp.SmoothX, filt,'sgolay',poly_degree); %smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp.SmoothY, filt,'sgolay',poly_degree); %smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    vecs = diff(subs);
    vecs = [vecs; vecs(end,:)];   % compensating for the size change/time shift after taking difference 
    angs = zeros(1,size(vecs,1));    
    
    %%% iterate through worms
    for dd = 2:length(angs)
        %%% angle function
        angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
    end
    
    %%% conditioning
    if temp1(end,1)>500 && temp1(end,1)<2500 && temp1(end,2)>500 && temp1(end,2)<2000 && ...
       M(floor(temp1(1,2)), floor(temp1(1,1))) - M(floor(temp1(end,2)), floor(temp1(end,1))) < 0. % thresholding tracks
%         M(floor(temp1(1,2)), floor(temp1(1,1))) - M(floor(temp1(end,2)), floor(temp1(end,1))) > 0. % thresholding tracks
        kk = kk+1;
    %%% LED signal
    stim = Tracks(ii).LEDVoltages;
    stim = stim(1:bin:end);
    pos = find(diff(stim)>1);% & diff(stim)>1);  %0.01 1 3 thresholding to find impulse
    pos = randi([1,length(stim)],1,length(pos)*1);%randi(length(stim));  %%%% important: randomized here  %%%
%     pos = find(diff(stim)<-1);
    turn_i = zeros(1,length(floor((1-acst):(1+windt))));
    if isempty(pos)~=1
        for pp = 1:length(pos)
            if (pos(pp)+windt)<length(stim) && pos(pp)-acst>1
                wind = floor((pos(pp)-acst):(pos(pp)+windt));
                tempa = angs(wind);
%                 turn_i(find(abs(temp)>turn_thr)) = 1;
                trig_angs = [trig_angs; tempa]; %turn_i]; %
             
%                 plot(subs(wind,1),subs(wind,2)); hold on

                %%% record turn direction upon impulse
                if M(floor(subs(wind(1)+10,2)), floor(subs(wind(1)+10,1))) - M(floor(subs(wind(end),2)), floor(subs(wind(end),1)))<0 && ...
                        max(abs(tempa(10:20)))>turn_thr
                    n_up = n_up+1;
                elseif max(abs(tempa(10:20)))>turn_thr
                    n_down = n_down+1;
%                 else 
%                     break
                end
                %%%
            end
        end
    end
    
%     plot(temp1(:,1), temp1(:,2)); hold on;
    end
end

figure;
% plot([-acst:windt]*((bin*fr)), mean(abs(trig_angs)) / (fr*bin))
t_vec = [-acst:windt]*((bin*fr));
mean_ang = mean(abs(trig_angs)) / (fr*bin);
std_ang = std(abs(trig_angs) / (fr*bin)) / sqrt(size(trig_angs,1));
plot(t_vec, mean_ang, 'k', 'LineWidth',3)
hold on
xArea = [t_vec, fliplr(t_vec)];
yArea = [mean_ang + std_ang, fliplr(mean_ang - std_ang)];
fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%%% CI upon impulse!
(n_up-n_down)/(n_up+n_down)

%% triggered plot
figure
title('aversive; up gradient')
patch([0 5 5 0], [20 20, 60 60], [0.7 0.7 0.9])
hold on
% plot([-acst:windt]*((bin*fr)), mean(abs(trig_angs)) / (fr*bin))
plot(t_vec, mean_ang, 'k', 'LineWidth',3)
hold on
xArea = [t_vec, fliplr(t_vec)];
yArea = [mean_ang + std_ang, fliplr(mean_ang - std_ang)];
fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('time (s)')
ylabel('angular speed (|deg|/s)')
set(gca,'FontSize',20); set(gcf,'color','w');
ylim([20,60])

%% bar analysis
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/app_trigs.mat')
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/nai_trigs.mat')
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/ave_trigs.mat')
all_ang_trigs = {app_up_trigs, app_down_trigs, nai_up_trigs, nai_down_trigs, ave_up_trigs, ave_down_trigs};
dang_means = zeros(1,6);
dang_stds = zeros(1,6);
temp_stats = cell(1,6);

figure
for ii = 1:6
    temp_trigs = all_ang_trigs{ii};
    temp = mean(abs(temp_trigs(:,acst:acst+14)),2) - mean(abs(temp_trigs(:,1:acst)),2);  % post - pre absolute average
    temp = temp/(fr*bin);   % per time unit
    dang_means(ii) = mean(temp);
    dang_stds(ii) = std(temp)/sqrt(length(temp));
    temp_stats{ii} = temp;
    bar(ii, dang_means(ii)); hold on
    errorbar(ii, dang_means(ii), dang_stds(ii), 'k.')
end

%% for joint analyis without context
joint_ang_trigs = {[app_up_trigs; app_down_trigs], [nai_up_trigs; nai_down_trigs], [ave_up_trigs; ave_down_trigs]};
jang_means = zeros(1,3);
jang_stds = zeros(1,3);
temp_stats = cell(1,3);

figure
for ii = 1:3
    temp_trigs = joint_ang_trigs{ii};
    temp = mean(abs(temp_trigs(:,acst:acst+14)),2) - mean(abs(temp_trigs(:,1:acst)),2);  % post - pre absolute average
    temp = temp/(fr*bin);   % per time unit
    temp_stats{ii} = temp;
    jang_means(ii) = mean(temp);
    jang_stds(ii) = std(temp)/sqrt(length(temp));
    
    bar(ii, jang_means(ii)); hold on
    errorbar(ii, jang_means(ii), jang_stds(ii), 'k.')
end

%% Raster stateanalysis
test = [app_down_trigs;app_up_trigs];
bin_raster = test;
bin_raster = test*0;bin_raster(abs(test)>50) = 1;
figure;imagesc(bin_raster)
figure;plot(mean(bin_raster,1))

%%
load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/state_fit/two_states.mat')
yy2 = abs(yy(1,:));
pos = find(yy2>100);
durs = diff(pos);
figure; hist(durs,100)
[aa,bb] = hist(durs,100);
figure;plot(bb*5/14,aa)

%%
figure
patch([0 5 5 0], [0.03 .03, .1 .1], [0.7 0.7 0.9])
hold on;
plot(t_vec,mean(bin_raster,1)')
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)')
ylabel('P(turn)')
ylim([0.03,0.1])
xlim([-4 14]);