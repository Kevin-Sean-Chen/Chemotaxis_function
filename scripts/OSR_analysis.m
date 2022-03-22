clear; clc;
%% OSR_analysis
%%%
% triggered behavior after periodic input
%%%
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker-new\Experimental Analysis')
fields_to_load = {'Path','Time','Behaviors','LEDVoltages','Frames'};
% fields_to_load = {'Path','Time','Behaviors','AngSpeed','ProjectedEigenValues'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);%(1:3)

allt = 0; for ii = 1:length(Tracks); temp = Tracks(ii).Time; allt = allt+ temp(end)-temp(1); end  %animal-hours

%% load LED
LED_txt = '/projects/LEIFER/Kevin/Data_odor_flow_equ/20220207_custom_stimulus_AML612_periodic_test/Data20220207_171935/LEDVoltages.txt';
LED_txt = '/projects/LEIFER/Kevin/Data_odor_flow_equ/20220212_custom_stimulus_AML612_periodic_test/Data20220212_180328/LEDVoltages.txt';
fileID = fopen(LED_txt, 'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

%% find transitions
absdiff = (diff(A));
local_peaks = find(absdiff>3.5);
local_peaks2 = find(diff(local_peaks)>500);
plot(local_peaks2,'-o')
stim_locs = local_peaks(local_peaks2);
plot(stim_locs,'-o')

%%% testing
stim_locs = stim_locs - 14*0;
%% looping
psth = [];
window = 14*10;
beh = 7;

for tt = 1:length(stim_locs)-1
    [sub_tracks, sub_id] = FilterTracksByTime(Tracks, stim_locs(tt)+1,stim_locs(tt)+window, false);  %Mochi's code to select tracks within a time window
    
    for ii = 1:length(sub_tracks)
        temp = sub_tracks(ii).Behaviors(beh,:);
        if sum(temp) > 0
            pos = find(temp>0);
            offset = stim_locs(tt) - sub_tracks(ii).Frames(1);
            if offset<0
                pos = pos + abs(offset) - 1;
            end
            temp2 = zeros(1,window);
            temp2(pos) = 1;
            psth = [psth; temp2];
        end
    end
end

%%
figure;
plot(mean(psth))

%% real PSTH!!
absdiff = (diff(A));
local_peaks = find(absdiff>3.5);
local_peaks2 = find(diff(local_peaks)>500);
% plot(local_peaks2,'-o')
stim_locs = local_peaks(local_peaks2);

%%% testing
stim_locs = stim_locs - 640;

psth = [];
window = 800;
beh = 7;

select_trials = [3 6 8 10 12 24];
% select_trials = [1 2 4 5 13 14 15 22 23 25];
% select_trials = [7 9 11 16:22];
select_trials = 1:26;
for ti = 1:length(select_trials)  %length(stim_locs)-1
    tt = select_trials(ti);
    [sub_tracks, sub_id] = FilterTracksByTime(Tracks, stim_locs(tt)+1,stim_locs(tt)+800, false);  %Mochi's code to select tracks within a time window
    
    for ii = 1:length(sub_tracks)
        temp = sub_tracks(ii).Behaviors(beh,:);
        if sum(temp) > 0
            pos = find(temp>0);
            offset = stim_locs(tt) - sub_tracks(ii).Frames(1);
            if offset<0
                pos = pos + abs(offset) - 1;
            end
            temp2 = zeros(1,window);
            temp2(pos) = 1;
            psth = [psth; temp2];
        end
    end
end
%%
figure;
plot(mean(psth))

%% separate trials
temp = zeros(1,84);
temp(1:24) = 1;
temp2 = repmat(temp,1,8);
cov_stim = conv(A,temp2,'same');
figure; plot(cov_stim)