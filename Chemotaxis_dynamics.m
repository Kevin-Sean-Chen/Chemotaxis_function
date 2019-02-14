%Chemotaxis_dynamics
%021419
clear
clc
%%
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%load data from one folder
% Tracks = [];
% load('Path.mat')
% paths = values;
% load('Time.mat')
% times = values;
% for ii = 1:length(paths); Tracks(ii).Path = paths{ii}; end
% for ii = 1:length(pat0hs); Tracks(ii).Time = times{ii}; end

%or batch analysis
fields_to_load = {'Path','Time','Frames','SmoothSpeed'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%%% criteria %%%
nn = length(Tracks); %number of worms selected
mint = 60;%60*1; %minimum time in seconds
disth = 500;  %radius of pixels from target
target = [2517,975];%[950,1100];%  %position of target/sourse of odorant (approximated from images)

%visualize all paths and check criteria
cand = [];
figure;
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        plot(Tracks(i).Path(:,1),Tracks(i).Path(:,2));% pause()
        hold on
        cand = [cand i];
    end
end

%% First-passage time measurement
crossing = 1800;
timing = [];
figure;
for i = 1:nn
    pos = find(Tracks(i).Path(:,1)>crossing);
    if size(pos) ~= 0
        timing = [timing Tracks(i).Time(pos(1))];
    end
end
%% Smoothed speed
speed = [];
for i = 1:nn
    speed = [speed Tracks(i).SmoothSpeed];
end

%% approximate diffusion coefficient thrugh time
partt = 5;
Mt = 30*60;
bin_times = [0:Mt/(partt):Mt];

Ds = [];
for tt = 1:length(bin_times)-1
    sub_tracks = FilterTracksByTime(Tracks,bin_times(tt)*14,bin_times(tt+1)*14);  %Mochi's code to select tracks within a time window
    xsqu = 0;
    for ii = 1:length(sub_tracks)
        xsqu = xsqu + mean(diff(sub_tracks(ii).Path(:,1)).^2);
    end
    Ds = [Ds xsqu/length(sub_tracks)];
end

%% approximate diffusion coefficient thrugh space
partt = 6;
Mt = 30*60;
%frame is 1944 X 2592
bin_space = [2592/2: (2592/2)/partt: 2592];

Ds = [];
for xx = 1:length(bin_space)-1
    xsqu = 0;
    tr = 0;
    for ii = 1:length(Tracks)
        pos = find(Tracks(ii).Path(:,1)>bin_space(xx) & Tracks(ii).Path(:,1)<bin_space(xx+1));
        if length(pos)>1
            subpath = Tracks(ii).Path(pos,1);
            xsqu = xsqu + mean(diff(subpath).^2);
            tr = tr+1;
        end
    end
    Ds = [Ds xsqu/tr];
end

