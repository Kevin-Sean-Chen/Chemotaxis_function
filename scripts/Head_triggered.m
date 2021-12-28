clear; clc;
%% Head_triggered
%%%
%tragger analysis directly on heading, eigen-worms, and path directions
%%%
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker-new\Experimental Analysis')
% fields_to_load = {'Path','Time','Behaviors','LEDPower','Centerlines','AngSpeed','ProjectedEigenValues'};
fields_to_load = {'Path','Time','Behaviors','AngSpeed','ProjectedEigenValues'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);%(1:3)

allt = 0; for ii = 1:length(Tracks); temp = Tracks(ii).Time; allt = allt+ temp(end)-temp(1); end  %animal-hours

%% visualization
for w = 1:length(Tracks)
    temp = Tracks(w).Centerlines;
    for t = 1:size(temp,3)
        plot(squeeze(temp(:,1,t)),squeeze(temp(:,2,t)),'-o'); pause();
    end
end
%%  head_angles
nn = length(Tracks);
tip = 1;
base = 5;  %head is from base to tip (3 or 5)
endp = 20;  %ending for body-center reference
tip = 20;
base = 15;  %head is from base to tip (3 or 5)
endp = 1;  %ending for body-center reference
for w = 1:nn
    temp = Tracks(w).Centerlines;
    stim = Tracks(w).LEDPower;
    
    head_angs = zeros(1,size(temp,3));
    for t = 1:size(temp,3)
        heada = squeeze(temp(tip,:,t)) - squeeze(temp(base,:,t));
        centl = squeeze(temp(base,:,t)) - squeeze(temp(endp,:,t));
        head_angs(t) = angles(heada,centl);
    end
    
    %plot(head_angs); pause();
    plot(diff(head_angs)); pause();
    %autocorr(head_angs,100); pause();
    %spectrogram(head_angs,'yaxis');  pause();
%     [r,lags] = xcorr(head_angs,stim,'unbiased');
%     mid = floor(length(r)/2);  %zero time lag
%     win = 100;  %looking back and forth
%     plot(lags(mid-win:mid+win),r(mid-win:mid+win)); pause();
    
end

%% triggered head angle
nn = length(Tracks);
tip = 1;
base = 5;  %head is from base to tip
endp = 20;  %ending for body-center reference

% tip = 20;
% base = 15;  %head is from base to tip (3 or 5)
% endp = 1;  %ending for body-center reference

thr = 60;  %quite arbitrary for sudden head angle change (30 or 50 (typical head swings <30 according to data))
win = 150;
trigs = [];  %triggering stimuli
headyn = [];  %looking at the head dynamics
for w = 1:5000%nn
    temp = Tracks(w).Centerlines;
    stim = Tracks(w).LEDPower;
    
    head_angs = zeros(1,size(temp,3));
    for t = 1:size(temp,3)
        heada = squeeze(temp(tip,:,t)) - squeeze(temp(base,:,t));
        centl = squeeze(temp(base,:,t)) - squeeze(temp(endp,:,t));
        head_angs(t) = angles(heada,centl);
    end
    
    pos = find(diff(head_angs)<thr & diff(head_angs)>30);  %thresholding conditiond!!!
    for ii = 1:length(pos)
        if pos(ii)-win>0 && pos(ii)+win<length(head_angs)
            trigs = [trigs ; stim(pos(ii)-win:pos(ii)+win)];
            headyn = [headyn; head_angs(pos(ii)-50:pos(ii)+50)];
        end
    end
    
end

%% triggered eigen-values!!
moden = 2;  %the eigen-mode of interest
nn = length(Tracks);
thre = 1;  %arbitrary... setting a trigger point larger than 2 std of the time-series

win = 150;  %window around the triggering point
trigs = [];  %triggering stimuli
eigdyn = [];  %looking at the eigen-value dynamics
for w = 1:nn
    temp = Tracks(w).ProjectedEigenValues;
    temp = temp(moden,:);
    stim = Tracks(w).LEDPower;
    
    thre_std = zscore(temp);
    
    pos = find((thre_std)>thre);  %thresholding conditiond!!!
    for ii = 1:length(pos)
        if pos(ii)-win>0 && pos(ii)+win<length(thre_std)
            trigs = [trigs ; stim(pos(ii)-win:pos(ii)+win)];
            eigdyn = [eigdyn; thre_std(pos(ii)-50:pos(ii)+50)];
        end
    end
    
end
%% %%%%%%%%%%%%%%%%%%%%%%%% path-triggered analysis
%% visualize paths
smooth_n = 10;
for w = 1:100%length(Tracks)
%     w = I(w);
    w = randi(length(Tracks));
    temp = Tracks(w).Path;
    temp(:,1) = smooth(temp(:,1),smooth_n);  temp(:,2) = smooth(temp(:,2),smooth_n);  %%% remove tracking noise
    Pangs = zeros(1,size(temp,1)-1);
    for t = 2:size(temp,1)-1  %no angles for the first and last data point!!
        v1 = temp(t,:) - temp(t-1,:);
        v2 = temp(t+1,:) - temp(t,:);
        Pangs(t) = angles(v1,v2);
    end
    
    %%% thresholding conditiond!!! 
%     pos = find(Pangs<180 & Pangs>120);  %(by value)
    [pks,pos] = findpeaks(Pangs,'MinPeakDistance',10,'MinPeakHeight',60);  %(by local peak)
    pos2 = find(pks>120);
    pos(pos2) = [];  pks(pos2) = [];
%     plot(Pangs); hold on; plot(pos,pks,'*');  hold off;
    plot(temp(:,1),temp(:,2),'b-'); hold on;  plot(temp(pos,1),temp(pos,2),'r*');  hold off;
%     pause();
    hold on
end

%% path_angles
smooth_n = 10;
all_ang = [];
for w = 1:length(Tracks)
    temp = Tracks(w).Path;
    temp(:,1) = smooth(temp(:,1),smooth_n);  temp(:,2) = smooth(temp(:,2),smooth_n);
    Pangs = zeros(1,size(temp,1)-1);
    for t = 2:size(temp,1)-1  %no angles for the first and last data point!!
        v1 = temp(t,:) - temp(t-1,:);
        v2 = temp(t+1,:) - temp(t,:);
        Pangs(t) = angles(v1,v2);
    end
    %plot(Pangs); pause();
    if isreal(Pangs)==1
        all_ang = [all_ang sum(Pangs)];
    end
end

%% Path-Triggered analysis
[B,I] = sort(all_ang);  %some tracks are way too noisy

thr = 30;  %quite arbitrary for sudden head angle change (30 or 50 (typical head swings <30 according to data))
win = 150;
trigs = [];  %triggering stimuli
pathdyn = [];  %looking at the head dynamics
for w = 1:5000%%%nn
    w = I(w);  %starting with smoother paths
    stim = Tracks(w).LEDPower;
    
    temp = Tracks(w).Path;
    temp(:,1) = smooth(temp(:,1),smooth_n);  temp(:,2) = smooth(temp(:,2),smooth_n);
    Pangs = zeros(1,size(temp,1)-1);
    for t = 2:size(temp,1)-1  %no angles for the first and last data point!!
        v1 = temp(t,:) - temp(t-1,:);
        v2 = temp(t+1,:) - temp(t,:);
        Pangs(t) = angles(v1,v2);
    end
    
%     pos = find(Pangs<180 & Pangs>100);%find(Pangs>150);%  %thresholding conditiond!!!
    [pks, pos] = findpeaks(Pangs,'MinPeakDistance',10,'MinPeakHeight',120);  %(by local peak)
    pos2 = find(pks>180);
    pos(pos2) = [];  pks(pos2) = [];
    for ii = 1:length(pos)
        if pos(ii)-win>0 && pos(ii)+win<length(Pangs)
            trigs = [trigs ; stim(pos(ii)-win:pos(ii)+win)];
            pathdyn = [pathdyn; Pangs(pos(ii)-50:pos(ii)+50)];
        end
    end
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%% Angular speed triggered analysis
%%%using Angspeed variable from Mochi's behavioral analysis code
%%%adding a check for the speed of angle of change
%% visualize ang_speed (compare with angle traces)
smooth_n = 10;
for w = 1:length(Tracks)
    temp = Tracks(w).AngSpeed;
    agspd = smooth(temp,1);  %%% remove tracking noise
    plot(agspd); hold on;
    temp = Tracks(w).Path;
    temp(:,1) = smooth(temp(:,1),smooth_n);  temp(:,2) = smooth(temp(:,2),smooth_n);
    Pangs = zeros(1,size(temp,1));
    for t = 2:size(temp,1)-1  %no angles for the first and last data point!!
        v1 = temp(t,:) - temp(t-1,:);
        v2 = temp(t+1,:) - temp(t,:);
        Pangs(t) = angles(v1,v2);
    end
    plot(Pangs)
    pause();  hold off;
end

%% building phase portrait
smooth_n = 10;
phase = [];
for w = 1:length(Tracks)
    temp = Tracks(w).AngSpeed;
    agspd = smooth(temp,1);  %%% remove tracking noise
    
    temp = Tracks(w).Path;
    temp(:,1) = smooth(temp(:,1),smooth_n);  temp(:,2) = smooth(temp(:,2),smooth_n);
    Pangs = zeros(1,size(temp,1));
    for t = 2:size(temp,1)-1  %no angles for the first and last data point!!
        v1 = temp(t,:) - temp(t-1,:);
        v2 = temp(t+1,:) - temp(t,:);
        Pangs(t) = angles(v1,v2);
    end
    
    lump = [Pangs; agspd'];
    %pos = find(lump(1,:)<30);
    %lump(:,pos) = [];  %arbitrary cutoff
    pos = find(abs(lump(2,:))<30 | abs(lump(2,:))>360);
    lump(:,pos) = [];  %arbitrary cutoff
    phase = [phase lump];
end

%% path_angles
smooth_n = 10;
all_ang = [];
for w = 1:length(Tracks)
    temp = smooth( Tracks(w).AngSpeed,smooth_n);  %%% remove tracking noise
    all_ang = [all_ang sum(temp)];
end
plot(sort(abs(all_ang)))

%% Ang-speed triggered average
[B,I] = sort(all_ang);  %some tracks are way too noisy

thr = 30;  %quite arbitrary for sudden head angle change (30 or 50 (typical head swings <30 according to data))
win = 250;
trigs = [];  %triggering stimuli
pathdyn = [];  %looking at the head dynamics
for w = 1:nn%8000%%%
    w = I(w);  %starting with smoother paths
    
    agspd = smooth( Tracks(w).AngSpeed,smooth_n)';  %%% angle change speed
    agspd = abs(agspd);  %turning direction does not matter for now
    stim = Tracks(w).LEDPower;   %%% LED stimuli time series
    temp = Tracks(w).Path;   %%% path of center-of-mass
    temp(:,1) = smooth(temp(:,1),smooth_n);  temp(:,2) = smooth(temp(:,2),smooth_n);
    Pangs = zeros(1,size(temp,1));
    for t = 2:size(temp,1)-1  %no angles for the first and last data point!!
        v1 = temp(t,:) - temp(t-1,:);
        v2 = temp(t+1,:) - temp(t,:);
        Pangs(t) = angles(v1,v2);
    end
    
    pos = find(Pangs<120 & Pangs>60 & abs(agspd) > 150);%find(Pangs>150);%  %thresholding conditiond!!!
    for ii = 1:length(pos)
        if pos(ii)-win>0 && pos(ii)+win<length(Pangs)
            trigs = [trigs ; stim(pos(ii)-win:pos(ii)+win)];
            pathdyn = [pathdyn; Pangs(pos(ii)-50:pos(ii)+50)];
        end
    end
    
end


