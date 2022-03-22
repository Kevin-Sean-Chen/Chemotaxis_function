%%% BTA_flow

clear; clc;
%% File_loading
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
fields_to_load = {'Path','Time','Speed','Behaviors','LEDVoltages','Pirouettes','Runs'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);%(1:3)

%% params
beh = 8;  %state of interest
win = 70;
acs = 70;

%% BTA
alltrigs = [];
n_events = 0;
for ww = 1:length(Tracks)
    
    %%% for behavior modes
    bi = Tracks(ww).Behaviors(beh,:);
    pos = find(bi == 1);
    
    %%% for Pirouttes
    pos = Tracks(ww).Pirouettes;
    if isempty(pos)~=1
        pos = pos(:,1);
    end
    
    trigs = [];%zeros(length(pos),win+acs+1);
    for ii = 1:length(pos)
        if pos(ii)-win>0 && pos(ii)+acs<length(bi)
            trigs = [trigs ; Tracks(ww).LEDVoltages(pos(ii)-win:pos(ii)+acs)];
        end
    end
    n_events = n_events+length(pos);
    alltrigs = [alltrigs; trigs];
end

figure;
plot([-acs:win]*(1/14), mean(alltrigs))

%% conditioned BTA
OS_threshold = 14000;

beh = 1:3;  %state of interest
win = 100;
acs = 100;

alltrigs = [];
n_events = 0;
for ww = 1:length(Tracks)
    
    path_i = Tracks(ww).Path;
    finalc = Fcon(path_i(1,end),path_i(1,end));
    
    if finalc > OS_threshold
        
        %%% behavior mode
        bi = Tracks(ww).Behaviors(beh,:);
        pos = find(bi == 1);
        
        %%% for Pirouttes
        pos = Tracks(ww).Runs;%Pirouettes;
        if isempty(pos)~=1
            pos = pos(:,1);
        end
    
        trigs = [];%zeros(length(pos),win+acs+1);
        for ii = 1:length(pos)
            if pos(ii)-win>0 && pos(ii)+acs<length(bi)
                trigs = [trigs ; Tracks(ww).LEDVoltages(pos(ii)-win:pos(ii)+acs)];
            end
        end
        n_events = n_events+length(pos);
        alltrigs = [alltrigs; trigs];
    end
end

figure;
plot([-acs:win]*(1/14), mean(alltrigs))
n_events
