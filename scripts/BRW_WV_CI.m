%% BRW_WV_CI
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% load odor landscape
test = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/Landscape.mat');
Cmap =load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/Landscape.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/OdorFx.mat');

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20220113_GWN_app+_MEK110mM_gasphase_30ml_200air/Landscape.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20220113_GWN_app+_MEK110mM_gasphase_30ml_200air/OdorFx.mat');

% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');

Fcon = Fcon.F;

%% test with RBW and WV index vs. CI (biased-random walk, weathervaning, and chemotaxis index)
% chemotaxis = struct('brw','wv','ci');
dC_window = 14*5;  %time window for dC measurement for turns
time_wind = 60*15;  %first few minutes

% initializing counts
run_up = 0;  %recording for runs
run_dn = 0;
turn_up = 0;  %recording for turns
turn_dn = 0;
ci_low = 0;  %recording for chemotaxis performance
ci_high = 0;

for c = 1:length(Tracks)
    %%% load observations
    prs = Tracks(c).Pirouettes;
    paths = Tracks(c).Path;
    runs = Tracks(c).Runs;
    times = Tracks(c).Time;
    
    %%%%% conditioning tracks (in space or in time)
    if max(times)<=time_wind
    
    %%% computing runs (biased-random walk)
    if isempty(runs)==0
    for rr = 1:size(runs,1)
        path_i = paths(runs(rr,1):runs(rr,2),:);  %segment of running
        dC = Fcon(path_i(1,1) , path_i(1,2)) - Fcon(path_i(end,1) , path_i(end,2));  %delta C
        if dC>=0  %up gradient  %%%%%%%%% hack for now to invert concentration~~~
            run_up = [run_up length(path_i)*dC];  %record run length
        elseif dC<0
            run_dn = [run_dn length(path_i)*abs(dC)]; 
        end      
    end
    end
    
    %%% computing turns (weathervaning)
    if isempty(prs)==0
    for tt = 1:size(prs,1)
        c_i = max((prs(tt,1)-dC_window), 1);  %initiation time 
        c_f = min((prs(tt,2)+dC_window), length(paths));  %exiting time
        dC = Fcon(paths(c_i,1) , paths(c_i,2)) - Fcon(paths(c_f,1) , paths(c_f,2));  %delta C
        if dC>=0  %up gradient
            turn_up = turn_up + dC;  %record run length
        elseif dC<0
            turn_dn = turn_dn + abs(dC); 
        end   
    end
    end
    
    %%% quantify chemotaxis index
    dC = Fcon(paths(1,1) , paths(1,2)) - Fcon(paths(end,1) , paths(end,2)); 
    if dC>=0
        ci_high = ci_high+dC;  %recording for chemotaxis performance
    elseif dC<0
        ci_low = ci_low+abs(dC);
    end
    
%     c
    end
    
end

% %% analysis
brw_index = (median(run_up)-median(run_dn)) / (median(run_up)+median(run_dn));
wv_index = ((turn_up)-(turn_dn)) / ((turn_up)+(turn_dn));
ci_ = (ci_high - ci_low) / (ci_high+ci_low);

disp(folder_names);
disp(['BRW: ',num2str(brw_index)]);
disp(['WV: ',num2str(wv_index)]);
disp(['CI: ',num2str(ci_)]);

%% groups
app_ = [0.09  0.21  0.35 ;
        0.15  0.16  0.34;
        0.09  0.16  0.32]';
nai_ = [0.08  0.09  0.15;
        0.13  0.07  0.26;
        0.13  0.069  0.21;
        0.09  0.1  0.25]';
ave_ = [0.04  0.14  0.19;
        0.01  0.18  0.25;
        0.03  0.21  0.25]';

    
%%% low C
app_ = [0.12, 0.16, 0.39;
        0.08, 0.15, 0.32;
        0.14  0.11,  0.38]';
nai_ = [0.08,0.07,0.26;
       0.10  0.071  0.20;
       0.08, 0.065, 0.21]';
%         0.09, 0.13, 0.24]';
ave_ = [0.05,-0.0008,0.08;
        0.02 , 0.02, 0.1;
        0.06, 0.0, 0.04;
        0.07, 0.04, 0.23]';

figure();
plot(app_,'b-o','linewidth',2); hold on
plot(nai_,'k-o','linewidth',2); hold on
plot(ave_,'r-o','linewidth',2); hold on
% plot([0.083  0.08  0.21],'k-*','linewidth',2)
plot([0.06  0.1  0.22],'r-*','linewidth',2)
names = {'BRW'; 'WV'; 'CI'};
ylabel('Index')
set(gca,'xtick',[1:3],'xticklabel',names,'FontSize',20)
set(gcf,'color','w');

%% bar plot
figure
m_ap = mean(app_');
m_na = mean(nai_');
m_av = mean(ave_');
s_ap = std(app_');
s_na = std(nai_');
s_av = std(ave_');

hBar = bar([m_ap;m_na;m_av]');
for k1 = 1:3
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');     
    ydt(k1,:) = hBar(k1).YData;                    
end
% set(hBar, {'DisplayName'}, {'App','Naive','Ave'}')
hold on
errorbar(ctr, ydt, [s_ap;s_na;s_av]', '.r')                  
hold off
ylabel('Index')
% set(gca,'linewidth',2,'FontSize',20)
set(gca,'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
legend([hBar(1), hBar(2),hBar(3)], 'Appetitive','Naive','Aversive')