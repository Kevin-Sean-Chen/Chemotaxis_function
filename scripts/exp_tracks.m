%%% exp_tracks
% good example tracks for three learning conditions
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothSpeed','AngSpeed','SmoothX','SmoothY'};%,'Behaviors'};

main_folder=('/projects/LEIFER/Kevin/Data_learn');
cd(main_folder)

folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% odor pre-processing visualization
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');

Fcon = Fcon.F;
M = Cmap.vq1;
M = fliplr(flipud(M));  %flipped camera

%% conditioning on track criteria
poly_degree = 3;  %polynomial fitting for moving window
filt = 7;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 60*10;% 5 10 %minimum time in seconds
minx = 450;% 450 250 500;  %minimum displacement (in terms of pixels)
endingt = 60*30;  %only taking the first few minutes
pix2mm = 1/31.5;
targ_track = 7;
conc_thre = 30;%0.8*max(max(M));

figure();
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
hold on
cand = [];  %index of tracks as candidates for analysis
alldists = [];
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); %pixel displacement
        alldists = [alldists displace*pix2mm^2];  %all dispacements in mm
        if displace > minx^2  %space cutoff
%             pos = find(Tracks(i).Path(1,2)<1000 | Tracks(i).Path(1,2)>1500);
            pos = find(Tracks(i).Time<endingt);  %time window cutoff (the later time points are less correct...)
            if isempty(pos)~=1
                x_smooth = smooth(Tracks(i).Path(:,1), filt,'sgolay',poly_degree);
                y_smooth = smooth(Tracks(i).Path(:,2), filt,'sgolay',poly_degree);
                    if M(floor(y_smooth(1)), floor(x_smooth(1))) < conc_thre
                        plot(x_smooth*pix2mm, y_smooth*pix2mm,'k','LineWidth',1); hold on;
                        plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.', 'MarkerSize',15)
                        plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.', 'MarkerSize',15)
                        cand = [cand i];
                    end
            end
        end
    end
end

xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');
% set ( gca, 'xdir', 'reverse' )

%% example tracks across conditions
cand_app = [591  593  616  642  3509 3521   3645   3649   3662  3704 3849 3855  3871];  % 10t 500x
cand_nai = [155  264  358   465   2389 2394 2397  2447 2482 2543   2672 2698  2708  2788 2789  2834 ];  % 5t 250x
cand_ave = [1225 1227  1437 1494 1596 1954 2063 2102 2206 2234 2864  2918  2968 3052 3101 3119];  % 10t 450x
cand = cand_app;
figure();
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
hold on
for j = 1:length(cand)
    i = cand(j);  %using the pre-selected candidate
    x_smooth = smooth(Tracks(i).Path(:,1), filt,'sgolay',poly_degree);
    y_smooth = smooth(Tracks(i).Path(:,2), filt,'sgolay',poly_degree);
    plot(x_smooth*pix2mm, y_smooth*pix2mm,'k','LineWidth',1); hold on;
    plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.', 'MarkerSize',15)
    plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.', 'MarkerSize',15)
end

xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');
%%