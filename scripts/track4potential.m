%%% track4potential
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')

%batch analysis
fields_to_load = {'Path','SmoothX','SmoothY','Time','Runs','Pirouettes'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% visualized tracks
poly_degree = 3;  %polynomial fitting for moving window
filt = 7;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 30*1; %minimum time in seconds
minx = 10;  %minimum displacement (in terms of pixels)
pix2mm = 1/31.5;
conc_thre = 0.8*max(max(M));

figure();
cand = [];  %index of tracks as candidates for analysis
alldists = [];
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); %pixel displacement
        alldists = [alldists displace*pix2mm^2];  %all dispacements in mm
        if displace > minx^2  %space cutoff
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

xlabel('x (mm)'); ylabel('y (mm)');
set(gca,'Fontsize',20); set(gcf,'color','w');
set ( gca, 'ydir', 'reverse' )

%% save tracks (for python loading)
save_dir = '/projects/LEIFER/Kevin/';
Track_data = Tracks(cand);
save([save_dir,'tracks_ave2.mat'],'Track_data')

