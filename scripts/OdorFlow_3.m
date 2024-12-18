%%% OdorFlow_3
% simple visualization of navigation path from the raw data
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothSpeed','AngSpeed','SmoothX','SmoothY','LEDVoltages'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% pre-processing visualization

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623.mat')
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat')

% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_low.mat');
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_low.mat');
% % 
% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_110mM.mat');
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_110mM.mat');

% Cmap = load('/projects/LEIFER/Kevin/Publications/Chen_flow_2022/Landscape_cone_agar.mat');

Fcon = Fcon.F;
M = Cmap.vq1;
M = fliplr(flipud(M));  %flipped camera

%% for salt
[rows, cols] = size(M);
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
% [x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M = (y_*0+1) .* gradient_x;  figure; imagesc(M)

%% pand a,b 
poly_degree = 3;  %polynomial fitting for moving window
filt = 14;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 60*2.;%3 %minimum time incand seconds
minx = 100*1;  %minimum displacement (in terms of pixels) %3 for app?
endingt = 60*30;  %only taking the first few minutes
pix2mm = 1/31.5;
targ_track = 7;
conc_thre = .9*max(max(M));

figure();
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
% colormap(gray)
hold on
cand = [];  %index of tracks as candidates for analysis
alldists = [];
ttt= 0;
for i = 1:nn %200
    ttt = ttt+ (Tracks(i).Time(end)-Tracks(i).Time(1));
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); %pixel displacement
        alldists = [alldists displace*pix2mm^2];  %all dispacements in mm
        if displace > minx^2  %space cutoff
            pos = find(Tracks(i).Time<endingt);  %time window cutoff (the later time points are less correct...)
            if isempty(pos)~=1
                x_smooth = smooth(Tracks(i).Path(:,1), filt,'sgolay',poly_degree);
                y_smooth = smooth(Tracks(i).Path(:,2), filt,'sgolay',poly_degree);
                if M(floor(y_smooth(1)), floor(x_smooth(1))) < conc_thre %& x_smooth(end)>100 & x_smooth(end)<2900 & y_smooth(end)>100 & y_smooth(end)<2400
                        
                    plot(x_smooth*pix2mm, y_smooth*pix2mm,'k','LineWidth',1); hold on;
                    plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.', 'MarkerSize',30/2)
                    plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.', 'MarkerSize',30/2)
                    cand = [cand i];

                end
            end
        end
    end
end

hold on
pos = find(Tracks(targ_track).Time<endingt);
x_smooth = smooth(Tracks(targ_track).Path(pos,1), filt,'sgolay',poly_degree);
y_smooth = smooth(Tracks(targ_track).Path(pos,2), filt,'sgolay',poly_degree);
plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.', 'MarkerSize',15)
plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.', 'MarkerSize',15)

xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');

%%  scan trough tracks
figure();
for j = 1:length(cand)
    i = cand(j);
    pos = find(Tracks(i).Time<endingt);  %time window cutoff (the later time points are less correct...)
    x_smooth = smooth(Tracks(i).Path(pos,1), filt,'sgolay',poly_degree);
    y_smooth = smooth(Tracks(i).Path(pos,2), filt,'sgolay',poly_degree);
    subplot(121)
    imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on;
    plot(x_smooth*pix2mm, y_smooth*pix2mm,'k');
    plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'r*')
    title(num2str(i))
    subplot(122)
    yyaxis right;  plot(Tracks(i).AngSpeed); yyaxis left
    plot(Tracks(i).SmoothSpeed); 
    pause();
end

%% example tracks
cand_list = [94,65,87];
figure();
imagesc(M','XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on;
for ll = 1:length(cand_list)
    i = cand_list(ll);
    x_smooth = smooth(Tracks(i).Path(:,1), filt,'sgolay',poly_degree);
    y_smooth = smooth(Tracks(i).Path(:,2), filt,'sgolay',poly_degree);
    plot(x_smooth*pix2mm, y_smooth*pix2mm,'k','LineWidth',2); hold on
    plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.', 'MarkerSize',15)
    plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.', 'MarkerSize',15)
end
axis([20,100 -10,100]); set(gca,'YDir','normal')

%% example gradient climbing
i = 7; %picking an example
wind = 1100:2650;

speeds = Tracks(i).SmoothSpeed;
paths = Tracks(i).Path;
angspeed_i = Tracks(i).AngSpeed;
runs = Tracks(i).Runs;
x_smooth = smooth(paths(:,1), filt,'sgolay',poly_degree);
y_smooth = smooth(paths(:,2), filt,'sgolay',poly_degree);
path_i = [x_smooth   y_smooth];  

vec_i = diff(path_i,1);
grad_i = zeros(1,length(vec_i)-1);
angs = zeros(1, length(vec_i)-1);
for ii = 1:length(angs)-1
    angs(ii) = angles(vec_i(ii,:)/norm(vec_i(ii,:)),vec_i(ii+1,:)/norm(vec_i(ii+1,:)));
    grad_i(ii) = M(floor(y_smooth((ii))),floor(x_smooth((ii))));
end

xx = [1:length(wind)]*1/14;  %seconds
cts = zeros(1,length(wind));
for tt = 1:length(wind);  cts(tt) = M(floor(y_smooth(wind(tt))),floor(x_smooth(wind(tt)))); end
figure;
subplot(211);
scatter(x_smooth(wind),y_smooth(wind),[],cts, 'filled'); hold on;
for pp = 3:6
scatter( x_smooth(Tracks(i).Pirouettes(pp,1):Tracks(i).Pirouettes(pp,2)+15), y_smooth(Tracks(i).Pirouettes(pp,1):Tracks(i).Pirouettes(pp,2)+15) ,'k');  
end
plot(x_smooth(wind(1)), y_smooth(wind(1)),'g.', 'MarkerSize',15)
plot(x_smooth(wind(end)), y_smooth(wind(end)),'r.', 'MarkerSize',15)
axis off; colorbar();
subplot(212)
plot(xx, Tracks(i).AngSpeed(wind));  hold on
for pp = 3:6
x_points = [[Tracks(i).Pirouettes(pp,1), Tracks(i).Pirouettes(pp,1), Tracks(i).Pirouettes(pp,2)+15, Tracks(i).Pirouettes(pp,2)+15]-wind(1)]/14;  
y_points = [-400, 400, 400, -400]; %[-.1, .1, .1, -.1];
color = [0, 0, 1];
hold on;
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.1;
end

set(gca,'Fontsize',20); set(gcf,'color','w');
xlabel('time (s)')
ylabel('angular veolocity')

%% show individual d_C and d_theta data
track = Tracks(244).Path(1:end,:);
bin = 7;
x_smooth = smooth(track(:,1), filt,'sgolay',poly_degree);
y_smooth = smooth(track(:,2), filt,'sgolay',poly_degree);
temp = [x_smooth'; y_smooth']';
subs = temp(1:bin:end,:);
vecs = diff(subs);
Ct = zeros(1,length(vecs));
dtht = Ct*1;
for pp = 2:length(vecs)
    Ct(pp) = M(floor(subs(pp,2)), floor(subs(pp,1)));
    dtht(pp) = angles(vecs(pp-1,:)/norm(vecs(pp-1,:)),vecs(pp,:)/norm(vecs(pp,:)));
end
time_ = [1:length(vecs)]/14*bin;
figure;
subplot(211); plot(time_(2:end-1), Ct(2:end-1),'b', 'Linewidth',3); ylabel('ppm');
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);
xticks([]);
subplot(212); plot(time_(2:end-1), dtht(2:end-1), 'k', 'Linewidth',1); ylabel('d\theta'); xlabel('time (s)')
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);
