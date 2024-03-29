%%% Statistical_description_flow
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% pre-processing visualization
% test = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/Landscape.mat');
% Cmap =load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/Landscape.mat');
% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/OdorFx.mat');
% 
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20220113_GWN_app+_MEK110mM_gasphase_30ml_200air/Landscape.mat');
% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20220113_GWN_app+_MEK110mM_gasphase_30ml_200air/OdorFx.mat');

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');

% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_110mM.mat');
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_110mM.mat');

% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_low.mat');
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_low.mat');

Fcon = Fcon.F;
M = Cmap.vq1;
M = (flipud(M));  %flipped camera

%%
poly_degree = 3;  %polynomial fitting for moving window
filt = 7;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 60*1; %minimum time in seconds
minx = 100;  %minimum displacement (in terms of pixels)
target = [2517,975];  %position of target/sourse of odorant (approximated from images)
endingt = 60*30;  %only taking the first few minutes
pix2mm = 1/31.5;

figure();
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
hold on
cand = [];  %index of tracks as candidates for analysis
alldists = [];
for i = 1:nn
%     if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); %pixel displacement
        alldists = [alldists displace*pix2mm^2];  %all dispacements in mm
        if displace > minx^2  %space cutoff
            pos = find(Tracks(i).Time<endingt);  %time window cutoff (the later time points are less correct...)
%             if isempty(pos)~=1
                x_smooth = smooth(Tracks(i).Path(pos,1), filt,'sgolay',poly_degree);
                y_smooth = smooth(Tracks(i).Path(pos,2), filt,'sgolay',poly_degree);
                plot(x_smooth*pix2mm, y_smooth*pix2mm,'k','LineWidth',1); hold on;
                plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.')
                plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.')
                cand = [cand i];
            end
%         end
%     end
end

xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');
set ( gca, 'xdir', 'reverse' )
%%
figure();
for i = 1:nn
    pos = find(Tracks(i).Time<endingt);  %time window cutoff (the later time points are less correct...)
    x_smooth = smooth(Tracks(i).Path(pos,1), filt,'sgolay',poly_degree);
    y_smooth = smooth(Tracks(i).Path(pos,2), filt,'sgolay',poly_degree);
    imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on;
    plot(x_smooth*pix2mm, y_smooth*pix2mm,'k');
    plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'r*')
    title(num2str(i))
    pause();
end

%%
track = Tracks(77).Path(1:end,:);
bin = 7;
x_smooth = smooth(track(:,1), filt,'sgolay',poly_degree);
y_smooth = smooth(track(:,2), filt,'sgolay',poly_degree);
temp = [x_smooth'; y_smooth']';
subs = temp(1:bin:end,:);
vecs = diff(subs);
Ct = zeros(1,length(vecs));
dtht = Ct*1;
for pp = 2:length(vecs)
%     Ct(pp) = Fcon(subs(pp,1), subs(pp,2));
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


%% calculate aboservables
figure;
poly_degree = 3;  %polynomial fit for the tracks
bin = 4;  %temporal binning
filt = 7;  %filtering tracks
l_window = 2;  %lag time
allas = []; %angles
alldcp = [];  %perpendicular dC
alldC = [];  %tangential dC
alldis = [];  %displacements
alldata = struct();

for c = 1:nn %length(del_cls)%(Paths)
    
    %%% process path to angles and distances
%     temp = values{c};
%     temp = deleted_tracks(c).Path;
    temp = Tracks(c).Path;
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    vecs = diff(subs);
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    

    %%% for time step calculations
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    
    %%% iterate through worms
    for dd = l_window+1:length(angs)
        %%% angle function
%         angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        
        %%% perpendicular concentration change
        perp_dir = [-vecs(dd-l_window,2), vecs(dd-l_window,1)];
        perp_dir = perp_dir/norm(perp_dir);
        dCp(dd) = Fcon(subs(dd-l_window,1)+perp_dir(1)*1., subs(dd-l_window,2)+perp_dir(2)*1)...
             - Fcon(subs(dd-l_window,1)-perp_dir(1)*1, subs(dd-l_window,2)-perp_dir(2)*1);
%         dCp(dd) = Est_con(subs(dd-l_window,1)+perp_dir(1)*1., subs(dd-l_window,2)+perp_dir(2)*1, target(1), target(2), 50)...
%                   -Est_con(subs(dd-l_window,1)-perp_dir(1)*1, subs(dd-l_window,2)-perp_dir(2)*1, target(1), target(2), 50);
        
        %forward concentration change
        dCs(dd) = Fcon(subs(dd,1),subs(dd,2)) - Fcon(subs(dd-l_window,1),subs(dd-l_window,2));
        
        %%%check displacement
        dds(dd) = norm(vecs(dd,:));
    
    end
    
    %remove zeros
    dCp = dCp(l_window+1:end);
    dCs = dCs(l_window+1:end);
    dds = dds(l_window+1:end);
    angs = angs(l_window+1:end);
    
    allas = [allas angs];
    alldcp = [alldcp dCp];
    alldC = [alldC dCs];
    alldis = [alldis dds];  %displacement (effective velocity)

%%% storing data
% alldata(c).dth = angs;
% alldata(c).dcp = dCp;
% alldata(c).dC = dCs;
% alldata(c).dis = dds; 

end

%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% weathervaning analysis
bins = 30;
dC_ = alldcp;
dA_ = allas;

avang = zeros(1,bins);
stdang = zeros(1,bins);
% h = histogram(dC_,bins);
% cts = h.Values;
% bis = h.BinEdges(1:end-1);
[cts,bis] = hist(dC_,bins);
for bi = 2:length(bis)
    pos = intersect(find(bis(bi-1)<dC_),find((bis(bi)>=dC_)));
    if isempty(pos)~=1
        avang(bi) = mean(dA_(pos));
        stdang(bi) = std(dA_(pos));
    end
end
figure()
set(gca,'FontSize',20);
errorbar(bis,avang,stdang)
figure()
plot(bis,avang,'-o')
xlabel('\delta C^{\perp}', 'Fontsize',20)
ylabel('<\delta \theta>','Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);

%% random turn analysis
bins = 30;
p_threshold = 120;  %hard threshold for random turns
dC_ = alldC;
dA_ = allas;

cnts = ones(1,bins);
p_cnts = zeros(1,bins);
[cts,bis] = hist(dC_,bins);
for bi = 2:length(bis)
    pos = intersect(find(bis(bi-1)<dC_),find((bis(bi)>=dC_)));
    if isempty(pos)~=1
        cnts(bi) = length(pos);  %number of points in this bin
        p_cnts(bi) = length(find(abs(dA_(pos)) > p_threshold));  %number of turns
    end
end
figure()
errorbar(bis,avang,stdang)
plot(bis(2:end),p_cnts(2:end)./cnts(2:end),'-o')
xlabel('dC')
ylabel('P(sharp turn)')

%% adaptive binning test
bins = 25;
p_threshold = 120;
num_points = floor(length(allas)/bins);  %number of points in one bin, now that we fix this value
[val,pos] = sort(dC_);

cnts = zeros(1,bins);
p_cnts = zeros(1,bins);
for bi = 1:bins
    tempb = pos((bi-1)*num_points+1:(bi)*num_points);  %index of for this bin
    cnts(bi) = median(dC_(tempb));  %take the median for bin center
    p_cnts(bi) = length(find(abs(dA_(tempb)) > p_threshold)) / num_points;  %probability of turns in this bin
end

figure()
plot(-cnts, p_cnts, '-o')
xlabel('\delta C', 'Fontsize',20)
ylabel('P(turn|\delta C)', 'Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);