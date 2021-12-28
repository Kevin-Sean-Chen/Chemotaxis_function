%Weather-Vaning Analysis
%%%
%Taking Runs from the pipeline and
%%%
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% pre-processing
% M = imread('Z:\Kevin\20190817\Data20190817_165031\Frame_000000.jpg');
% M = imread('Z:\Kevin\20191113_GWN_N2_naive\Data20191113_140408\Frame_000000.jpg');
test = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211013_biased_110mM/Landscape.mat');
M = test.vq1;
pix2mm = 1/31.5;  %pixel to mm (camera position before the flow chanber setup) %%%16.5 for new camera and 31.5 for old one
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);

poly_degree = 3;  %polynomial fitting for moving window
filt = 15;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 60*3; %minimum time in seconds
minx = 100;  %minimum displacement (in terms of pixels)
target = [2517,975];  %position of target/sourse of odorant (approximated from images)
endingt = 60*15;  %only taking the first few minutes

hold on
cand = [];  %index of tracks as candidates for analysis
alldists = [];
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); %pixel displacement
        alldists = [alldists displace*pix2mm^2];  %all dispacements in mm
        if displace > minx^2  %space cutoff
            pos = find(Tracks(i).Time>endingt);  %time window cutoff (the later time points are less correct...)
            if isempty(pos)~=1
                x_smooth = smooth(Tracks(i).Path(pos,1), filt,'sgolay',poly_degree);
                y_smooth = smooth(Tracks(i).Path(pos,2), filt,'sgolay',poly_degree);
                plot(x_smooth*pix2mm, y_smooth*pix2mm,'k'); hold on;
                cand = [cand i];
            end
        end
    end
end


%% %%% by time %%%
figure
delay_scan = 1:100:501;
for ll = 1:length(delay_scan)
ll = delay_scan(ll);
%%
bin = 5;
filt = 30;
l_window = ll;  %lag time (real lag is bin*l_window*fr seconds)
ct = 0;
allas = [];
alldB = [];
alldis = [];

timelim = 60*10;

for c = 1:length(cand)
    
    id = cand(c);
    temp0 = Tracks(id).Path;
    %%%only look at runs
    runs = Tracks(id).Runs;
    temp = [];
    for rr = 1:size(runs,1)
        temp = [temp; temp0(runs(rr,1):runs(rr,2),:)];  %intervals for runs
    end
    %%%
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = diff(subs);
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
        
    %%%for time step calculations
    dBs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    for dd = l_window+1:length(angs)
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/(norm(vecs(dd,:)))*pix2mm);  %angle between past and now
        dBs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),(subs(dd-l_window,:)-target)/norm(subs(dd-l_window,:)-target));  %bear angle back then
        dds(dd) = norm(vecs(dd,:));
    end
    %remove zeros
    angs = angs(l_window+1:end);
    dBs = dBs(l_window+1:end);
    dds = dds(l_window+1:end);
    
    allas = [allas angs];
    alldB = [alldB dBs];
    alldis = [alldis dds];
 
end

% figure
% plot(alldB,allas,'o')

% %% removing small jitterings (test)
% vec_threshold = 0.1;
% vec_marks = find(alldis > vec_threshold);
% %% removing large anlge changes (test)
% ang_threshold = 20;
% ang_marks = find(abs(allas) < ang_threshold);
% va_marks = intersect(vec_marks, ang_marks);

%% plot adaptive binning
bins = 60;
% dC_ = alldB;
% dA_ = allas;
dC_ = alldB(va_marks);
dA_ = allas(va_marks);

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
errorbar(bis,avang,stdang)

hold on
end
legendInfo = {};
for d = 1:length(delay_scan)
    legendInfo{d} = ['delay=', num2str(delay_scan(d)*fr*bin),'s'];
end
legend(legendInfo)

