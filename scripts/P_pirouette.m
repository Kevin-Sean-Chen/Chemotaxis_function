%%% P_Pirouettes
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
Fcon = Fcon.F;

%% compute dc before pirouettes
dc_window = 14*2.;  %7bins for 0.5s in real time
fr2sec = 14/dc_window;  %frames to seconds in real time
all_dC = [];
all_pr = [];

for c = 1:length(Tracks)
    %%% load observations
    prs = Tracks(c).Pirouettes;
    paths = Tracks(c).Path;
    
    %%% compute dC
    dCi = [];
    for dd = dc_window+1:dc_window:size(paths,1)
        dCi = [dCi  (Fcon(paths(dd-dc_window,1),paths(dd-dc_window,2)) - Fcon(paths(dd,1),paths(dd,2)))*fr2sec];  %record dC for this track
    end
    
    %%% mark pirouettes
    time_i = dc_window+1:dc_window:size(paths,1);
    pr_event = zeros(1,length(time_i));
    if isempty(prs)==0
        for pp=1:size(prs,1)
            [aa,pos] = min(abs(prs(pp,1)-time_i));  %find closest time bin
            pr_event(pos) = 1;
        end
    end
    %%% recording
    all_dC = [all_dC dCi];
    all_pr = [all_pr pr_event];
end

%% adaptive binning
bins = 30;
num_points = floor(length(all_dC)/bins);  %number of points in one bin, now that we fix this value
[val,pos] = sort(all_dC);
[nnn,xout] = hist(all_dC,bins); 

dc_i = zeros(1,bins);  %median dC
cnt_i = zeros(1,bins);  %count
p_cnts = zeros(1,bins);  %pr/s
for bi = 1:bins
    tempb = pos((bi-1)*num_points+1:(bi)*num_points);  %index of for this bin
    dc_i(bi) = median(all_dC(tempb));  %take the median for bin center
    p_cnts(bi) = sum(all_pr(tempb)) / (length(tempb)*fr2sec);  %probability of turns in this bin
    cnt_i(bi) = sum(all_pr(tempb));
end

%%
%%% calculate error bars
cnt_b = nnn;
EE = ( ((cnt_i-1)./(cnt_b.^2)) + (cnt_i.^2.*(cnt_b-1)./(cnt_b.^4)) ).^0.5 *1/14;

figure()
plot(dc_i, p_cnts, '-o')
% errorbar(dc_i, p_cnts, EE,'-o','LineWidth',2)
xlabel('\delta C/s', 'Fontsize',20)
% ylabel('P(turn|\delta C)', 'Fontsize',20)
ylabel('P(pirouette) (event/s)', 'Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for runs, for weathervaning strategies
dcp_window = 14*2;  %7bins for 0.5s in real time
fr2sec = 14/dcp_window;  %frames to seconds in real time
pix2mm = 1/31.5;  %image to real scale
all_dcp = [];
all_ang = [];

for c=1:length(Tracks)
    %%% load observations
    runs = Tracks(c).Runs;
    paths = Tracks(c).Path;
    
    %%% loop through runs
    dcp = [];
    angs = [];
    for rr = 1:size(runs,1)
        path_i = paths(runs(rr,1):runs(rr,2),:);  %segment of running
        vecs = diff(path_i);  %vector of motions
        
        %%% compute dcp and d_angle
        for dd = dcp_window+1:dcp_window:size(vecs,1)
            disti = norm(path_i(dd-dcp_window,:)-path_i(dd,:));  %per distance
            
            perp_dir = [-vecs(dd-dcp_window,2), vecs(dd-dcp_window,1)];
            perp_dir = perp_dir/norm(perp_dir);
            dcpt = Fcon(path_i(dd-dcp_window,1)+perp_dir(1)*1., path_i(dd-dcp_window,2)+perp_dir(2)*1)...
                 - Fcon(path_i(dd-dcp_window,1)-perp_dir(1)*1, path_i(dd-dcp_window,2)-perp_dir(2)*1);  %computing normal direction dcp
            dcp = [dcp dcpt/1 ];
            
            angs = [angs  angles(vecs(dd-dcp_window,:)/norm(vecs(dd-dcp_window,:)),vecs(dd,:)/norm(vecs(dd,:))) /disti/pix2mm ];
        end
        
    end
    
    %%% recording
    all_dcp = [all_dcp dcp];
    all_ang = [all_ang angs];
end

%%
bins = 30;
dC_ = all_dcp;
dA_ = all_ang;

avang = zeros(1,bins);
stdang = zeros(1,bins);
[cts,bis] = hist(dC_,bins);
for bi = 2:length(bis)
    pos = intersect(find(bis(bi-1)<dC_),find((bis(bi)>=dC_)));
    if isempty(pos)~=1
        avang(bi) = nanmean(dA_(pos));
        stdang(bi) = nanstd(dA_(pos));
    end
end

% EE = ( ((cnt_i-1)./(cnt_b.^2)) + (cnt_i.^2.*(cnt_b-1)./(cnt_b.^4)) ).^0.5 *1/14;

figure()
set(gca,'FontSize',20);
errorbar((bis),avang,stdang)
figure()
plot((bis),avang,'-o')
xlabel('\delta C^{\perp}', 'Fontsize',20)
ylabel('<\delta \theta>','Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);

