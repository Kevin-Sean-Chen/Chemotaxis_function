%%% OdorFlow_3
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothSpeed'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% pre-processing visualization
% 
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');

Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_110mM.mat');
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_110mM.mat');

Fcon = Fcon.F;
M = Cmap.vq1;
M = (flipud(M));  %flipped camera

%%
poly_degree = 3;  %polynomial fitting for moving window
filt = 7;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 60*3; %minimum time in seconds
minx = 100;  %minimum displacement (in terms of pixels)
endingt = 60*20;  %only taking the first few minutes
pix2mm = 1/31.5;
targ_track = 7;
conc_thre = 200;

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

hold on
pos = find(Tracks(targ_track).Time<endingt);
x_smooth = smooth(Tracks(targ_track).Path(pos,1), filt,'sgolay',poly_degree);
y_smooth = smooth(Tracks(targ_track).Path(pos,2), filt,'sgolay',poly_degree);
plot(x_smooth*pix2mm, y_smooth*pix2mm,'W','LineWidth',1); hold on;
plot(x_smooth(1)*pix2mm, y_smooth(1)*pix2mm,'g.', 'MarkerSize',15)
plot(x_smooth(end)*pix2mm, y_smooth(end)*pix2mm,'r.', 'MarkerSize',15)

xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');
set ( gca, 'xdir', 'reverse' )

%%  scan trough tracks
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

%% show individual d_C and d_theta data
track = Tracks(7).Path(1:end,:);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis... summary statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% replace with real landscape later
bias = 0;
xx = linspace(0,15,3000);
yy = linspace(-7,7,2500);
[X,Y] = meshgrid(xx,yy);
v = 0.5;   %flow in chamber
D = 0.08;  %MEK
Def = D*(1 + 6.25^2/7.5);  %effective diffusion
C0 = 200;  %initial odor source
% Z = C0/2*(1-erf(Y.^2 ./ (2*sqrt(Def.*X./v)) ));  %1D
Z = C0/1*(1-erf(X ./ (2*sqrt(Def.*X./v)*5) )) .* exp(-(Y-bias).^2./(4*Def*X./v));  %combine
% Z = C0./(4*pi*Def.*X./v).^0.1 .*(1-erf(Y.^2./(4*Def*X./v)*1));   %solve with BC!!!

figure;
imagesc(Z)

%% contour and quiver plot of the landscape gradient
figure
contour(xx,yy,Z)
hold on
x_ = linspace(0,15,15);
y_ = linspace(-7,7,15);
z_ = C0/1*(1-erf(x_ ./ (2*sqrt(Def.*x_./v)*5) )) .* exp(-(y_-bias)'.^2./(4*Def*x_./v));
[px,py] = gradient(z_);
quiver(x_,y_,px,py)

%% analysis loop
nn = 1;
fixlen = 47;  %we only compute a 2d vector when it moves forward this much
[fx,fy] = gradient(fliplr((Z)),1);
% [fx,fy] = gradient((conv2(M,ones(100,100),'same')/10000),1);  %prepare 2D gradient field
all_grad = [];
all_angs = [];

figure()
for nn = 1:length(Tracks)
    if isempty(find(cand==nn))==0

vec_i = [];    %2d vectors with adaptive sampling rate
grad_i = [];   %at this location, what was the gradient direction

speed_i = Tracks(nn).SmoothSpeed;
path_i = Tracks(nn).Path;

ii = 1;
pos_temp = path_i(1,:);  %initial location
while ii<length(speed_i)
    delta_t = min(round(fixlen/speed_i(ii)), floor(length(speed_i)/2));  %discrete time window to update, adpated to the velocity
    vec_i = [vec_i; path_i(ii,:)-pos_temp];  %compute dispacement vector
    grad_i = [grad_i; [fx(floor(path_i(ii,2)),floor(path_i(ii,1))), fy(floor(path_i(ii,2)),floor(path_i(ii,1)))] ];  %gradident direction at this location
    pos_temp = path_i(ii,:);  %update postion
    ii = ii + delta_t+1;
end

angs = zeros(1, length(vec_i)-1);
grad = angs*1;
for pp = 1:length(vec_i)-1
    angs(pp) = angles(vec_i(pp,:)/norm(vec_i(pp,:)),vec_i(pp+1,:)/norm(vec_i(pp+1,:)));
    grad(pp) = angles(vec_i(pp,:)/norm(vec_i(pp,:)), grad_i(pp,:)/norm(grad_i(pp,:)));
    end

% subplot(121); plot(grad); hold on; plot(angs); hold off;
% subplot(122); plot(path_i(:,1),path_i(:,2))
% pause();

all_grad = [all_grad  grad];
all_angs = [all_angs  angs];
nn
    end
end


%% BRW
figure()

thre_pr = 100;   %threshold for a piroutte
thre_grad = 60;  %threshold for a gradient
cnt_i = []; cnt_b = [];
pos_0 = find(abs(all_grad)<thre_grad);
angs_0 = find(abs(all_angs(pos_0))>thre_pr);
aa = length(angs_0)/length(pos_0);  cnt_i = [cnt_i length(angs_0)]; cnt_b = [cnt_b length(pos_0)];
pos_0 = find(abs(all_grad)>thre_grad & abs(all_grad)<120);
angs_0 = find(abs(all_angs(pos_0))>thre_pr);
bb = length(angs_0)/length(pos_0);  cnt_i = [cnt_i length(angs_0)]; cnt_b = [cnt_b length(pos_0)];
pos_0 = find(abs(all_grad)>120);
angs_0 = find(abs(all_angs(pos_0))>thre_pr);
cc = length(angs_0)/length(pos_0);  cnt_i = [cnt_i length(angs_0)]; cnt_b = [cnt_b length(pos_0)];

EE = ( ((cnt_i-1)./(cnt_b.^2)) + (cnt_i.^2.*(cnt_b-1)./(cnt_b.^4)) ).^0.5 *1/1;

errorbar([1,2,3], [aa,bb,cc], EE,'ko','LineWidth',2); hold on
xlabel('gradient angles');  ylabel('turn probability')
bar([1,2,3],[aa,bb,cc])
xticks([1 2 3]); xticklabels({'<60','60~120','>120'})
set(gca,'Fontsize',20); set(gcf,'color','w');

%% WV!
nbs = 15;
figure()
pos_0 = find((all_grad)>-110 & all_grad<-70);
pos_60 = find(abs(all_grad)<30);
pos_120 = find((all_grad)>70 & all_grad<110);
H1 = histogram(all_angs(pos_0), nbs, 'Normalization', 'pdf'); hold on
H2 = histogram(all_angs(pos_60), nbs, 'Normalization', 'pdf'); hold on
H3 = histogram(all_angs(pos_120), nbs, 'Normalization', 'pdf');
% close(fig);

figure()
aa = H1.Values;  bb = H1.BinEdges;
plot(bb(2:end),aa)
mean(aa*bb(2:end)'); y = skewness(all_angs(pos_0))
hold on
aa = H2.Values;  bb = H2.BinEdges;
plot(bb(2:end),aa)
mean(aa*bb(2:end)'); y = skewness(all_angs(pos_60))
hold on
aa = H3.Values;  bb = H3.BinEdges;
plot(bb(2:end),aa)
mean(aa*bb(2:end)'); y = skewness(all_angs(pos_120))