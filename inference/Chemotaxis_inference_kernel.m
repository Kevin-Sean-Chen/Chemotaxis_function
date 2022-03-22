%% mGLM fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load tracks and maps
%%% given Paths as cells with all 2D-trajectories
%%% and Fcon as a function that maps 2D position to an odor value readout

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/Landscape.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211029_GWN_app+_MEK110mM_40ml/OdorFx.mat');
Fcon = Fcon.F;

Cmap = load('/home/kschen/github/OdorSensorArray/OSA_MFC_PID_scripts/Landscape_low.mat');
Fcon = load('/home/kschen/github/OdorSensorArray/OSA_MFC_PID_scripts/OdorFx_low.mat');
Fcon = Fcon.F;

Cmap = Cmap.vq1;
Cmap = fliplr(flipud(Cmap));  %%% for inverted camper!

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

% Paths = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/App+_tracks.mat'); Tracks = Paths.Tracks;
% Tracks = deleted_tracks;

%%
figure;
%%% pre-processing parameters
poly_degree = 3;  %polynomial fit for the tracks
bin = 4;  %temporal binning  %~0.5 s
filt = 7;  %filtering tracks
l_window = 2;  %lag time

%%% gather time series
allas = []; %angles
alldcp = [];  %perpendicular dC
alldC = [];  %tangential dC
alldis = [];  %displacements

for c = 1:length(Tracks) %length(Paths)
    
    %%% process path to angles and distances
%     temp = Paths{c};  %for loading deleted tracks
    temp = Tracks(c).Path;  %for loading saved tracks
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    vecs = diff(subs);
    angs = zeros(1,size(vecs,1));    

    %%% for time step calculations
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    
    %%% iterate through worms
    for dd = (l_window+1):length(angs)
        %%% angle function
%         angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        
        %%% perpendicular concentration change
        perp_dir = [-vecs(dd-l_window,2), vecs(dd-l_window,1)];
        perp_dir = perp_dir/norm(perp_dir);
        dCp(dd) = Fcon(subs(dd-l_window,1)+perp_dir(1)*1., subs(dd-l_window,2)+perp_dir(2)*1)...
             - Fcon(subs(dd-l_window,1)-perp_dir(1)*1, subs(dd-l_window,2)-perp_dir(2)*1);
        
        %forward concentration change
        dCs(dd) = Fcon(subs(dd,1),subs(dd,2)) - Fcon(subs(dd-l_window,1),subs(dd-l_window,2));
        
        %%%check displacement
        dds(dd) = norm(subs(dd,1)-subs(dd-l_window,1), subs(dd,2)-subs(dd-l_window,2));
        %norm(vecs(dd,:));
%         dd
%         dds(dd)
%         pause();

    end
    %remove zeros
    dCp = dCp((l_window+1):end);
    dCs = dCs((l_window+1):end);
    dds = dds((l_window+1):end);
    angs = angs((l_window+1):end);
    
    allas = [allas angs];
    alldcp = [alldcp dCp];
    alldC = [alldC dCs];
    alldis = [alldis dds];  %displacement (effective velocity)
    
    
end

%%
%     K_ = THETA(1);  %variance of von Mises
%     A_ = THETA(2);  %max turn probability
%     B_ = THETA(3:8);  %kernel for dC transitional concentration difference (weights on kerenl basis)
%     C_ = THETA(9);  %baseline turning rate
%     Amp = THETA(10);  %amplitude for kernel for dcp normal concentration difference
%     tau = THETA(11);  %time scale for dcp kernel
%     K2_ = THETA(12);  %vairance of the sharp turn von Mises

%% test with stats-model for kernels
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);
ang_fit = allas(1:end-1);
dcp_fit = alldcp(1:end-1);
ddc_fit = (alldC(1:end-1));
lfun = @(x)nLL_kernel_chemotaxis(x,ang_fit, dcp_fit, ddc_fit, cosBasis, 0);
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [1e-5, 1e-5, ones(1,nB)*-inf, 1e-5 -inf, 1e-5, 1e-5];
UB = [10, 0.04, ones(1,nB)*inf, 0.01, inf, 100, 30];
% prs0 = rand(1,10);
prs0 = [1, 0.01, randn(1,nB)*10, 0.01, -10, 100, 0.1] + randn(1,10)*0.;
[x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

x
fval

%% test parameter distribution
repeats = 10;
prs_rep = zeros(length(prs0), repeats);
for rr =1:repeats
    % load observation
    ang_fit = allas(1:end-1);
    dcp_fit = alldcp(1:end-1);
    ddc_fit = (alldC(1:end-1));
    lfun = @(x)nLL_kernel_chemotaxis(x,ang_fit, dcp_fit, ddc_fit, cosBasis);
    % fitting
    prs0 = [1, 0.1, randn(1,nB), 0.01, -1, 100, 1];  %initialization 
    prs0 = prs0+prs0.*randn(1,10)*.1;  %perturbation
    [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[]); 
    
    prs_rep(rr,:) = x;
end

figure()
plot(prs_naive(1:5,:)','-o')
%% analysis
var_names = {'K_1','A','W','B','\alpha','\tau','K_2'};
var_app = zeros(repeats,length(var_names));
var_nai = zeros(repeats,length(var_names));
for rr = 1:repeats
    var_app(rr,:) = [prs_app(rr,1:2),norm(prs_app(rr,3:6)),prs_app(rr,7:10)];
    var_nai(rr,:) = [prs_naive(rr,1:2),norm(prs_naive(rr,3:6)),prs_naive(rr,7:10)];
end

y = [mean(var_app) ; mean(var_nai)]';
err = [std(var_app) ; std(var_nai)]';
figure
bar(y)
set(gca, 'XTickLabel', var_names);
hold on
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), '.');
end
hold off

%% training and testing??
samp_i = 1000:4000;
ang_samp = allas(samp_i);
dcp_samp = alldcp(samp_i);
ddc_samp = alldC(samp_i);
pred_ang = P_dth_dC(mean(prs_app),dcp_samp,ddc_samp);

figure()
plot(ang_samp); hold on
plot(pred_ang+280)
legend({'data','prediction'})
xlabel('time steps')
ylabel('angle')

figure()
plot(ang_samp, pred_ang,'o')

%% parameter tuning
num_sets = 5;
set_len = floor(size(allas,2)/num_sets);
ang_tt = reshape(allas(1:num_sets*set_len),num_sets,set_len);  %%groups for train and test
dcp_tt = reshape(alldcp(1:num_sets*set_len),num_sets,set_len);
ddc_tt = reshape(alldC(1:num_sets*set_len),num_sets,set_len);

regs = [0,1,10,50,100,200];  %regularization
TT_temp = 1:num_sets;  %for train test indexing
cros_vals = zeros(num_sets,length(regs));
for rr = 1:length(regs)
    
for ii = 1:num_sets
    % load training data
    ang_fit = ang_tt(ii,:);
    dcp_fit = dcp_tt(ii,:);
    ddc_fit = ddc_tt(ii,:);
    lfun = @(x)nLL_kernel_chemotaxis(x,ang_fit, dcp_fit, ddc_fit, cosBasis,regs(rr));
    % fitting
    prs0 = [1, 0.1, randn(1,nB), 0.01, -1, 100, 1];  %initialization 
    prs0 = prs0+prs0.*randn(1,10)*.1;  %perturbation
    [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[]); 
    
    % test predictions
    this_tt = TT_temp;
    this_tt(ii) = [];
    ang_test = ang_tt(this_tt,:);
    dcp_test = dcp_tt(this_tt,:);
    ddc_test = ddc_tt(this_tt,:);
%     pred_ang = P_dth_dC(x,dcp_test(:)',ddc_test(:)');
%     M_ = corrcoef(pred_ang, ang_test(:)');
    %for training%
    pred_ang = P_dth_dC(x,dcp_tt(:)',ddc_tt(:)');
    M_ = corrcoef(pred_ang, ang_tt(:)');
    
    cros_vals(ii,rr) = M_(1,2);
    
    rr
    ii
    
end
    
end

%%
figure()
plot(regs,TEST,'k-*')
hold on
plot(regs,cros_vals)  %training
xlabel('\lambda')
ylabel('Corr')
title('black:test, color:train')

%% sufficient statistics
K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7)/1; Amp = x(8); tau = x(9); K2_ = x(10);

figure
subplot(2,2,1)
xx = 1:length(cosBasis);
plot(Amp*exp(-xx/tau))
title('\delta C^{\perp} kernel')
subplot(2,2,2)
plot(B_ * cosBasis')
title('\delta C kernel')

subplot(2,2,3)
filt_dcp = conv(dcp_fit, Amp*exp(-xx/tau), 'same');
[aa,bb] = hist((ang_fit - filt_dcp)*pi/180 , 300);
bar( bb, 1/(2*pi*besseli(0,K_^2)) * exp(K_^2*cos( bb )) , 100); hold on
bar( bb, 1/(2*pi*besseli(0,K2_^2)) * exp(K2_^2*cos( bb-pi )) , 100,'r');
title('von Mises for \delta C^{\perp}')
subplot(2,2,4)
filt_ddc = conv( ddc_fit, (B_ * cosBasis'), 'same' );
Pturns = A_ ./ (1 + exp( filt_ddc)) + C_;
plot(filt_ddc , Pturns,'o')
title('Logistic for \delta C')


%%% printing
disp(['K_=',num2str(K_),' K2_',num2str(K2_),'beta',num2str(sum(B_)),'alpha',num2str(Amp)])
%% Generative model, with inferred parameters~
figure;
imagesc(Cmap); hold on;

% reconstruct parameters 
kappa = (1/K_)^1*(180/pi);  
A = A_*1;
Kdcp = Amp*exp(-xx/tau);
Kddc = B_ * cosBasis';
Pturn_base = C_;
kappa2 = (1/K2_)^1*(180/pi);  
wind = size(cosBasis,1);

origin = [size(Cmap,2)/2,size(Cmap,1)/2];
alldths = [];
alldCs = [];
alldCps = [];
allths = [];

REP = 20;
T = 3000;
trackss = zeros(REP,T,2);

for rep = 1:REP
    
% T = 3000;
dt = 1;
vm = 0.2*33;  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
vs = .2;
tracks = zeros(T,2);
tracks(1,:) = [rand()*size(Cmap,2), rand()*size(Cmap,1)];%origin; %initial position
tracks(2,:) = tracks(1,:)+randn(1,2)*vm*dt;%origin+randn(1,2)*vm*dt;
ths = zeros(1,T);  ths(1:3) = randn(1,3)*360; %initial angle
dxy = randn(1,2);  %change in each step

dths = zeros(1,T); %"recordings"
dcs = zeros(1,T);
dcps = zeros(1,T);
dCv = zeros(1,wind);
dCpv = zeros(1,wind);
for t = 1+2:T
    %%% dC in a window
    dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) - Fcon(tracks(t-2,1),tracks(t-2,2)) ,  dCv];  % step-wise derivative
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) ,  dCv];  % absolute concentration
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) - dCv(end),   dCv];  % moving threshold
    dCv = dCv(1:end-1);  %remove out of window
    perp_dir = [-dxy(2) dxy(1)];
    perp_dir = perp_dir/norm(perp_dir);
    %%% dCp in a window
    dCpv = [Fcon(tracks(t-1,1)+perp_dir(1),tracks(t-1,2)+perp_dir(2)) - Fcon(tracks(t-1,1)-perp_dir(1),tracks(t-1,2)-perp_dir(2)) ,  dCpv];
    dCpv = dCpv(1:end-1);
    
    wv = -sum(Kdcp.*dCpv)/1 + kappa*randn;%length(wind)
    P_event = A/(1+exp( (sum(Kddc.*dCv)/1) *dt)) + Pturn_base;%length(wind)
    if rand < P_event*1
        beta = 1;
    else
        beta = 0;
    end
    rt = beta*[(randn*kappa2-180) + 0*(rand*360-180)];
    dth = wv+rt;
    if dth>180; dth = dth-180; end;  if dth<-180; dth = dth+360; end  %within -180~180 degree range
    
    dths(t) = dth;
    dcs(t) = dCv(1);
    dcps(t) = dCpv(1);
    
    vv = vm+vs*randn;
    vv = alldis(randi(length(alldis)));
    if vv==0
        vv = 0.1;
    end
    ths(t) = ths(t-1)+dth*dt;
%     if ths(t)>180; ths(t) = ths(t)-180; end;  if ths(t)<-180; ths(t) = ths(t)+360; end  %within -180~180 degree range
    e1 = [1,0];
    vec = [tracks(t-1,1)  tracks(t-1,2)]-origin; %current vector
    theta = acosd(max(-1,min((vec*e1')/norm(vec)/norm(e1),1)));  %current angle
    dd = [vv*sin(ths(t)*pi/180) vv*cos(ths(t)*pi/180)];
    R = [cos(theta*pi/180) sin(theta*pi/180); -sin(theta*pi/180) cos(theta*pi/180)];
    dxy = (R)*dd';

    tracks(t,1) = tracks(t-1,1)+dxy(1)*dt;
    tracks(t,2) = tracks(t-1,2)+dxy(2)*dt;
        
end

plot(tracks(:,1),tracks(:,2))
hold on
plot(tracks(1,1), tracks(1,2),'ko') %(origin(1),origin(2),'ko')%

alldths = [alldths dths];
alldCs = [alldCs dcs];
alldCps = [alldCps dcps];
allths = [allths ths];

trackss(rep,:,:) = tracks;

end

%% Check statistics
figure()
subplot(3,2,1); hist(allas,100);xlim([-180,180]);  ylabel('\delta \theta'); title('data'); subplot(3,2,2); hist(alldths,100); xlim([-180,180]);title('model');
subplot(3,2,3); hist(alldC,100);  ylabel('\delta C');      subplot(3,2,4); hist(alldCs,100)
subplot(3,2,5); hist(alldcp,100); ylabel('\delta C^{\perp}');     subplot(3,2,6); hist(alldCps,100)

%%
figure()
for ii =1:REP
    plot(trackss(ii,:,1)-trackss(ii,1,1),trackss(ii,:,2)-trackss(ii,1,2));
    hold on
end
plot([0,0],[0,0],'k.','MarkerSize',40)
hold off
set(gca,'Ydir','reverse')


%% Reconstruction from inferred kernel
bins = 30;
dC_ = dcp_fit;
dA_ = filt_dcp*14*l_window;

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

%% 
figure;
Pturns = A_ ./ (1 + exp( filt_ddc)) + C_;
plot(filt_ddc/length((B_ * cosBasis')) , Pturns,'o')

xlabel('\delta C', 'Fontsize',20)
ylabel('P(turn|\delta C)', 'Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);
