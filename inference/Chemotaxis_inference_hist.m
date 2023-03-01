%% mGLM fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load tracks and maps
%%% given Paths as cells with all 2D-trajectories
%%% and Fcon as a function that maps 2D position to an odor value readout

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');

% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_110mM.mat');
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_110mM.mat');

% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_cone_low.mat');
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_cone_low.mat');

Fcon = Fcon.F;
M = Cmap.vq1;
M = fliplr(flipud(M));  %flipped camera
H = fspecial('average',700);
M = imfilter(M, H, 'replicate');

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothX','SmoothY'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

cand = 1:length(Tracks);  % if all

%%
figure;
%%% pre-processing parameters
poly_degree = 3;  %polynomial fit for the tracks
bin = 5;  %temporal binning  %~0.5 s
filt = 14*1;  %filtering tracks
l_window = 1;  %lag time
perp_dist = 1;
fr2sec = 1/14;

%%% gather time series
allas = []; %angles
alldcp = [];  %perpendicular dC
alldC = [];  %tangential dC
alldis = [];  %displacements
alltrials = [];  %mark for tracks
allxys = [];  %all track location in 2D
allang_loc = [];
alldeltaC = [];  % finial-initial concentration

for c = 1:length(Tracks) %length(Paths)
    
    if isempty(find(cand==c))==0
    %%% process path to angles and distances
%     temp = Paths{c};  %for loading deleted tracks
    temp = Tracks(c);%.Path;  %for loading saved tracks
    temp1 = zeros(round(size(temp.Path,1)/1),2);
    temp1(:,1) = smooth(temp.SmoothX, filt,'sgolay',poly_degree); %smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp.SmoothY, filt,'sgolay',poly_degree); %smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    vecs = diff(subs);
    vecs = [vecs; vecs(end,:)];   % compensating for the size change/time shift after taking difference 
    angs = zeros(1,size(vecs,1));    
    ang_loc = zeros(1,size(vecs,1));

    %%% for time step calculations
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    trials = ones(1,size(vecs,1));
    xys = zeros(2,size(vecs,1));
    
    %%% iterate through worms
    for dd = (l_window+1):length(angs)
        %%% angle function
%         angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        ang_loc(dd) = angles([1,0],vecs(dd,:)/norm(vecs(dd,:)));
        
        %%% perpendicular concentration change
        perp_dir = [-vecs(dd-0,2), vecs(dd-0,1)];
        perp_dir = perp_dir/norm(perp_dir)*perp_dist;
%         dCp(dd) = Fcon(subs(dd-l_window,1)+perp_dir(1)*1., subs(dd-l_window,2)+perp_dir(2)*1)...
%              - Fcon(subs(dd-l_window,1)-perp_dir(1)*1, subs(dd-l_window,2)-perp_dir(2)*1);
%         dCp(dd) = (M(floor(subs(dd-l_window,2)+perp_dir(2)), floor(subs(dd-l_window,1)+perp_dir(1)))...
%                  - M(floor(subs(dd-l_window,2)-perp_dir(2)), floor(subs(dd-l_window,1)-perp_dir(1)))) / 1;  %computing normal direction dcp
        [xil, yil] = plate_bound(M, subs(dd-0,1)+perp_dir(1), subs(dd-0,2)+perp_dir(2));
        [xir, yir] = plate_bound(M, subs(dd-0,1)-perp_dir(1), subs(dd-0,2)-perp_dir(2));
        dCp(dd) = (M(yil,xil) - M(yir,xir));
%         norm(vecs(dd,:))
        %%% forward concentration change
%         dCs(dd) = Fcon(subs(dd,1),subs(dd,2)) - Fcon(subs(dd-l_window,1),subs(dd-l_window,2));
%         dCs(dd) = -(M(floor(subs(dd-l_window,2)), floor(subs(dd-l_window,1))) - M(floor(subs(dd,2)), floor(subs(dd,1)))) / 1;
        dCs(dd) = M(floor(subs(dd-0,2)), floor(subs(dd-0,1)));
%         (bin*fr2sec)
%         dCs(dd) = (log(M(floor(subs(dd-l_window,2)), floor(subs(dd-l_window,1)))) - log(M(floor(subs(dd,2)), floor(subs(dd,1)))))*fr2sec;
        
        %%% check displacement
        dds(dd) = norm(subs(dd,1)-subs(dd-l_window,1), subs(dd,2)-subs(dd-l_window,2));
        %norm(vecs(dd,:));
        
        %%% record location in 2D\
        xys(:,dd) = subs(dd-0,:)';
%         dd
%         dds(dd)
%         pause();

    end
    %remove zeros
    dCp = dCp((l_window+1):end);
    dCs = dCs((l_window+1):end);
    dds = dds((l_window+1):end);
    angs = angs((l_window+1):end);
    ang_loc = ang_loc(l_window+1:end);
    trials = trials((l_window+1):end);
    xys = xys(:, (l_window+1):end);
    
    allas = [allas angs];
    allang_loc = [allang_loc ang_loc];
    alldcp = [alldcp dCp];
    alldC = [alldC dCs];
    alldis = [alldis dds];  %displacement (effective velocity)
    trials(1:35) = NaN; trials(end-35:end) = NaN;
    alltrials = [alltrials trials];
    allxys = [allxys xys];
    
    alldeltaC = [alldeltaC (dCs(end)-dCs(2))];
    
    end
end

pos = find(allas==0 | alldC==0 | alldcp==0 | alldis==0);% | isnan(alltrials));
allas(pos) = [];
alldC(pos) = [];
alldcp(pos) = [];
alldis(pos) = [];
alltrials(pos) = [];
allxys(:,pos) = [];

%%  test with mask
pos = find(alldis<0.5);
alltrials(pos) = NaN;

%%
%     K_ = THETA(1);     % variance of von Mises
%     A_ = THETA(2);     % max turn probability
%     B_ = THETA(3:8);   % kernel for dC transitional concentration difference (weights on kerenl basis)
%     C_ = THETA(9);     % baseline turning rate
%     Amp = THETA(10);   % amplitude for kernel for dcp normal concentration difference
%     tau = THETA(11);   % time scale for dcp kernel
%     K2_ = THETA(12);   % vairance of the sharp turn von Mises
%     gamma = THETA(13); % uniform angle weight

%% test with stats-model for kernels
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
ang_fit = allas(1:100000);
dcp_fit = alldcp(1:100000);
ddc_fit = (alldC(1:100000));
trials_fit = alltrials(1:100000);
lfun = @(x)nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
% lfun = @(x)nLL_randomwalk(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);  % negative control
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.1    -inf -5];
UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1    inf 5];
% prs0 = rand(1,10);
prs0 = [50, 0.5, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.    0 0]; 
prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
x_MLE = x
fval
% K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7)/1; Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); sb = x(12); K2_ = x(13);  gamma = x(14);

%% test with stats-model for kernels, adding history of WV angles (for 54.5 strain with biased turns)
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
ang_fit = allas(1:100000);
dcp_fit = alldcp(1:100000);
ddc_fit = (alldC(1:100000));
trials_fit = alltrials(1:100000);
lfun = @(x)nLL_kernel_hist3(x, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.1    -inf 2  -10 .1];
UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1    inf 4  10 200];
% prs0 = rand(1,10);
prs0 = [50, 0.5, randn(1,nB)*1, 0.01, -1, 25, 1, 25, 5, 1.    -5 3 .1 1]; 
prs0 = prs0 + prs0.*randn(1,length(UB))*0.1;
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
x_MLE = x
fval

% K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7)/1; Amp = x(8); tau = x(9);
% Amp_h = x(10); tau_h = x(11); K2_ = x(13);  gamma = x(13); base_dc=x(14);
% base_dcp=x(15); Amp_h_wv = x(16); tau_h_wv = x(17);

%% test parameter distribution
% repeats = 10;
% prs_rep = zeros(length(prs0), repeats);
% for rr =1:repeats
%     % load observation
%     ang_fit = allas(1:end-1);
%     dcp_fit = alldcp(1:end-1);
%     ddc_fit = (alldC(1:end-1));
%     lfun = @(x)nLL_kernel_hist(x,ang_fit, dcp_fit, ddc_fit, cosBasis);
%     % fitting
%     prs0 = [1, 0.1, randn(1,nB), 0.01, -1, 100, -1, 100, 1];  %initialization 
%     prs0 = prs0+prs0.*randn(1,length(prs0))*.5;  %perturbation
%     [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[]); 
%     
%     prs_rep(:,rr) = x;
% end
% 
% figure()
% plot(prs_rep(:,:)','-o')

%% sufficient statistics
K_ = x(1); A_ = x(2)*1; B_ = x(3:6); C_ = x(7); Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); K2_ = x(12);  gamma = x(13); base_dc = x(14); base_dcp = x(15); 
% Amp_h_wv=x(16); tau_h_wv=x(17);
% base_dc = 0; base_dcp = 0;

figure
subplot(2,2,1)
xx = 0:length(cosBasis)-1;
yyaxis left; plot(Amp*exp(-xx/tau)); hold on
yyaxis right; plot(Amp_h_wv*exp(-xx/tau_h_wv))
title('\delta C^{\perp}, \delta \theta kernel')
subplot(2,2,2)
yyaxis left; plot(B_ * cosBasis'); hold on
yyaxis right; plot(Amp_h*exp(-xx/tau_h))
title('\delta C, \delta \theta kernel')

subplot(2,2,3)
K_dcp_rec = Amp*exp(-xx/tau);
filt_dcp = conv_kernel(dcp_fit, K_dcp_rec);%conv(dcp_fit, fliplr(Amp*exp(-xx/tau)), 'same');
[aa,bb] = hist((ang_fit - filt_dcp - base_dcp*1)*pi/180 , 1000);
bar( bb, 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb )) , 30); hold on
bar( bb, 1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi) , 30,'r');
title('von Mises for \delta C^{\perp}')
subplot(2,2,4)
K_dc_rec = B_*cosBasis';
filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
xx_h = 1:length(xx)*1;
K_h_rec = Amp_h*exp(-xx_h/tau_h);
filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
dc_dth = filt_ddc + 1*filt_dth;
Pturns = (A_-C_)./ (1 + exp( -(dc_dth + base_dc) )) + C_; %+sb
plot(dc_dth/length(K_dcp_rec) , Pturns,'o')
title('Logistic for \delta C')

%%% printing
disp(['K_=',num2str(K_),' K2_',num2str(K2_),'beta',num2str(sum(B_)),'alpha',num2str(Amp)])

%% joint density
n_brw = sum(Pturns);
n_wv = sum(1-Pturns);
p_z = n_brw + n_wv;
p_brw = n_brw/p_z;
p_wv = n_wv/p_z;

[aa,bb] = hist(allas*pi/180,500);
figure;
H = bar(bb,aa/sum(aa),30);
x_data = H.XData;  y_data = H.YData

[aa,bb] = hist(-(ang_fit - filt_dcp - base_dcp)*pi/180 , 500);
n_vm_wv = 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb )); 
n_vm_brw = 1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi);

figure
bar( bb, n_vm_wv/sum(n_vm_wv)*p_wv , 3); hold on
bar( bb, n_vm_brw/sum(n_vm_brw)*p_brw , 3,'r');
plot(x_data, y_data, '--')

%% autocorr check
[acf,lags] = autocorr(abs(allas),100); %(allang_loc,100);
f = fit(lags',acf','exp1');
disp(['ACF_tau= ',num2str(-1/f.b),', kernel_tau= ',num2str(tau)])
figure; plot(lags, acf)

%% Hessian computation
plfun = @(x)-nLL_kernel_hist2(x, ang_fit, dcp_fit, ddc_fit, cosBasis, 1, trials_fit);
[H, g] = compHess(plfun, x, 0.1);
[V,D] = eig(H);
figure; plot(diag(real(D)),'-o')

%% Generative model, with inferred parameters~
% reconstruct parameters '
% kappa = ((1/K_)^2*(180/pi)^1)^0.5;
% kappa2 = ((1/K2_)^2*(180/pi)^1)^0.5;
kappa = ((K_)^1)^1;
kappa2 = ((K2_)^1)^1;
if kappa2<0.05; kappa2=0.1; end  %simplify calculation
A = A_*1;
Kdcp = Amp*exp(-xx/tau)*14/5;
Kdth = Amp_h*exp(-xx_h/tau_h)*14/5;
Kddc = B_ * cosBasis'*14/5;
Kdth_wv = Amp_h_wv*exp(-xx_h/tau_h_wv);%.1*exp(-xx_h/50);%
Pturn_base = C_;
wind = size(cosBasis,1);

origin = [size(Cmap,2)/2,size(Cmap,1)/2];
alldths = [];
alldCs = [];
alldCps = [];
allths = [];

REP = 20;
T = floor(20*60*14/5);%2000;
dt = 1.;
% trackss = zeros(REP,T,2);
trackss = {};

%% simulations
clear simdata
simdata(length(REP)) = struct();

%%% gather time series
simdata(1).dth = []; %angle change
simdata(1).dcp = [];  %perpendicular dC
simdata(1).dc = [];  %tangential dC
simdata(1).dis = [];  %displacements
simdata(1).mask = [];  %mark for tracks
simdata(1).xy = [];  %all track location in 2D
simdata(1).theta = [];  %local angles

figure;
imagesc((M)); hold on;
% M = fliplr(M);
[Mx, My] = size(M);

for rep = 1:REP
    rep
vm = 0.2*32*(bin/14);  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
vs = .2;
tracks = zeros(T,2);
%%% from middle
tracks(1,:) = [size(M,2)*1/2 + randn()*300, randn()*300 + size(M,1)*1/2]*1. + 0.*[size(M,2)*rand(), size(M,1)*5/8];%origin; %initial position
%%% from random track
% temp = Tracks(randi(length(Tracks))).Path; tracks(1,:) = temp(1,:);  %data intials
%%% from collected structure
% temp = Data(rep).xy'; tracks(1,:) = temp(1,:);  % initial data point

tracks(2,:) = tracks(1,:)+randn(1,2)*vm*dt;%origin+randn(1,2)*vm*dt;
ths = zeros(1,T);  ths(1:3) = randn(1,3)*360; %initial angle
dxy = randn(1,2);  %change in each step

%"recordings"
dths = zeros(1,T); 
dcs = zeros(1,T);
dcps = zeros(1,T);

%placer for filtering vectors
dCv = zeros(1,wind);
dCpv = zeros(1,wind);
dthv = zeros(1,length(xx_h));
for t = 2:T
 
    %%% dC in a window
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) - Fcon(tracks(t-2,1),tracks(t-2,2)) ,  dCv];  % step-wise derivative
%     dCv = [ (M(floor(tracks(t-1,2)),floor(tracks(t-1,1))) - M(floor(tracks(t-2,2)), floor(tracks(t-2,1)))) / (1) ,  dCv];
    [xi, yi] = plate_bound(M, tracks(t-1,1), tracks(t-1,2));
    dCv = [ M(yi, xi) , dCv];
%     dCv = [M(floor(tracks(t-1,2)),floor(tracks(t-1,1))) - log(M(floor(tracks(t-2,2)), floor(tracks(t-2,1)))) ,  dCv];
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) ,  dCv];  % absolute concentration
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) - dCv(end),   dCv];  % moving threshold
    dCv = dCv(1:end-1);  %remove out of window
    perp_dir = [-dxy(2) dxy(1)];
    perp_dir = perp_dir/norm(perp_dir)*perp_dist;
    %%% dCp in a window
%     dCpv = [Fcon(tracks(t-1,1)+perp_dir(1),tracks(t-1,2)+perp_dir(2)) - Fcon(tracks(t-1,1)-perp_dir(1),tracks(t-1,2)-perp_dir(2)) ,  dCpv];
    [xil, yil] = plate_bound(M, tracks(t-1,1)+perp_dir(1), tracks(t-1,2)+perp_dir(2));
    [xir, yir] = plate_bound(M, tracks(t-1,1)-perp_dir(1), tracks(t-1,2)-perp_dir(2));
    dCpv = [(M(yil,xil) - M(yir,xir))/ 1 ,  dCpv];
    dCpv = dCpv(1:end-1);
    
    wv = (1*sum(Kdcp.*dCpv) + base_dcp*1 + sum(Kdth_wv.*(dthv))*1) + (vmrand(0,kappa))*180/pi;%kappa^1*randn;%length(wind)
    P_event = (A - Pturn_base) / (1+exp( -(sum(Kddc.*dCv)*1. + (sum(Kdth.*abs(dthv)*1)) *dt + 1*base_dc)+0) ) + Pturn_base;%length(wind)
    if rand < P_event*1
        beta = 1;
    else
        beta = 0;
    end
    
    if rand < gamma
        rt = beta*(vmrand(pi,1*kappa2)*180/pi);
    else
        rt = beta*(vmrand(0,0)*180/pi);
    end
%     rt = beta*[(randn*kappa2-180)*gamma + (1-gamma)*(rand*360-180)];
%     rt = beta*[(vmrand(pi,kappa2)*180/pi)*(gamma) + (1-gamma)*(vmrand(0,0)*180/pi)];
    dth = wv+rt;
%     if dth>180; dth = dth-180; end;  if dth<-180; dth = dth+360; end  %within -180~180 degree range
%     dth = sign(dth)*mod(abs(dth),180);  %periodic boundary
    dth = wrapToPi(dth*pi/180)*180/pi;
    
    %%% angular or turning history
    dthv = [dth*pi/pi , dthv]; %[dth  beta];%
    dthv = dthv(1:end-1);
    
    dths(t) = dth;
    dcs(t) = dCv(1);
    dcps(t) = dCpv(1);
    
%     vv = vm+vs*randn;
    vv = alldis(randi(length(alldis)));
    if vv<1
        vv = 1;
    end
    ths(t) = ths(t-1)+dth*dt;
%     if ths(t)>180; ths(t) = ths(t)-180; end;  if ths(t)<-180; ths(t) = ths(t)+360; end  %within -180~180 degree range
    e1 = [1,0];
    vec = [tracks(t-1,1)  tracks(t-1,2)]-origin; %current vector
    theta = acosd(max(-1,min((vec*e1')/norm(vec)/norm(e1),1)));  %current angle
    dd = [vv*sin(ths(t)*pi/180) vv*cos(ths(t)*pi/180)];
    R = [cos(theta*pi/180) sin(theta*pi/180); -sin(theta*pi/180) cos(theta*pi/180)];
    dxy = dd';%(R)*dd';

    tracks(t,1) = tracks(t-1,1)+dxy(1)*dt;
    tracks(t,2) = tracks(t-1,2)+dxy(2)*dt;
    
    %%% within the odor environment
    if tracks(t,1)>size(M,2)-1 | tracks(t,1)<3 | tracks(t,2)>size(M,1)-1 | tracks(t,2)<3
        tracks_ = zeros(t,2);
        tracks_(:,1) = tracks(1:t,1);
        tracks_(:,2) = tracks(1:t,2);  %ending track
        tracks = tracks_;
        break;
    end
end

% plot(tracks(:,1),tracks(:,2))
plot(tracks(:,1),tracks(:,2),'Color',[0 0 0]+0.8)
hold on
plot(tracks(1,1), tracks(1,2),'g.', 'MarkerSize',15) %(origin(1),origin(2),'ko')%
plot(tracks(end,1), tracks(end,2),'r.', 'MarkerSize',15)

alldths = [alldths dths];
alldCs = [alldCs dcs];
alldCps = [alldCps dcps];
allths = [allths ths];

% store as structure
    simdata(rep).dth = dths; %angle change
    simdata(rep).dcp = dcps;  %perpendicular dC
    simdata(rep).dc = dcs;  %tangential dC
    simdata(rep).xy = tracks';  %all track location in 2D
    simdata(rep).theta = ths;  %local angles

% trackss(rep,:,:) = tracks;
trackss{rep} = tracks;

end

set(gcf,'color','w'); %set( gca, 'xdir', 'reverse' )

pos = find(alldths==0 | alldCs==0 | alldCps==0 | allths==0);
alldths(pos) = [];
alldCs(pos) = [];
alldCps(pos) = [];
allths(pos) = [];

%% logP distribution
%%%
% using the mixture of VM model for navigation, we simulated with MLE paramerer and check the distribution of predictive logP
%%%
LLs = zeros(1,REP);  % tracks used for fitting
for ii = 1:REP
    LLs(ii) = nLL_kernel_hist2(x, simdata(ii).dth, simdata(ii).dcp, simdata(ii).dc, cosBasis, 1., ones(1,length(simdata(ii).dth))) / length(simdata(ii).dth);
end
normV = [(LLs-min(LLs))./(max(LLs)-min(LLs))]';
% blue to red. 
C = [normV normV normV];%[normV zeros(size(normV)) 1-normV];

figure();
imagesc(M); hold on;
for ii = 1:length(LLs) %10:60 %200:300
    plot(simdata(ii).xy(1,:),simdata(ii).xy(2,:), 'Color', C(ii,:)); 
    hold on
%     plot(Data_fit(ii).xy(1,1),Data_fit(ii).xy(2,1),'g.', 'MarkerSize',15)
%     plot(Data_fit(ii).xy(1,end),Data_fit(ii).xy(2,end),'r.', 'MarkerSize',15)
end

figure;
hist(LLs)

%% Check statistics
figure()
subplot(3,2,1); hist(allas,100);xlim([-180,180]);  ylabel('\delta \theta'); title('data'); subplot(3,2,2); hist(alldths,100); xlim([-180,180]);title('model');
subplot(3,2,3); hist(alldC,100);  ylabel('\delta C');      subplot(3,2,4); hist(alldCs,100)
subplot(3,2,5); hist(real((alldcp)),100); ylabel('\delta C^{\perp}');     subplot(3,2,6); hist(real((alldCps)),100)

%% summary stats
[a_data,b_data] = hist(allas,100);
norm_a_data = a_data/sum(a_data);
rad_b_data = b_data*pi/180;
figure;
bar(rad_b_data, norm_a_data);
hold on
% [aa,bb] = hist(-(ang_fit - filt_dcp - base_dcp)*pi/180 , 1000);
bb = rad_b_data;
K_ = 8;
K2_ = 2;
gamma = 0.1;
p_wv = 1/(2*pi*besseli(0,K_^2)) * exp(K_^2*cos( bb ));
p_rt= 1/(2*pi*besseli(0,K2_^2)) * exp(K2_^2*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi);
tr_rate = sum(Pturns>0.999)/length(Pturns);
pp = p_wv + p_rt*tr_rate*0.6;
bar(bb, pp/sum(pp),'FaceAlpha',0.5);
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('d\theta'); ylabel('P(d\theta)')

%%
figure();
[dcp_data_a,dcp_data_b] = hist(real((alldC(1:100000))),30);
[dcp_model_a,dcp_model_b] = hist(real((alldCs(:))),30);
bar(dcp_data_b, dcp_data_a/sum(dcp_data_a));
hold on
bar(dcp_model_b, dcp_model_a/sum(dcp_model_a),'FaceAlpha',0.5);
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('dC^{\perp}'); ylabel('P(dC^{\perp})')

%%
figure()
for ii =1:REP
    plot(trackss{ii}(:,1)-trackss{ii}(1,1),trackss{ii}(:,2)-trackss{ii}(1,2));
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
% figure()
% plot(bis,avang,'-o')
xlabel('\delta C^{\perp}', 'Fontsize',20)
ylabel('<\delta \theta>','Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);

%% 
figure;
dc_dth = (filt_ddc + 1*filt_dth);
Pturns = A_ ./ (1 + exp( (dc_dth))) + C_;
plot(dc_dth/length((B_ * cosBasis')) , Pturns,'o')

xlabel('filtered \delta C', 'Fontsize',20)
ylabel('P(turn|\delta C)', 'Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);

%%
figure;
Pturns = A_ ./ (1 + exp( (filt_ddc + 0*filt_dth))) + C_;
plot(filt_ddc/length(Kddc) , Pturns,'o')

xlabel('filtered \delta \theta', 'Fontsize',20)
ylabel('P(turn|\delta \theta)', 'Fontsize',20)
set(gcf,'color','w');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
xAX = get(gca,'YAxis');
set(xAX,'FontSize', 15);
