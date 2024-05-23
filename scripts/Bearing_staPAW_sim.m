% Bearing_staPAW_sim
%%% given the fitted parameters, simulate trajectoriews and then compute bearing 
%%% this is different from the forward smoothing method as it take full simulations
%%% here we compare to different model predictions

%% loading data
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_staPAW.mat')
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt100_50_staPAW.mat')

%% load trained models
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat')
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240520_023653_cv_staPWA.mat') %%% new Data_salt0_50_
rng(1) %37 42 1

%% asign models
rr = 3; %3
dPAW_fit = all_record(rr,1,1).params;  % single state
staPAW_fit = all_record(rr,2,1).params;  % two-state model
temp_ws = dPAW_fit.wts;
temp_ws([2:5]) = zeros(1,4);
temp_ws([14:17]) = ones(1,4);
% null_fit = dPAW_fit;
% null_fit.wts = temp_ws;  % ad-hoc ablation!

%%% one at a time
model_choice = dPAW_fit;  % time step=5/14, tune wv strength, remove states, compute pirouette with window

% model_choice = staPAW_fit; % test with data in the same way (z-state)
% model_choice.wts_state = model_choice.wts_state*0;  % same setup as staPAW but not transition

%% environement
mmhat = model_choice; %all_record(rr, 2, 3).params;  % two-state model
% [rows, cols] = size(M);
rows=2500; cols=3000;
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
% [x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M = (y_*0+1) .* gradient_x;

%% generative from fitted GLM-HMM
A_inf = mmhat.A;  % hidden Markov transition
w_inf = mmhat.wts;  % inferred state eimission parameters
w_state_inf = mmhat.wts_state;  % inferred state-transition kernels
nstates = size(w_inf,3);  % number of states
basis = mmhat.basis;  % basis function for the kernels
wind = size(basis,1);  % time window for kernel operation
betaT = 1.;  % testing soft-max for now
Basis = mmhat.basis;
perp_dist = 1;

%% initialize vectors
lt = 3000;  % length of simulation
REP = 200;  % repeating simulation tracks

alldths = [];
alldCs = [];
alldCps = [];
allths = [];
allstate = [];
allvs = [];
allxy = [];
% trackss = {};
trackss = struct();
T = lt;
dt = 1;
t_step = 5/14;  % simulation steps  5, 14

%% chemotaxis dynamics
figure;
imagesc((M)); hold on;% M = fliplr(M);
[Mx, My] = size(M);

for rep = 1:REP
    rep

%initialize tracks
tracks = zeros(T,2);
tracks(1,:) = [size(M,2)*1/2 + randn()*200, randn()*200 + size(M,1)*1/2]*1. + 0.*[size(M,2)*rand()*1/2, size(M,1)*1/2];%origin; %initial position
% temp = Data(rep).xy'; tracks(1,:) = temp(1,:);  %data intials
tracks(2,:) = tracks(1,:)+randn(1,2)*dt;%origin+randn(1,2)*vm*dt;
ths = zeros(1,T);  ths(1:3) = randn(1,3)*360; %initial angle
dxy = randn(1,2);  %change in each step
origin = [size(M,2)/2, size(M,1)/2];

%"recordings"
dths = zeros(1,T); 
dcs = zeros(1,T);
dcps = zeros(1,T);
vs = zeros(1,T);

%placer for filtering vectors
dCv = zeros(1,wind);
dCpv = zeros(1,wind);
dthv = zeros(1,wind);

%latent state dynamics
kt = randi([1 2],1,T); 

for t = 2:T
    
    %%% state-transitions
%     Tij = ones(nstates, nstates);  % state transition probability
%     for ii = 1:nstates
%         for jj = 1:nstates
%             if ii == jj
%                 Tij(ii,jj) = 0 + A_inf(ii,jj);  
%             else
%                 K_ij_dc = (squeeze(w_state_inf(ii,jj,:))'*Basis');  % reconstruct kernel
%                 Tij(ii,jj) = (0+exp(sum(K_ij_dc.*dCv)*betaT))*1 + 0*A_inf(ii,jj);
%             end
%         end
%     end
%     Tij = Tij ./ (sum(Tij,2) + 1e-8);  % normalize rows
%     temp_state = [find(rand < cumsum(Tij(kt(t-1),:)))];  % Markov transition
%     if length(temp_state)==0
%         kt(t) = randi(numel([2,1]));
%     else
%         kt(t) = temp_state(1);
%     end
    kt(t) = 1;  %%% for dPAW w/o states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% unwrap state parameter
    [Khat] = wts2params(mmhat.wts(:,:,floor(kt(t))), mmhat.basis');
    Kdcp = Khat.K_dcp;  Kddc = Khat.K_dc; Kdth = Khat.K_h; gamma = Khat.gamma; kappa = Khat.kappa1; kappa2 = Khat.kappa2;  % emissions
    A = Khat.A*1; B = Khat.B;  % sigmoid
    base_dcp = Khat.b_dcp;  base_dc = Khat.b_dc;  % bias
%     base_dcp = 0; base_dc = 0;
    k_ = Khat.k;  theta_ = Khat.theta;  % velcoity
    
    %%% dC in a window
    [xi, yi] = plate_bound(M, tracks(t-1,1), tracks(t-1,2));
    dCv = [ M(yi, xi) , dCv];
    dCv = dCv(1:end-1);  %remove out of window
    perp_dir = [-dxy(2) dxy(1)];
    perp_dir = perp_dir/norm(perp_dir)*perp_dist;
    %%% dCp in a window
    [xil, yil] = plate_bound(M, tracks(t-1,1)+perp_dir(1), tracks(t-1,2)+perp_dir(2));
    [xir, yir] = plate_bound(M, tracks(t-1,1)-perp_dir(1), tracks(t-1,2)-perp_dir(2));
    dCpv = [(M(yil,xil) - M(yir,xir))/ 1 ,  dCpv];
    dCpv = dCpv(1:end-1);
    
    wv = (-1*sum(Kdcp.*dCpv) + base_dcp*0) + (vmrand(0,kappa))*180/pi;%kappa^1*randn;%length(wind)
    P_event = (A-B) / (1+exp( (sum(Kddc.*dCv)*-1. + 1*(sum(Kdth.*abs(dthv)*180/pi)) *dt + 0*base_dc)+0) ) + B;%Pturn_base;%length(wind)
    if rand < P_event*t_step;
        beta = 1;
        vv = gamrnd(k_(1), theta_(1))/1;
    else
        beta = 0;
        vv = gamrnd(k_(2), theta_(2))/1;
    end
    
    if rand < gamma*1;
        rt = beta*(vmrand(pi,1*kappa2)*180/pi);
    else
        rt = beta*(vmrand(0,0)*180/pi);
    end

    dth = wv+rt;
    dth = wrapToPi(dth*pi/180)*180/pi;
    
    %%% angular or turning history
    dthv = [dth*pi/180 , dthv]; %[dth  beta];%
    dthv = dthv(1:end-1);
    
    dths(t) = dth;
    dcs(t) = dCv(1);
    dcps(t) = dCpv(1);
    
    %%% displacement from gamma
    vs(t) = vv;
    
    ths(t) = ths(t-1)+dth*dt;
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

plot(tracks(:,1),tracks(:,2),'Color',[0 0 0]+0.2)
hold on
plot(tracks(1,1), tracks(1,2),'g.', 'MarkerSize',15) %(origin(1),origin(2),'ko')%
plot(tracks(end,1), tracks(end,2),'r.', 'MarkerSize',15)

alldths = [alldths dths(1:t)];
alldCs = [alldCs dcs(1:t)];
alldCps = [alldCps dcps(1:t)];
allths = [allths ths(1:t)];
allstate = [allstate kt((1:t))];
allvs = [allvs vs(1:t)];
allxy = [allxy tracks'];

% trackss{rep} = tracks;
trackss(rep).dth = dths(1:t);
trackss(rep).vs = vs(1:t);
trackss(rep).state = kt(1:t);
trackss(rep).dc = dcs(1:t);
trackss(rep).dcp = dcps(1:t);

end

set(gcf,'color','w'); set( gca, 'xdir', 'reverse' )

pos = find(alldths==0 | alldCs==0 | alldCps==0 | allths==0);
alldths(pos) = [];
alldCs(pos) = [];
alldCps(pos) = [];
allths(pos) = [];
allstate(pos) = [];
allvs(pos) = [];
allxy(:,pos) = [];

%% processing data
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:100000];
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(model_choice,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end
xys = allxys(:,wind_test);

target = [3000, 1250];%2500];  % gradient vector

%% make a mini generative process for xy position with the fitted staPAW... given states!
% xy_sim = xys*0;
% dth_sim = zeros(1,length(xys));
% wind = 50;  % should be larger than the kernel size
% xy_sim(:,wind) = xys(:,wind);
% 
% for tt = wind:length(xx)-1
%     [dth_pred, dr_pred] = model_predit_PAW(model_choice, xx(:,tt-wind+1:tt-0), yy(:,tt-wind+1:tt-0), gams_(:,tt-wind+1:tt-0));  % forward model
%     dr = dr_pred(end);
%     dth = dth_pred(end);
%     xy_sim(:,tt+1) = xys(:,tt) + ([dr*sin(dth*pi/180) dr*cos(dth*pi/180)])';  % make predictions about location (function of dr, dv) for one step
% %     xy_sim(:,tt+1) = xy_sim(:,tt) + ([dr*sin(dth*pi/180) dr*cos(dth*pi/180)])';  % without data
%     dth_sim(tt+1) = dth;
% end

%% replace with generative model
xy_sim = (allxy);
dth_sim = alldths;

%% replacing data back
% xy_sim = xys;
% dth_sim = yy(1,:);

%% defining states (change for dPAW comparison)
pre_t = 28; %8, 2, 28
post_t = 28;
state_vec = allstate*0;%gams_(1,:)*0; %gams_*0; %
%%% staPAW states!
% pos_state_stapaw = find(gams_(2,:)>0.8); %<gams_(2,:));
%%% mock states!
% pos_state_turn = find(abs(yy(1,:))>50);  % use this for fare comparison across models
% pos_state_turn = find(abs(dth_sim)>50);  %18
pos_state_turn = find(abs(allstate)>1.5);  %18

%%% logics
% pos_state = setdiff(pos_state_turn, pos_state_stapaw);  % only turns not states
% pos_state = setdiff(pos_state_stapaw, pos_state_turn);  % only state-switches not turns

% pos_state = pos_state_stapaw; 
pos_state = pos_state_turn;
fix_t = 10;

%% test for classic pirouette! for data or dPAW
windp = 10; %10, 4
pos_state_turn = find(abs(dth_sim)>50);
% pos_state_turn = find(abs(yy(1,:))>50); 
state_vec = allstate*0;
state_vec(pos_state_turn) = ones(1,length(pos_state_turn));
state_vec = conv(state_vec,ones(1,windp), 'same');
state_vec(find(state_vec>1)) = 1;

%% iterations
% state_vec(pos_state) = ones(1,length(pos_state)); %%% for dPAW
trans_pos = diff(state_vec);
trans12 = find(trans_pos>0)-0;
trans21 = find(trans_pos<0)+0;
npairs = min([length(trans12), length(trans21)]);
Bpairs = zeros(2,npairs);
vecs = diff(xy_sim,1,2);
for bb = 10:npairs-10
    %%% pre vector
    v1 = (target - xy_sim(:, trans12(bb))');  %target
%     v2 = vecs(:,trans12(bb)-pre_t)';
    v2 = -(xy_sim(:,trans12(bb)) - xy_sim(:,trans12(bb)-pre_t))';
    v1 = [3000 v2(2)];
    
    Bpre = angles(v1, v2);

    %%% post vector
    v1 = (target - xy_sim(:, trans21(bb))'); %target
%     v2 = vecs(:,trans21(bb)+post_t)';
    v2 = -(xy_sim(:,trans21(bb)) - xy_sim(:,trans21(bb)+post_t))';

%     v1 = -(target - xys(:, trans12(bb)+fix_t)'); %target
%     v2 = (xys(:,trans12(bb)+fix_t) - xys(:,trans12(bb)+fix_t+post_t))';  % another timed-control!
    v1 = [3000 v2(2)];
    
    Bpost = angles(v1, v2);
    
    %%% recording pairs
    Bpairs(:,bb) = [Bpre; Bpost];
end

%%
figure()
% hist(Bpairs(1,:),20, 'FaceColor', 'r'); hold on
% hist(Bpairs(2,:),20, 'FaceColor', 'b')
nbins = 12;
subplot(131)
histogram(Bpairs(1,:), nbins, 'FaceColor', 'r', 'EdgeColor', 'none', 'Normalization', 'probability');  % Use 20 bins for the histogram
hold on; xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20);

% subplot(132)
% histogram(Bpairs(2,:), nbins, 'FaceColor', 'r', 'EdgeColor', 'none', 'Normalization', 'probability');

subplot(133)
shuffle_rep = 150;
dB = -diff(Bpairs,1);
Bc  = [];
for ss = 1:shuffle_rep
    Bc = [Bc dB(randperm(length(dB))) + Bpairs(1,:)];
%     Bc = dB(randperm(length(dB))) + Bpairs(1,:);
end
Bc = wrapToPi(Bc*pi/180)*180/pi;
% histogram(Bc, nbins, 'FaceColor', 'k', 'EdgeColor', 'none', 'Normalization', 'probability'); hold on
H = histogram(Bc, nbins, 'Normalization', 'probability','Visible', 'off');
binEdges = H.BinEdges;  counts = H.BinCounts;
midb = (binEdges(2:end)+binEdges(1:end-1))/2;
plot(midb, counts/sum(counts)); hold on
histogram(Bpairs(2,:), nbins, 'FaceColor', 'b', 'EdgeColor', 'none', 'Normalization', 'probability'); hold on; xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20);

subplot(132)
dB = wrapToPi(dB*pi/180)*180/pi;
histogram(dB, nbins, 'FaceColor', 'g', 'EdgeColor', 'none', 'Normalization', 'probability'); xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20);

test = Bpairs(2,:);
length(find(abs(test)<90))/length(test)  % compute fraction that aligns

%%
%%% NOTES
% should be slightly different from the old analysis
% the bearing pre-pirouette might be accounted by individual turns
% show that staPAW recap this

%%
% % Estimate PDFs using kernel density estimation
% [f1, x1] = ksdensity(Bpost_data); % PDF of data 1
% [f2, x2] = ksdensity(Bpost_stapaw); % PDF of data 2  %Bpost_stapaw
% % Compute KL divergence
% kl_divergence = sum(f2 .* log(f2 ./ f1))

%% wts2params
function [Khat] = wts2params(x,cosBasis)
%%%
% Given the weight vector from HMM inference, unpack it back to the
% staPAW parameters for the generative process
%%%
    
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A=x(12); B=x(13); b_dc=x(18); b_dcp=x(19);
    k_ = x(14:15);  theta_ = x(16:17);
    xx = 1:length(cosBasis);
    K_dcp = Amp*exp(-xx/tau);
    K_h = Amp_h*exp(-xx/tau_h);
    K_dc = B_ * cosBasis;
    
    Khat = struct('K_dcp',(K_dcp),'K_h',K_h,'K_dc',(K_dc),'gamma',gamma,'A',A,'B',B,'kappa1',K_,'kappa2',K2_, 'b_dc',b_dc,'b_dcp',b_dcp,'k',k_,'theta',theta_);
end
