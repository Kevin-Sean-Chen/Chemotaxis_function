% mVMG-HM generative model
% load the fitted mmhat, with transition matrix and the weight parameters
% generative chemotaxis trajectories on the odor landscape M
% this is an updated version that includes velocity and input-driven
% transition probabilities

%% generative from fitted GLM-HMM
A_inf = mmhat.A;  % hidden Markov transition
w_inf = mmhat.wts;  % inferred state eimission parameters
w_state_inf = mmhat.wts_state;  % inferred state-transition kernels
nstates = size(w_inf,3);  % number of states
basis = mmhat.basis;  % basis function for the kernels
wind = size(basis,1);  % time window for kernel operation
betaT = 1;  % testing soft-max for now
Basis = mmhat.basis;
perp_dist = 1;

%% initialize vectors
lt = 10000;  % length of simulation
REP = 10;  % repeating simulation tracks

alldths = [];
alldCs = [];
alldCps = [];
allths = [];
allstate = [];
allvs = [];
trackss = {};
T = lt;
dt = 1.;

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
origin = [size(Cmap,2)/2,size(Cmap,1)/2];

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
    Tij = ones(nstates, nstates);  % state transition probability
    for ii = 1:nstates
        for jj = 1:nstates
            if ii == jj
                Tij(ii,jj) = 0 + A_inf(ii,jj);  
            else
                K_ij_dc = (squeeze(w_state_inf(ii,jj,:))'*Basis');  % reconstruct kernel
                Tij(ii,jj) = (0+exp(sum(K_ij_dc.*dCv)*betaT))*1 + 0*A_inf(ii,jj);
            end
        end
    end
    Tij = Tij ./ (sum(Tij,2) + 1e-8);  % normalize rows
    temp_state = [find(rand < cumsum(Tij(kt(t-1),:)))];  % Markov transition
    kt(t) = temp_state(1);
    
    %%% unwrap state parameter
    [Khat] = wts2params(mmhat.wts(:,:,floor(kt(t))), mmhat.basis');
    Kdcp = Khat.K_dcp;  Kddc = Khat.K_dc; Kdth = Khat.K_h; gamma = Khat.gamma; kappa = Khat.kappa1; kappa2 = Khat.kappa2;  % emissions
    A = Khat.A*1; B = Khat.B;  % sigmoid
    base_dcp = Khat.b_dcp;  base_dc = Khat.b_dc;  % bias
    base_dcp = 0; base_dc = 0;
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
    
    wv = (1*sum(Kdcp.*dCpv) + base_dcp*1) + (vmrand(0,kappa))*180/pi;%kappa^1*randn;%length(wind)
    P_event = (A-B) / (1+exp( -(sum(Kddc.*dCv)*1. + 1*(sum(Kdth.*abs(dthv)*180/pi)) *dt + 1*base_dc)+0) ) + B;%Pturn_base;%length(wind)
    if rand < P_event*14/14;
        beta = 1;
    else
        beta = 0;
    end
    
    if rand < gamma
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
    vv = gamrnd(k_, theta_);
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

trackss{rep} = tracks;

end

set(gcf,'color','w'); set( gca, 'xdir', 'reverse' )

pos = find(alldths==0 | alldCs==0 | alldCps==0 | allths==0);
alldths(pos) = [];
alldCs(pos) = [];
alldCps(pos) = [];
allths(pos) = [];
allstate(pos) = [];
allvs(pos) = [];

%% compare to data (if available)
nbins = 50;
figure
subplot(121)
hh = histogram(yyf(1,:), nbins, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
hh = histogram(alldths, nbins, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('d\theta'); set(gca,'Fontsize',20); set(gca, 'YScale', 'log');

temp_hist = yy(2,:)*1/30/(5/14);
temp_hist(temp_hist>.7) = [];
subplot(122)
hh1 = histogram(temp_hist, nbins, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
hh = histogram(allvs*1/30/(5/14), hh1.BinEdges, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('dr'); set(gcf,'color','w'); set(gca,'Fontsize',20);

%% analyze dwell time distribution
bins = 50;
rescale_t = 5/14;
[dwell_times] = compute_dwell_time(allstate, []);

max_dt = max(max([dwell_times{1} , dwell_times{2}]));
bin_edges = linspace(0,1,bins).*max_dt*rescale_t;
figure;
H2 = histogram(dwell_times{1}*rescale_t, bin_edges,'Normalization', 'pdf');hold on
H1 = histogram(dwell_times{2}*rescale_t, bin_edges,'Normalization', 'pdf'); 
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('pdf'); title('dwell time')

%% wts2params
function [Khat] = wts2params(x,cosBasis)
%%%
% Given the weight vector from HMM inference, unpack it back to the
% mixtue-VM emission model parameters for the generative process
%%%
    
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A=x(12); B=x(13); b_dc=x(16); b_dcp=x(17);
    k_ = x(14);  theta_ = x(15);
    xx = 1:length(cosBasis);
    K_dcp = Amp*exp(-xx/tau);
    K_h = Amp_h*exp(-xx/tau_h);
    K_dc = B_ * cosBasis;
    
    Khat = struct('K_dcp',K_dcp,'K_h',K_h,'K_dc',K_dc,'gamma',gamma,'A',A,'B',B,'kappa1',K_,'kappa2',K2_, 'b_dc',b_dc,'b_dcp',b_dcp,'k',k_,'theta',theta_);
end
