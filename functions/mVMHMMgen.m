% mVM-HM generative model
% load the fitted mmhat, with transition matrix and the weight parameters
% generative chemotaxis trajectories on the odor landscape M

%% generative from fitted GLM-HMM
A_inf = mmhat.A;  % hidden Markov transition
w_inf = mmhat.wts;  % inferred state eimission parameters
nstates = size(w_inf,3);  % number of states
basis = mmhat.basis;  % basis function for the kernels
wind = size(basis,1);  % time window for kernel operation

%% initialize vectors
lt = 10000;  % length of simulation
mc = dtmc(A_inf);
X = simulate(mc, lt);  % Markov chain for hidden states
REP = 50;  % repeating simulation tracks

vm = 0.2*32*(bin/14);  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
vs = .2;

alldths = [];
alldCs = [];
alldCps = [];
allths = [];
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
% tracks(1,:) = [size(M,2)*1/2 + randn()*200, randn()*200 + size(M,1)*1/2]*0. + 1.*[size(M,2)*rand(), size(M,1)*5/8];%origin; %initial position
temp = Tracks(rep).Path; tracks(1,:) = temp(1,:);  %data intials
tracks(2,:) = tracks(1,:)+randn(1,2)*vm*dt;%origin+randn(1,2)*vm*dt;
ths = zeros(1,T);  ths(1:3) = randn(1,3)*360; %initial angle
dxy = randn(1,2);  %change in each step
origin = [size(Cmap,2)/2,size(Cmap,1)/2];

%"recordings"
dths = zeros(1,T); 
dcs = zeros(1,T);
dcps = zeros(1,T);

%placer for filtering vectors
dCv = zeros(1,wind);
dCpv = zeros(1,wind);
dthv = zeros(1,length(xx_h));

%latent state dynamics
kt = simulate(mc, lt);  % integer states
% kt(kt==1)=0; kt(kt==2)=1; kt(kt==0)=2; %swapping state numbers
% kt(1:5000) = 2; %test for different initial conditions!

for t = 2:T
    
    %%% unwrap state parameter
    [Khat] = wts2params(mmhat.wts(:,:,floor(kt(t))), mmhat.basis');
%     ('K_dcp',K_dcp,'K_h',K_h,'K_dc',K_dc,'gamma',gamma,'A',A,'B',B,'kappa1',K_,'kappa2',K2_)
    Kdcp = Khat.K_dcp;  Kddc = Khat.K_dc; Kdth = Khat.K_h; gamma = Khat.gamma; kappa = Khat.kappa1; kappa2 = Khat.kappa2;
    A = Khat.A*1; B = Khat.B;
    base_dcp = 0;  base_dc = 0;
    
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
    
    wv = (1*sum(Kdcp.*dCpv) + base_dcp*0) + (vmrand(0,kappa))*180/pi;%kappa^1*randn;%length(wind)
    P_event = A / (1+exp( -(sum(Kddc.*dCv)*1. + 1*(sum(Kdth.*abs(dthv)*180/pi)) *dt + 0*base_dc)+0) ) + 0;%Pturn_base;%length(wind)
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

    dth = wv+rt;
    dth = wrapToPi(dth*pi/180)*180/pi;
    
    %%% angular or turning history
    dthv = [dth*pi/180 , dthv]; %[dth  beta];%
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

plot(tracks(:,1),tracks(:,2),'Color',[0 0 0]+0.8)
hold on
plot(tracks(1,1), tracks(1,2),'g.', 'MarkerSize',15) %(origin(1),origin(2),'ko')%
plot(tracks(end,1), tracks(end,2),'r.', 'MarkerSize',15)

alldths = [alldths dths];
alldCs = [alldCs dcs];
alldCps = [alldCps dcps];
allths = [allths ths];

trackss{rep} = tracks;

end

set(gcf,'color','w'); set( gca, 'xdir', 'reverse' )

pos = find(alldths==0 | alldCs==0 | alldCps==0 | allths==0);
alldths(pos) = [];
alldCs(pos) = [];
alldCps(pos) = [];
allths(pos) = [];

%% wts2params
function [Khat] = wts2params(x,cosBasis)
%%%
% Given the weight vector from HMM inference, unpack it back to the
% mixtue-VM emission model parameters for the generative process
%%%

    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11);
    xx = 1:length(cosBasis);
    K_dcp = Amp*exp(-xx/tau);
    K_h = Amp_h*exp(-xx/tau_h);
    K_dc = B_ * cosBasis;
    
    A = 1; B = 0;  % use sigmoid for now but change if we also fit this!
    
    Khat = struct('K_dcp',K_dcp,'K_h',K_h,'K_dc',K_dc,'gamma',gamma,'A',A,'B',B,'kappa1',K_,'kappa2',K2_);
end
