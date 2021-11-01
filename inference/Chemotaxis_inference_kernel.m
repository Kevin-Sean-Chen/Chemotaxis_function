%% mGLM fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load tracks and maps
%%% given Paths as cells with all 2D-trajectories
%%% and Fcon as a function that maps 2D position to an odor value readout

%%
figure;
%%% pre-processing parameters
poly_degree = 3;  %polynomial fit for the tracks
bin = 7;  %temporal binning
filt = 5;  %filtering tracks
l_window = 5;  %lag time

%%% gather time series
allas = []; %angles
alldcp = [];  %perpendicular dC
alldC = [];  %tangential dC
alldis = [];  %displacements
for c = 1:length(Paths)
    
    %%% process path to angles and distances
    temp = Paths{c};
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
        
        %forward concentration change
        dCs(dd) = Fcon(subs(dd-l_window,1),subs(dd-l_window,2));
        
        %%%check displacement
        dds(dd) = norm(vecs(dd,:));

    %remove zeros
    dCp = dCp(l_window+1:end);
    dCs = dCs(l_window+1:end);
    dds = dds(l_window+1:end);
    angs = angs(l_window+1:end);
    
    allas = [allas angs];
    alldcp = [alldcp dCp];
    alldC = [alldC dCs];
    alldis = [alldis dds];  %displacement (effective velocity)
    
    end
    
    
end

%% test with stats-model for kernels
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(6, [0, 10], 1.3);
ang_fit = allas(1:30000-1);
dcp_fit = alldcp(1:30000-1);
ddc_fit = (alldC(1:30000-1));
lfun = @(x)nLL_kernel_chemotaxis(x,ang_fit, dcp_fit, ddc_fit, cosBasis);
% [x,fval] = fminunc(f,randn(1,10));  %random initiation
[x,fval] = fminunc(lfun,[100, 0.001, randn(1,6), -1, 50]+randn(1,10)*0.);  %a closer to a reasonable value

% opts = optimset('display','iter');
% LB = [1e-5, 1e-5, ones(1,6)*-inf, -inf, 1e-5];
% UB = [1000, 1, ones(1,6)*inf, inf, 100];
% prs0 = rand(1,10);
% fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

x
fval

%% sufficient statistics
K_ = x(1); A_ = x(2); B_ = x(3:8); Amp = x(9); tau = x(10);
figure
subplot(2,2,1)
xx = 1:length(cosBasis);
plot(Amp*exp(-xx/tau))
title('\delta C_p kernel')
subplot(2,2,2)
plot(B_ * cosBasis')
title('\delta C kereel')
subplot(2,2,3)
filt_dcp = conv(dcp_fit, Amp*exp(-xx/tau), 'same');
[aa,bb] = hist((ang_fit - filt_dcp)*pi/180 , 300);
bar( bb, 1/(2*pi*besseli(0,K_)) * exp(K_*cos( bb )) , 100);
title('von Mises for \delta C_p')
subplot(2,2,4)
filt_ddc = conv( ddc_fit, B_, 'same' );
plot(filt_ddc , A_ ./ (1 + exp( filt_ddc)) ,'o')
title('Logistic for \delta C')

%% Generative model, with inferred parameters~
figure;
imagesc(Cmap); hold on;

% reconstruct parameters 
kappa = (1/K_)^0.5*(180/pi);  
A = A_;
Kddc = Amp*exp(-xx/tau);
Kdcp = B_ * cosBasis';
wind = size(cosBasis,1);

origin = [size(Cmap,2)/2,size(Cmap,1)/2];
alldths = [];
alldCs = [];
alldCps = [];
allths = [];

for rep = 1:30
    
T = 2000;
dt = 1;
vm = 2.5;  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
vs = 1.;
tracks = zeros(T,2);
tracks(1,:) = origin; %initial position
tracks(2,:) = origin+randn(1,2)*vm*dt;
ths = zeros(1,T);  ths(1:3) = randn(1,3)*360; %initial angle
dxy = randn(1,2);  %change in each step

dths = zeros(1,T); %"recordings"
dcs = zeros(1,T);
dcps = zeros(1,T);
dCv = zeros(1,wind);
dCpv = zeros(1,wind);
for t = 1+2:T
    %%% dC in a window
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) - Fcon(tracks(t-2,1),tracks(t-2,2)) ,  dCv];  % step-wise derivative
    dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) ,  dCv];  % absolute concentration
%     dCv = [Fcon(tracks(t-1,1),tracks(t-1,2)) - dCv(end),   dCv];  % moving threshold
    dCv = dCv(1:end-1);  %remove out of window
    perp_dir = [-dxy(2) dxy(1)];
    perp_dir = perp_dir/norm(perp_dir);
    %%% dCp in a window
    dCpv = [Fcon(tracks(t-1,1)+perp_dir(1),tracks(t-1,2)+perp_dir(2)) - Fcon(tracks(t-1,1)-perp_dir(1),tracks(t-1,2)-perp_dir(2)) ,  dCpv];
    dCpv = dCpv(1:end-1);
    
    wv = -sum(Kdcp.*dCpv)/length(wind) + kappa*randn;
    P_event = A/(1+exp( (sum(Kddc.*dCv)/length(wind)) *dt));
    if rand < P_event
        beta = 1;
    else
        beta = 0;
    end
    rt = beta*(rand*360-180);
    dth = wv+rt;
%     if dth>180; dth = dth-180; end;  if dth<-180; dth = dth+360; end  %within -180~180 degree range
    
%     dths(t) = dth;
%     dcs(t) = dC;
%     dcps(t) = dCp;
    
    vv = vm+vs*randn;
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
plot(origin(1),origin(2),'ko')%(target2(1)/2,target2(2),'ko')

% alldths = [alldths dths];
% alldCs = [alldCs dcs];
% alldCps = [alldCps dCp];
% allths = [allths ths];

end
