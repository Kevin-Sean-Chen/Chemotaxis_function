%% mGLM fitting with one-step back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% point-wise fitting without considering kernels
%%
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211013_biased_110mM/Landscape.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20211013_biased_110mM/OdorFx.mat');
Fcon = Fcon.F;
Cmap = Cmap.vq1;
Cmap = fliplr(flipud(Cmap));  %%% for inverted camper!
Paths = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Naive_tracks.mat');
Paths = Paths.Tracks;

%%
figure;
poly_degree = 3;  %polynomial fit for the tracks
bin = 7;  %temporal binning
filt = 7;  %filtering tracks
l_window = 7;  %lag time
allas = []; %angles
alldcp = [];  %perpendicular dC
alldC = [];  %tangential dC
alldis = [];  %displacements


for c = 1:length(Paths)
    
    %%% process path to angles and distances
    temp = Paths{c};
%     temp = del_cls(c).Path;
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
end

%% test with stats-model
f = @(x)nLL_chemotaxis(x,allas(1:10000-1),alldcp(1:10000-1),(alldC(1:10000-1)));
[x,fval] = fminunc(f,rand(1,4));%[0,5,0.1,0.1]);  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(f,[-1, 100, 0.0, 10]+0*rand(1,4));  %a closer to a reasonable value
x
fval

%% Generative model
figure;
imagesc(Cmap); hold on;
%dth = N(alpha*dcp,K) + (A/(1+exp(B*dC)))*U[0,2*pi];
alpha = x(1);  kappa = (1/x(2))^0.5*(180/pi);  A = x(3);  B = x(4);
% alpha = 1.1;  kappa = 0.1*(180/pi);  A = 0.4;  B = 20;  %some kind of model
% alpha = rand(1); kappa = rand(1); A = rand(1); B = rand(1);  %initial random control

origin = [size(Cmap,2)/2,size(Cmap,1)/2];

alldths = [];
alldCs = [];
alldCps = [];
allths = [];
for rep = 1:30
    
T = 3000;
dt = 5.5;
vm = .5;  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
vs = 0.5;
tracks = zeros(T,2);
tracks(1,:) = origin; %initial position
tracks(2,:) = origin+randn(1,2)*vm*dt;
ths = zeros(1,T);  ths(1) = randn(1); %initial angle
dxy = randn(1,2);  %change in each step

dths = zeros(1,T); %"recordings"
dcs = zeros(1,T);
dcps = zeros(1,T);
for t = 1+2:T
    dC = Fcon(tracks(t-1,1),tracks(t-1,2)) - Fcon(tracks(t-2,1),tracks(t-2,2));
    perp_dir = [-dxy(2) dxy(1)];
    perp_dir = perp_dir/norm(perp_dir);
    dCp = Fcon(tracks(t-1,1)+perp_dir(1),tracks(t-1,2)+perp_dir(2)) - Fcon(tracks(t-1,1)-perp_dir(1),tracks(t-1,2)-perp_dir(2));
    
    wv = -alpha*dCp + kappa*randn;
    P_event = A/(1+exp(B*dC*dt));
    if rand < P_event
        beta = 1;
    else
        beta = 0;
    end
    rt = beta*(rand*360-180);
    dth = wv+rt;
%     if dth>180; dth = dth-180; end;  if dth<-180; dth = dth+360; end  %within -180~180 degree range
    
    dths(t) = dth;
    dcs(t) = dC;
    dcps(t) = dCp;
    
%     vv = vm+vs*randn;
    vv = alldis(randi(length(alldis)));
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

alldths = [alldths dths];
alldCs = [alldCs dcs];
alldCps = [alldCps dcps];
allths = [allths ths];

end