%Chemotaxis_dynamics
%021419
clear
clc
%%
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Frames','SmoothSpeed'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% % criteria %%%
nn = length(Tracks); %number of worms selected
mint = 60*1;%60*1; %minimum time in seconds
minx = 50*1;  %minimum displacement
disth = 500;  %radius of pixels from target
target = [2517,975];%[950,1100];%  %position of target/sourse of odorant (approximated from images)
endingt = 60*15;  %only taking the first few minutes

%visualize all paths and check criteria
cand = [];
figure;
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        if (Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2 > minx^2
            pos = find(Tracks(i).Time<endingt);
            if isempty(pos)~=1
                %plot(Tracks(i).Path(pos,1),Tracks(i).Path(pos,2)); %axis([1500,2500,0,2000]); pause();
                plot(Tracks(i).Path(:,1),Tracks(i).Path(:,2)); %pause();%axis([1500,2500,0,2000]); 
                hold on
                cand = [cand i];
            end
        end
    end
end

%% First-passage time measurement
crossing = 1800;
timing = [];
figure;
for c = 1:length(cand)
    pos = find(Tracks(cand(c)).Path(:,1)>crossing);
    if size(pos) ~= 0
        timing = [timing Tracks(cand(c)).Time(pos(1))];
    end
end
hist(timing)
%% Smoothed speed
speed = [];
for i = 1:nn
    speed = [speed Tracks(i).SmoothSpeed];
end

%% approximate diffusion coefficient through time
partt = 5;
Mt = 30*60;
bin_times = [0:Mt/(partt):Mt];

Ds = [];
for tt = 1:length(bin_times)-1
    sub_tracks = FilterTracksByTime(Tracks,bin_times(tt)*14,bin_times(tt+1)*14);  %Mochi's code to select tracks within a time window
    xsqu = 0;
    for ii = 1:length(sub_tracks)
        xsqu = xsqu + mean(diff(sub_tracks(ii).Path(:,1)).^2);
    end
    Ds = [Ds xsqu/length(sub_tracks)];
end

%% approximate diffusion coefficient through space
partt = 6;
Mt = 30*60;
%frame is 1944 X 2592
bin_space = [2592/2: (2592/2)/partt: 2592];

Ds = [];
for xx = 1:length(bin_space)-1
    xsqu = 0;
    tr = 0;
    for ii = 1:length(Tracks)
        pos = find(Tracks(ii).Path(:,1)>bin_space(xx) & Tracks(ii).Path(:,1)<bin_space(xx+1));
        if length(pos)>1
            subpath = Tracks(ii).Path(pos,1);
            xsqu = xsqu + mean(diff(subpath).^2);
            tr = tr+1;
        end
    end
    Ds = [Ds xsqu/tr];
end

%% %%%%% Peurette
%% Estimating concentration
figure;
bin = 5;
filt = 7;
p_thr = 60;
ct = 0;
allas = [];
alldC = [];

timelim = 60*10;

for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1),filt);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2),filt);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [[0,0]; diff(subs)];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    for dd = 1:length(dists)
        dists(dd) = distance(temp1(dd,:),target);
    end
    %if distance(subs(1,:),p2) > 300%disth  &&  distance(subs(1,:),p2) > disth 
    %if isempty(find(dists<disth)) ~= 1  %&&   min(newtime) < 600
        ct = ct +1;
    %if sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) < 300 && sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) > 300%dist
    %if newtime(1)<60*30  &&  distance(p1,p2) > 300  %newtime(1)> 60*5
    if 1==1
%     if length(find(temp1(:,2)>1000)) < length(find(temp1(:,2)<1000))  %approximately the gradient front~~
        

    %%%for angle
    dCs = zeros(1,size(vecs,1));  
    for dd = 2:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        %angs(dd) = angles(vecs(dd,:),target-subs(dd,:)); %ThetaInDegrees;%
        angs(dd) = angles(vecs(dd-1,:),vecs(dd,:));
        dCs(dd) = Est_con(temp1(dd-1,1),temp1(dd-1,2),target(1),target(2),50);
    end
    
    %%%for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(real(angs),p_thr,vecs);
%     dC = [];
%     for tt = 1:length(timestamps)
%         dC = [dC Est_con(temp1(timestamps(tt),1),temp1(timestamps(tt),2),target(1),target(2),100)];
%     end
    allas = [allas angs];
    alldC = [alldC dCs];
    
    end
    
    
end

plot(diff(alldC),allas(1:end-1),'o')

%% plot adaptive binning
bins = 30;
threshold = 60;
dC_ = diff(alldC);
dA_ = allas;
dA_(isnan(allas)) = 0;
dA_ = dA_(1:end-1);
pos = find(abs(dC_)<0.001);
dC_(pos) = [];  dA_(pos) = [];

% [va,id] = sort(dC_);
% avang = zeros(1,bins);
% dCbin = zeros(1,bins);
% numinbin = floor(length(id)/bins);
% for ii = 1:bins
%     %avang(ii) = mean(dA_(id((ii-1)*numinbin+1:ii*numinbin)));
%     avang(ii) = length(find(dA_(id((ii-1)*numinbin+1:ii*numinbin))>threshold));
%     dCbin(ii) = va(ii*numinbin);
% end
% plot(dCbin,avang,'-o')

avang = zeros(1,bins);
[cts,bis] = hist(dC_,bins);
for bi = 2:length(bis)
    pos = find(bis(bi-1)<dC_ & bis(bi)>dC_);
    if isempty(pos)~=1
        avang(bi) = length(dA_(pos)>threshold)/(length(pos)+1);%%mean(dA_(pos));%
    end
end
plot(bis(avang~=NaN),avang(avang~=NaN),'-o')
%plot(bis(avang~=0),avang(avang~=0),'-o')

% b = glmfit(dC_,dA_,'binomial','link','logit')
% plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)
%% Logistic
CC = bis;%dCbin; %bis(avang~=0);%
AA = avang%/max(avang); %avang(avang~=0);%
pos = CC<-0.01;
CC(pos) = [];
AA(pos) = [];

[b,dev,stats] = glmfit(CC, AA', 'binomial', 'link', 'logit')
xx = linspace(min(CC), max(CC), 30);
yfit = glmval(b,xx,'logit');
plot(xx,yfit,'-')
hold on
plot(CC,AA,'ro')

%% %%%%% Weathervaning
%% %%%%%
figure;
bin = 5;
filt = 7;
p_thr = 50;
ct = 0;
allas = [];
alldB = [];
alldcp = [];
alldC = [];

timelim = 60*10;

for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1),filt);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2),filt);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [[0,0]; diff(subs)];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    for dd = 1:length(dists)
        dists(dd) = distance(temp1(dd,:),target);
    end
    %if distance(subs(1,:),p2) > 300%disth  &&  distance(subs(1,:),p2) > disth 
    %if isempty(find(dists<disth)) ~= 1  %&&   min(newtime) < 600
        ct = ct +1;
    %if sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) < 300 && sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) > 300%dist
    %if newtime(1)<60*30  &&  distance(p1,p2) > 300  %newtime(1)> 60*5
    if 1==1
%     if length(find(temp1(:,2)>1000)) < length(find(temp1(:,2)<1000))  %approximately the gradient front~~
        

    %%%for angle
    dBs = zeros(1,size(vecs,1));
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    for dd = 2:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        %angs(dd) = angles(vecs(dd,:),target-subs(dd,:)); %ThetaInDegrees;
        %angs(dd) = angles(vecs(dd-1,:),vecs(dd,:)) / norm(vecs(dd-1));
        %%%complex plan method
        v1 = vecs(dd-1,:)/norm(vecs(dd-1,:));
        v2 = vecs(dd,:)/norm(vecs(dd,:));
        tempz1 = [v1(1)+v1(2)*(-1)^0.5];
        tempz2 = [v2(1)+v2(2)*(-1)^0.5];  %complex plane representation
        angs(dd) =  sign(angle(tempz1)-angle(tempz2))*(angles(vecs(dd-1,:),vecs(dd,:)) / norm(vecs(dd-1)));  %
        %%%asin method
        %angs(dd) = ...%sign( asin( dot( vecs(dd-1,:)/norm(vecs(dd-1,:)) , (vecs(dd-1,:)-vecs(dd,:))/norm(vecs(dd-1,:)-vecs(dd,:)) ) ) )*...
        %    (angles(vecs(dd-1,:),vecs(dd,:)) / norm(vecs(dd-1)));  %curving rate
        %if isreal(angs(dd))~=1
        %    break
        %end
        v1 = vecs(dd-1,:)/norm(vecs(dd-1,:));
        v2 =(subs(dd,:) - target)/norm(subs(dd,:) - target);
        tempz1 = [v1(1)+v1(2)*(-1)^0.5];
        tempz2 = [v2(1)+v2(2)*(-1)^0.5];  %complex plane representation
        dBs(dd) = sign(angle(tempz1)-angle(tempz2))*(angles(vecs(dd-1,:),(subs(dd,:) - target)));  %+/- angle to the targetarcsin(1/sqrt(2))
        
        %perpendicular concentration change
        perp_dir = [-vecs(dd-1,2), vecs(dd-1,1)];
        perp_dir = perp_dir/norm(perp_dir);
        dCp(dd) = Est_con(subs(dd-1,1)+perp_dir(1)*1., subs(dd-1,2)+perp_dir(2)*1, target(1), target(2), 50)...
            -Est_con(subs(dd-1,1)-perp_dir(1)*1, subs(dd-1,2)-perp_dir(2)*1, target(1), target(2), 50);
        %(temp1(dd-1,1),temp1(dd-1,2),target(1),target(2),50)
        %perp_dir = np.array([-dxy[1], dxy[0]])
        %perp_dir = perp_dir/np.linalg.norm(perp_dir)
        %perp_dC = gradient(C0, xx+perp_dir[0], yy+perp_dir[1]) - gradient(C0, xx-perp_dir[0], yy-perp_dir[1])
        
        %forward concentration change
        dCs(dd) = Est_con(temp1(dd-1,1),temp1(dd-1,2),target(1),target(2),50);
    end
    
    %%%for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(abs(real(angs)),p_thr,vecs);
    %%%remival
%     angs(timestamps) = [];
%     dBs(timestamps) = [];
%     dCp(timestamps) = [];
%     dCs(timestamps) = [];
    %angs(isreal(angs)~=1) = [];
    %dBs(isreal(angs)~=1) = [];
%     dC = [];
%     for tt = 1:length(timestamps)
%         dC = [dC Est_con(temp1(timestamps(tt),1),temp1(timestamps(tt),2),target(1),target(2),100)];
%     end
    allas = [allas angs];
    alldB = [alldB dBs];
    alldcp = [alldcp dCp];
    alldC = [alldC dCs];
    
    end
    
    
end

%%%
rem = find(abs(allas)>360);
allas(rem) = [];
alldB(rem) = [];
alldcp(rem) = [];
alldC(rem) = [];
plot(alldB,allas,'o')
figure; plot(alldcp,allas,'o')
%% plot adaptive binning
bins = 100;
dC_ = alldB;
dA_ = allas;

avang = zeros(1,bins);
stdang = zeros(1,bins);
[cts,bis] = hist(dC_,bins);
for bi = 2:length(bis)
    pos = find(bis(bi-1)<dC_ & bis(bi)>dC_);
    if isempty(pos)~=1
        avang(bi) = mean(dA_(pos));%length(dA_(pos)>threshold);%/(length(pos)+1);%
        stdang(bi) = std(dA_(pos));
    end
end
%plot(bis(avang~=NaN),avang(avang~=NaN),'-o')
%plot(bis(avang~=0),avang(avang~=0),'-o')
errorbar(bis,avang/0.3571,stdang)

% b = glmfit(dC_,dA_,'binomial','link','logit')
% plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)

%% Log-Likelihood method
%objective function to minimize: nLL_chemotaxis(THETA, dth, dcp, dc)
%a_,k_,A_,B_ = THETA
f = @(x)nLL_chemotaxis(x,allas,alldcp,alldC);
%[x,fval] = fminunc(f,rand(1,4));%[0,5,0.1,0.1]);
[x,fval] = fminunc(f,[0.5, 100, 0.5, 50]+rand(1,4));
x
fval

%% Generative model
figure;
%dth = N(alpha*dcp,K) + (A/(1+exp(B*dC)))*U[0,2*pi];
alpha = x(1);  kappa = (1/x(2))^0.5*(180/pi);  A = x(3);  B = x(4);
%alpha = 1.1;  kappa = 0.1*(180/pi);  A = 0.4;  B = 20;
test = [];
C0 = 50000000;
origin = [target(1)/2,target(2)];%[0,0];
target2 = target-[100,0];%origin+[500,0];%-[target(1)/2,target(2)];
for rep = 1:20
    
T = 1000;
dt = 1;
vm = 2.5;  %should be adjusted with the velocity statistics~~
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
    dC = Est_con(tracks(t-1,1),tracks(t-1,2),target2(1),target2(2),C0) - Est_con(tracks(t-2,1),tracks(t-2,2),target2(1),target2(2),C0);
    perp_dir = [-dxy(2) dxy(1)];
    perp_dir = perp_dir/norm(perp_dir);
    dCp = Est_con(tracks(t-1,1)+perp_dir(1),tracks(t-1,2)+perp_dir(2),target2(1),target2(2),C0) - Est_con(tracks(t-1,1)-perp_dir(1),tracks(t-1,2)-perp_dir(2),target2(1),target2(2),C0);
    
    wv = -alpha*dCp + kappa*randn;
    P_event = A/(1+exp(B*dC*dt));
    if rand < P_event
        beta = 1;
    else
        beta = 0;
    end
    rt = beta*(rand*360-180);
    dth = wv+rt;
    if dth>180; dth = dth-180; end;  if dth<-180; dth = dth+360; end  %within -180~180 degree range
    
    dths(t) = dth;
    dcs(t) = dC;
    dcps(t) = dCp;
    
    vv = vm+vs*randn;
    ths(t) = ths(t-1)+dth*dt;
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
plot(target2(1),target2(2),'ro')
plot(origin(1),origin(2),'ko')%(target2(1)/2,target2(2),'ko')

end


