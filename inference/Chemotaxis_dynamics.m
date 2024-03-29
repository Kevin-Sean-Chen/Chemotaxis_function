%Chemotaxis_dynamics
%%%
% test with inference on turning rate and bearing angle in this script
%%%
%083019
clear
clc
%%
%addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
%addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')
addpath('/home/kschen/github/leifer-Behavior-Triggered-Averaging-Tracker-new')
addpath('/home/kschen/github/Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Frames','SmoothSpeed'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% % criteria %%%
figure();
% M = imread('Z:\Kevin\20190817\Data20190817_165031\Frame_000000.jpg');
M = imread('/tigress/LEIFER/Kevin/20190817/Data20190817_165031/Frame_000000.jpg');
% M = imread('Z:\Kevin\20191113_GWN_N2_app+\Data20191113_152718\Frame_000000.jpg');
% pix2mm = 1/16.5;
pix2mm = 1/31.5;  %pixel to mm (camera position before the flow chanber setup) %%%16.5 for new camera and 31.5 for old one
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
nn = length(Tracks); %number of worms selected
mint = 60*3;%60*1; %minimum time in seconds
minx = 100*1;  %minimum displacement (in terms of pixels)
disth = 300;  %radius of pixels from target
target = [2517,975];%[2000 750];%[950,1100];%  %position of target/sourse of odorant (approximated from images)%[250,1750];%
endingt = 60*30;  %only taking the first few minutes

%visualize all paths and check criteria
cand = [];
% figure;
hold on
alldists = [];
filt = 30;  %filtering time window
poly_degree = 3;  %polynomial fit for the tracks
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); 
        alldists = [alldists displace*pix2mm^2];
        if displace > minx^2  %space cutoff
            pos = find(Tracks(i).Time<endingt);
            if isempty(pos)~=1
                temp = Tracks(i).Path;
                temp1 = zeros(round(size(temp,1)/1),2);
                temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
                temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
                plot(temp1(:,1)*pix2mm, temp1(:,2)*pix2mm,'k'); %pause();%axis([1500,2500,0,2000]); 
%                 plot(Tracks(i).Path(pos,1)*pix2mm,Tracks(i).Path(pos,2)*pix2mm,'k');  %pause();%axis([1500,2500,0,2000]); 
%                 plot(Tracks(i).Path(:,1)*pix2mm,Tracks(i).Path(:,2)*pix2mm,'k');  %pause();%axis([1500,2500,0,2000]); 
%                 hold on
                cand = [cand i];
            end
        end
    end
end

set(gca,'fontsize',20)
xlabel('X (mm)')
ylabel('Y (mm)')
%% First-passage time measurement
crossing = 250;  %criteria in pixel space
timing = [];
figure;
for c = 1:length(cand)
    pos = find(Tracks(cand(c)).Path(:,1)>crossing);
    if size(pos) ~= 0
        timing = [timing Tracks(cand(c)).Time(pos(1))];  %record the time spent to the first crossing
    end
end
hist(timing)
%% Smoothed speed
speed = [];
for i = 1:nn
    speed = [speed Tracks(i).SmoothSpeed];  %speed distribution
end
hist(speed)
%% approximate diffusion coefficient through time
partt = 5;
Mt = 30*60;
bin_times = [0:Mt/(partt):Mt];   %time bins

Ds = [];
for tt = 1:length(bin_times)-1
    sub_tracks = FilterTracksByTime(Tracks,bin_times(tt)*14,bin_times(tt+1)*14);  %Mochi's code to select tracks within a time window
    xsqu = 0;
    for ii = 1:length(sub_tracks)
        xsqu = xsqu + mean(diff(sub_tracks(ii).Path(:,1)).^2);  %mean square displacement
    end
    Ds = [Ds xsqu/length(sub_tracks)];
end

hist(Ds)
%% approximate diffusion coefficient through space
partt = 6;
Mt = 30*60;
bin_space = [2592/2: (2592/2)/partt: 2592];  %spatial bins with frame as 1944 X 2592

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
hist(Ds)
%% %%%%% Peurette
%% Estimating concentration
figure;
bin = 10;  %bin for subsampling
filt = 30;  %filtering time window
poly_degree = 3;  %polynomial fit for the tracks
p_thr = 80;  %threshold of turning angle to define discrete events
l_window = 2;  %lag in time window to detect concentration
ct = 0;
allas = [];
alldC = [];

timelim = 60*10;

for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;
    temp1 = zeros(round(size(temp,1)/1),2);
%     temp1(:,1) = smooth(temp(1:round(length(temp)/1),1),filt);
%     temp1(:,2) = smooth(temp(1:round(length(temp)/1),2),filt);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [[0,0]; diff(subs)];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    for dd = 1:length(dists)
        dists(dd) = distance(subs(dd,:),target);
    end
    ct = ct +1;
    
    %%%condition on paths
    %if distance(subs(1,:),p2) > 300%disth  &&  distance(subs(1,:),p2) > disth 
    %if isempty(find(dists<disth)) ~= 1  %&&   min(newtime) < 600    
    %if sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) < 300 && sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) > 300%dist
    %if newtime(1)<60*30  &&  distance(p1,p2) > 300  %newtime(1)> 60*5
    %if length(find(temp1(:,2)>1000)) < length(find(temp1(:,2)<1000))  %approximately the gradient front~~
    if 1==1  

    %%%for angle
    dCs = zeros(1,size(vecs,1));  
    for dd = l_window+1:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        %angs(dd) = angles(vecs(dd,:),target-subs(dd,:)); %ThetaInDegrees;%
        angs(dd) = angles(vecs(dd-l_window,:),vecs(dd,:));
        dCs(dd) = Est_con(subs(dd-l_window,1),subs(dd-l_window,2),target(1),target(2),50);
    end
    angs = angs(l_window+1:end);
    dCs = dCs(l_window+1:end);
    
    %%%for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(real(angs),p_thr,vecs);
%     dC = [];
%     for tt = 1:length(timestamps)
%         dC = [dC Est_con(temp1(timestamps(tt),1),temp1(timestamps(tt),2),target(1),target(2),100)];
%     end
    allas = [allas angs(2:end)];
    alldC = [alldC dCs(2:end)];
    
    end
    %subplot(121);plot(subs(:,1),subs(:,2)); subplot(122);plot(dCs); pause();
    %subplot(121);scatter(subs(:,1),subs(:,2),[],1:size(subs,1)); xlim([1250,2500]);ylim([0,1800]); subplot(122);scatter(2:length(dCs),dCs(2:end),[],2:size(dCs,2)); pause();
    
end

% plot(diff(alldC),allas(1:end-1),'o')
plot((alldC),allas(1:end),'o')


%% test with design matrix
% win = 50;
% allC = zeros(length(alldC)-win,win+1);
% for ii = 1:10000%size(allC,1)
%     allC(ii,:) = alldC(ii:ii+win);%diff(alldC(ii:ii+win))/(alldC(ii)+0.000001);
% end
% imagesc(cov(allC))  %covariance of the concentration change in a sliding window
%% plot adaptive binning
bins = 50;  %number of adaptive bins to begin with
threshold = 60;  %threshold for a random sharp turn
dC_ = alldC;%diff(alldC);
dA_ = allas;
dA_ = dA_;%(1:end-1);
pos = find(abs(dC_)<1e-9);  %remove super small dC due to jittering
dC_(pos) = [];  dA_(pos) = [];

avang = [];
adapt_b = [];
[sortedC,sortedI] = sort(dC_);
rough_in_bins = floor(length(dC_)/bins);

for bi = 1:bins-1
    pos = sortedI((bi-1)*rough_in_bins+1:bi*rough_in_bins);
    avang = [avang sum(abs(dA_(pos))>threshold)/(length(pos)+1)];
    adapt_b = [adapt_b mean(sortedC(pos))];
end
pos = sortedI(bi*rough_in_bins:end);
avang = [avang sum(abs(dA_(pos))>threshold)/(length(pos)+1)];
adapt_b = [adapt_b mean(sortedC(pos))];
%[cts,bis] = hist(dC_,bins);
% for bi = 2:length(bis)
%     pos = find(bis(bi-1)<dC_ & bis(bi)>dC_);
%     if isempty(pos)~=1 && cts(bi-1)~=0
% %         avang(bi) = length(dA_(pos)>threshold)/(length(pos)+1);%%mean(dA_(pos));%
%         avang = [avang length(dA_(pos)>threshold)/(length(pos)+1)];
%         adapt_b = [adapt_b bis(bi)];
%     end
% end
plot(adapt_b,avang,'o')
% plot(bis(find(isnan(avang)==0)),avang(find(isnan(avang)==0)),'-o')
%plot(bis(avang~=0),avang(avang~=0),'-o')

% b = glmfit(dC_,dA_,'binomial','link','logit')
% plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)

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
%% custom Logistic
%test with doulbe exponent
Fexponent = fittype('a/(1+exp(b*x)) + c','dependent',{'y'},'independent',...
{'x'},'coefficients',{'a', 'b', 'c'});  %'e'
xVals = adapt_b-mean(adapt_b);
rVals = avang;
%all the fitting options required
minWindow = 1;
fitOptions = fitoptions(Fexponent);
fitOptions.Lower = [0,0.0001,0];
fitOptions.Upper = [1,100000,1];
fitOptions.StartPoint=[1,0.01,0.01];%[range(rVals(rVals~=0)),-.0001,range(rVals(rVals~=0)),-.0001];
fitOptions.Weights=zeros(size(rVals));
fitOptions.Weights(minWindow:end-minWindow)=1;
        
%do exponential fitting
[f,fout] = fit(xVals',rVals',Fexponent,fitOptions);
f
hold on
plot(adapt_b,f.a./(1+exp(f.b*(adapt_b-mean(adapt_b))))+f.c)
%% Logistic
CC = adapt_b*10000;%bis;%dCbin; %bis(avang~=0);%
AA = avang%/max(avang); %avang(avang~=0);%
% pos = CC<-0.01;
% CC(pos) = [];
% AA(pos) = [];

[b,dev,stats] = glmfit(CC, AA', 'binomial', 'link', 'logit')
xx = linspace(min(CC), max(CC), 30);
yfit = glmval(b,xx,'logit');
plot(xx,yfit,'-')
hold on
plot(CC,AA,'ro')

%% or ignore the binning part and fit directly to binary data
[b,dev,stats] = glmfit(diff(alldC), temp_turn(1:end-1)', 'binomial', 'link', 'logit')
xx = linspace(min(CC), max(CC), 30);
yfit = glmval(b,xx,'logit');
plot(xx,yfit,'-')
hold on
plot(diff(alldC),temp_turn(1:end-1),'ro')

%% %%%%% Weathervaning
%% %%%%%
figure;
bin = 5;
filt = 30;
p_thr = 60;
l_window = 20;  %lag time
ct = 0;
allas = [];
alldB = [];
alldcp = [];
alldC = [];
alldis = [];

timelim = 60*10;

for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = diff(subs);%[[0,0]; diff(subs)];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    for dd = 1:length(dists)
        dists(dd) = distance(subs(dd,:),target);
    end
    %if distance(subs(1,:),p2) > 300%disth  &&  distance(subs(1,:),p2) > disth 
    %if isempty(find(dists<disth)) ~= 1  %&&   min(newtime) < 600
        ct = ct +1;
    %if sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) < 300 && sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) > 300%dist
    %if newtime(1)<60*30  &&  distance(p1,p2) > 300  %newtime(1)> 60*5
    if 1==1
%     if length(find(temp1(:,2)>1000)) < length(find(temp1(:,2)<1000))  %approximately the gradient front~~
        

    %%%for time step calculations
    dBs = zeros(1,size(vecs,1));
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    for dd = l_window+1:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        %angs(dd) = angles(vecs(dd,:),target-subs(dd,:)); %ThetaInDegrees;
        %angs(dd) = angles(vecs(dd-1,:),vecs(dd,:)) / norm(vecs(dd-1));
        %%%complex plan method
%         v1 = vecs(dd-1,:)/norm(vecs(dd-1,:));
%         v2 = vecs(dd,:)/norm(vecs(dd,:));
%         tempz1 = [v1(1)+v1(2)*(-1)^0.5];
%         tempz2 = [v2(1)+v2(2)*(-1)^0.5];  %complex plane representation
%         angs(dd) =  sign(angle(tempz1)-angle(tempz2))*(angles(vecs(dd-1,:),vecs(dd,:)) / norm(vecs(dd-1)));  %
        %%%angle function
%         angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
%         angs(dd) = angles(sum(vecs(dd-l_window:dd,:))/norm(sum(vecs(dd-l_window:dd,:))),vecs(dd,:)/norm(vecs(dd,:)));
        %%%asin method
        %angs(dd) = ...%sign( asin( dot( vecs(dd-1,:)/norm(vecs(dd-1,:)) , (vecs(dd-1,:)-vecs(dd,:))/norm(vecs(dd-1,:)-vecs(dd,:)) ) ) )*...
        %    (angles(vecs(dd-1,:),vecs(dd,:)) / norm(vecs(dd-1)));  %curving rate
        %if isreal(angs(dd))~=1
        %    break
        %end
%         v1 = vecs(dd-1,:)/norm(vecs(dd-1,:));
%         v2 =(subs(dd,:) - target)/norm(subs(dd,:) - target);
%         tempz1 = [v1(1)+v1(2)*(-1)^0.5];
%         tempz2 = [v2(1)+v2(2)*(-1)^0.5];  %complex plane representation
%         dBs(dd) = sign(angle(tempz1)-angle(tempz2))*(angles(vecs(dd-1,:),(subs(dd,:) - target)));  %+/- angle to the targetarcsin(1/sqrt(2))
%         dBs(dd) = angles(vecs(dd,:)/norm(vecs(dd,:)),(target-subs(dd-1,:))/norm(subs(dd-1,:)-target));
        dBs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),(subs(dd-l_window,:)-target)/norm(subs(dd-l_window,:)-target));
%         dBs(dd) = angles(sum(vecs(dd-l_window:dd,:))/norm(sum(vecs(dd-l_window:dd,:))),(subs(dd-l_window,:)-target)/norm(subs(dd-l_window,:)-target));
        
        %perpendicular concentration change
        perp_dir = [-vecs(dd-l_window,2), vecs(dd-l_window,1)];
        perp_dir = perp_dir/norm(perp_dir);
        dCp(dd) = Est_con(subs(dd-l_window,1)+perp_dir(1)*1., subs(dd-l_window,2)+perp_dir(2)*1, target(1), target(2), 50)...
                  -Est_con(subs(dd-l_window,1)-perp_dir(1)*1, subs(dd-l_window,2)-perp_dir(2)*1, target(1), target(2), 50);
        %(temp1(dd-1,1),temp1(dd-1,2),target(1),target(2),50)
        %perp_dir = np.array([-dxy[1], dxy[0]])
        %perp_dir = perp_dir/np.linalg.norm(perp_dir)
        %perp_dC = gradient(C0, xx+perp_dir[0], yy+perp_dir[1]) - gradient(C0, xx-perp_dir[0], yy-perp_dir[1])
        
        %forward concentration change
        dCs(dd) = Est_con(subs(dd-l_window,1),subs(dd-l_window,2),target(1),target(2),50);
        
        %%%check displacement
        dds(dd) = norm(vecs(dd,:));
    end
    %remove zeros
    dBs = dBs(l_window+1:end);
    dCp = dCp(l_window+1:end);
    dCs = dCs(l_window+1:end);
    dds = dds(l_window+1:end);
    angs = angs(l_window+1:end);
    
    %%%for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(abs(real(angs)),p_thr,vecs);
    %%%removal
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
    alldis = [alldis dds];  %displacement (effective velocity)
    
    end
    
    
end

%%%
% rem = find(abs(allas)>360);
% allas(rem) = [];
% alldB(rem) = [];
% alldcp(rem) = [];
% alldC(rem) = [];
plot(alldB,allas,'o')
figure; plot(alldcp,allas,'o')
% histogram(alldis)

%% removing small jitterings (test)
vec_threshold = 0.1;%median(alldis);
vec_marks = find(alldis > vec_threshold);

%% removing large anlge changes (test)
ang_threshold = 20;%median(alldis);
ang_marks = find(abs(allas) < ang_threshold);

va_marks = intersect(vec_marks, ang_marks);
%% plot adaptive binning
bins = 60;
dC_ = alldB(va_marks);
dA_ = allas(va_marks)./(alldis(va_marks)*pix2mm);  %curving rate (deg/mm)

avang = zeros(1,bins);
stdang = zeros(1,bins);
h = histogram(dC_,bins);
cts = h.Values;
bis = h.BinEdges(1:end-1);
for bi = 2:length(bis)
    pos = intersect(find(bis(bi-1)<dC_),find((bis(bi)>=dC_)));
%     pos = find((bis(bi-1)<dC_) & (bis(bi)>dC_));
    if isempty(pos)~=1
        avang(bi) = mean(dA_(pos));%length(dA_(pos)>threshold);%/(length(pos)+1);%
        stdang(bi) = std(dA_(pos));
    end
end
%plot(bis(avang~=NaN),avang(avang~=NaN),'-o')
%plot(bis(avang~=0),avang(avang~=0),'-o')
errorbar(bis,avang,stdang)

% b = glmfit(dC_,dA_,'binomial','link','logit')
% plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)

%% Log-Likelihood method
%objective function to minimize: nLL_chemotaxis(THETA, dth, dcp, dc)
%a_,k_,A_,B_ = THETA
f = @(x)nLL_chemotaxis(x,allas(1:end-1),alldcp(1:end-1),diff(alldC));
[x,fval] = fminunc(f,rand(1,4));%[0,5,0.1,0.1]);  %random initiation
% [x,fval] = fminunc(f,[0.5, 10, 0.5, 50]+rand(1,4));  %a closer to a reasonable value
x
fval

%% Generative model
figure;
%dth = N(alpha*dcp,K) + (A/(1+exp(B*dC)))*U[0,2*pi];
alpha = x(1);  kappa = (1/x(2))^0.5*(180/pi);  A = x(3);  B = x(4);
%alpha = 1.1;  kappa = 0.1*(180/pi);  A = 0.4;  B = 20;
test = [];
C0 = 100000;
origin = [target(1)/2,target(2)];%[0,0];
target2 = target-[100,0];%origin+[500,0];%-[target(1)/2,target(2)];

alldths = [];
alldCs = [];
alldCps = [];
allths = [];
for rep = 1:30
    
T = 2000;
dt = 1;
vm = 5.0;  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
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
    if ths(t)>180; ths(t) = ths(t)-180; end;  if ths(t)<-180; ths(t) = ths(t)+360; end  %within -180~180 degree range
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

alldths = [alldths dths];
alldCs = [alldCs dcs];
alldCps = [alldCps dCp];
allths = [allths ths];

end


