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

%visualize all paths and check criteria
cand = [];
figure;
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        if (Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2 > minx^2
            plot(Tracks(i).Path(:,1),Tracks(i).Path(:,2)); %axis([1500,2500,0,2000]); pause();
            hold on
            cand = [cand i];
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

%% approximate diffusion coefficient thrugh space
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
threshold = 30;
dC_ = diff(alldC);
dA_ = allas;
dA_(isnan(allas)) = 0;
dA_ = dA_(1:end-1);
%pos = find(abs(dC_)<0.00001);
%dC_(pos) = [];  dA_(pos) = [];

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
%% Logitic
CC = dCbin; %bis(avang~=NaN);
AA = avang/max(avang); %avang(avang~=NaN);

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
p_thr = 30;
ct = 0;
allas = [];
alldB = [];

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
    end
    
    %%%for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(abs(real(angs)),p_thr,vecs);
    %%%remival
    angs(timestamps) = [];
    dBs(timestamps) = [];
    %angs(isreal(angs)~=1) = [];
    %dBs(isreal(angs)~=1) = [];
%     dC = [];
%     for tt = 1:length(timestamps)
%         dC = [dC Est_con(temp1(timestamps(tt),1),temp1(timestamps(tt),2),target(1),target(2),100)];
%     end
    allas = [allas angs];
    alldB = [alldB dBs];
    
    end
    
    
end

plot(alldB,allas,'o')

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



