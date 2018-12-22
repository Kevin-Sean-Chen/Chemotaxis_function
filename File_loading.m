clear; clc;
%% File_loading
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
fields_to_load = {'Path','Time','Speed','Behaviors','LEDPower'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);%(1:3)
%filtered_tracks = FilterTracksByTime(T

allt = 0; for ii = 1:length(Tracks); temp = Tracks(ii).Time; allt = allt+ temp(end)-temp(1); end  %animal-hours
%% Chemotaxis analysis
load('Path'); path = values;

%%%angale analysis
bins = 1:2:10;
for bb = 1:length(bins)   
bin = bins(bb);
ns = 60;
figure;

for ii = ns
    %trace = getfield(Tracks,{ii+5},'Path');
    traces = path{ii};
    subs = traces(1:bin:end,:);
    vecs = diff(subs);

alpha = zeros(size(vecs,1)-1);
for vv = 1:size(vecs,1)-1
    u = vecs(vv,:);
    v = vecs(vv+1,:);
    CosTheta = dot(u,v)/(norm(u)*norm(v));
    ThetaInDegrees = acosd(CosTheta);
    alpha(vv) = ThetaInDegrees;
end
%subplot(ns,1,ii)
plot(alpha)
end

end

%%%visualize paths
figure;
cc = hsv(size(traces,1));
for ii = 1:size(traces,1)
    plot(traces(ii,1),traces(ii,2),'*','color',cc(ii,:)); 
    hold on
end

%% angle analysis  
bin = 4;
ns = 10;
figure;

for ii = 1:ns
    %trace = getfield(Tracks,{ii+5},'Path');
    traces = path{ii};
    subs = traces(1:bin:end,:);
    vecs = diff(subs);

alpha = zeros(size(vecs,1)-1);
for vv = 1:size(vecs,1)-1
    u = vecs(vv,:);
    v = vecs(vv+1,:);
    CosTheta = dot(u,v)/(norm(u)*norm(v));
    ThetaInDegrees = acosd(CosTheta);
    alpha(vv) = ThetaInDegrees;
end

%subplot(ns,1,ii)
plot(alpha)
hold on
ii

end


%% 

%load data
load('LEDPower.mat')
allS = values;
load('Behaviors.mat')
allB = values;

%STA
num = 3;
beh = 3;
win = 50;
acs = 10;
pos = find(allB{num}(beh,:) == 1);
trigs = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(allB{num}(beh,:))
        trigs(ii,:) = allS{num}(pos(ii)-win:pos(ii)+acs);
    end
end

%% all STA
alltrigs = [];
for ww = 1:435
    
num = ww;
beh = 4;
win = 140;
acs = 140;
pos = find(allB{num}(beh,:) == 1);
%pos = find(sum(allB{num}) >= 1);
trigs = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(allB{num}(beh,:))
        trigs(ii,:) = allS{num}(pos(ii)-win:pos(ii)+acs);
    end
end
alltrigs = [alltrigs; trigs];
end
plot(mean(alltrigs))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partial STA calculation
figure
bb = 8;
Mt = 30;  %30 min recordings
beh = bb; % behavioral state
win = 140;
acs = 140;
%%%time windows
partt = 6; %6
bin_times = [0:Mt/(partt):Mt]*60;
trknum = zeros(1,partt);
trknum2 = zeros(1,partt);
trknum3 = cell(1,partt);

for bt = 1:length(bin_times)-1
alltrigs = [];

for ww = 1:size(Tracks,2)
    
num = ww;

temp2 = Tracks(num).Time;
overlap = find(temp2>bin_times(bt) & temp2<bin_times(bt+1));
%%%only analyze overlaping time points
if isempty(overlap)~=1
%     trknum(bt) = trknum(bt)+1;
            
temp = Tracks(ww).Behaviors;
stim = Tracks(ww).LEDPower;

temp = temp(:,overlap);
stim = stim(overlap);

pos = find(temp(beh,:) == 1);

trknum(bt) = trknum(bt)+length(pos);
trknum2(bt) = trknum2(bt)+(temp2(end)-temp2(1));
trknum3{bt} = [trknum3{bt} length(pos)./(temp2(end)-temp2(1))];

trigs = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(temp(beh,:))
        trigs(ii,:) =stim(pos(ii)-win:pos(ii)+acs);
    end
end
alltrigs = [alltrigs; trigs];

end

end
% figure; 
Ker = mean(alltrigs);
%Ker = Ker-mean(Ker(1:100));
%Ker = smooth(Ker,20)';
%Ker = Ker/norm(Ker); 
plot([-acs:win]*(1/14),Ker,'Linewidth',2)
set(gca,'Fontsize',20)
hold on

end

%% partial NL
Ker = mean(alltrigs);
Ker = Ker-mean(Ker(1:100));
Ker = smooth(Ker,20)';
Ker = Ker/norm(Ker);  %%%unit vector
LN = [];
dur = 0;
bb = 8;
for ww = 1:size(Tracks,2)
    
    temp2 = Tracks(ww).Time;
    overlap = find(temp2>bin_times(bt) & temp2<bin_times(bt+1));
%%%only analyze overlaping time points
if isempty(overlap)~=1
    
beh = bb; % behavioral state
temp = Tracks(ww).Behaviors;
stim = Tracks(ww).LEDPower;

temp = temp(:,overlap);
stim = stim(overlap);
%%%
xx = conv(Ker,stim-33);  %try to zero mean both
xx = xx(length(Ker):end);
btt = temp(beh,:);

LN = [LN [xx; btt]];
end

end

%%  plot NL
figure
%out = find(LN(1,:)<-2000);
LN2 = LN;
%LN2(:,out) = [];
[nnn,xout] = hist(LN2(1,:),10);  %filtered signal
%cts = cell(1,length(xout)); %
cts = zeros(1,length(xout));
for in = 1:length(LN2(1,:))
    [min_val, index] = min(abs(xout-LN2(1,in)));
    cts(index) = cts(index)+LN2(2,in);  %real transitions
    %cts{index} = [cts{index} LN2(2,in)];
end
ct2 = nnn;
tr = cts./ct2*14*60;%(allt/60);
plot(xout,tr,'-o')

%% 
% figure
load('ADV2.mat')
hold on
EE = ( ((cts-1)./(ct2.^2)) + (cts.^2.*(ct2-1)./(ct2.^4)) ).^0.5 *14*60;
errorbar(xout,tr,EE,'-o','LineWidth',2)

%lump
mu_v = mean(tr(5:7));
sig_v = ( ((sum(cts(5:7))-1)./(sum(ct2(5:7)).^2)) + (sum(cts(5:7)).^2.*(sum(ct2(5:7))-1)./(sum(ct2(5:7)).^4)) ).^0.5 *14*60;

% bar(1:3,[mu_a,mu_v,mu_n])
% hold on
% errorbar([1,2,3],[mu_a,mu_v,mu_n],[sig_a,sig_v,sig_n],'.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all STA for pooled data
for bb = 8%1:size(Tracks(1).Behaviors,1)
alltrigs = [];

for ww = 1:size(Tracks,2)
num = ww;
beh = bb; % behavioral state
win = 140;
acs = 140;
temp = Tracks(ww).Behaviors;
stim = Tracks(ww).LEDPower;
pos = find(temp(beh,:) == 1);
%%%if isempty(pos)~=1; pos = randi(length(temp(beh,:)),length(pos)); end
%pos = find(sum(allB{num}) >= 1);
trigs = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(temp(beh,:))
        trigs(ii,:) =stim(pos(ii)-win:pos(ii)+acs);
    end
end
alltrigs = [alltrigs; trigs];
end
%figure; plot([-acs:win]*(1/14),mean(alltrigs))
set(gca, 'XTickLabel', [],'XTick',[])
subplot(size(Tracks(1).Behaviors,1),1,bb); plot([-acs:win]*(1/14),mean(alltrigs))
set(gca,'Fontsize',20)

end

%% NL part!
Ker = mean(alltrigs);
Ker = Ker-mean(Ker(1:100));
Ker = smooth(Ker,20)';
LN = [];
dur = 0;
bb = 8;
for ww = 1:size(Tracks,2)
beh = bb; % behavioral state
temp = Tracks(ww).Behaviors;
stim = Tracks(ww).LEDPower;
tts = Tracks(ww).Time;
%%%
xx = conv(Ker,stim-33);  %try to zero mean both
xx = xx(length(Ker):end);
bt = temp(beh,:);
% pos = find(bt == 1);
% bts = zeros(length(pos),win+acs+1);
% xxs = zeros(length(pos),win+acs+1);
% for ii = 1:length(pos)
%     if pos(ii)-win>0 && pos(ii)+acs<length(bt)
%         bts(ii,:) = bt(pos(ii)-win:pos(ii)+acs);
%         xxs(ii,:) = xx(pos(ii)-win:pos(ii)+acs);
%     end
% end
% 
% LN = [LN [xxs(:)'; bts(:)']];
LN = [LN [xx; bt]];
% dur = dur + (tts(end)-tts(1));
end

%%  plot NL... removing large values
figure
%out = find(LN(1,:)<-2000);
LN2 = LN;
%LN2(:,out) = [];
[nnn,xout] = hist(LN2(1,:),10);  %filtered signal
%cts = cell(1,length(xout)); %
cts = zeros(1,length(xout));
for in = 1:length(LN2(1,:))
    [min_val, index] = min(abs(xout-LN2(1,in)));
    cts(index) = cts(index)+LN2(2,in);  %real transitions
    %cts{index} = [cts{index} LN2(2,in)];
end
ct2 = nnn;
tr = cts./ct2*14*60;%(allt/60);
plot(xout,tr,'-o')

%% 
save('NAI3.mat','allt','alltrigs','folder_names','num','bb','xout','tr','cts','ct2')
%% 
EE = ( ((cts-1)./(ct2.^2)) + (cts.^2.*(ct2-1)./(ct2.^4)) ).^0.5 *14*60;
errorbar(xout,tr,EE,'-o','LineWidth',2)
% for ii = 1:length(xout)
%     cts(ii) = cts(ii)/ww/(allt/60)/ct2(ii);
%     errorbar(xout(ii),mean(cts(ii)),std(cts(ii)))
%     plot(xout(ii),mean(cts(ii)),'o')
%     hold on
% end

%% Pulse induced behavior
win = 100;
acs = 100;
thr = 0.5;  %all pulse should be larger than 0.5

for bb = 8%1:size(Tracks(1).Behaviors,1)

alltrigs = [];
allpI = [];
for ww = 1:size(Tracks,2)
beh = bb; % behavioral state
temp = Tracks(ww).Behaviors;
stim = Tracks(ww).LEDPower;
pos = find(diff(stim) > thr);  %find onset of pulses
pI = stim(pos+1);
trigs = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(temp(beh,:))
        trigs(ii,:) = temp(beh,pos(ii)-win:pos(ii)+acs);
    end
end
alltrigs = [alltrigs; trigs];
allpI = [allpI, pI];
end
figure; plot(sum(alltrigs))

end

binedg = [0,2,6,15,50,80];
res = zeros(1,length(binedg)-1);
for ii = 1:length(binedg)-1
    pos = find(allpI>binedg(ii) & allpI<binedg(ii+1));
    partI = mean(alltrigs(pos,:));
    res(ii) = max(partI);
    figure; plot(sum(alltrigs(pos,:)))
    hold on
end
%plot(binedg(2:end),res,'-o')


%% Kernel statistics
ALL = cell(1,3);
load('NAI_new.mat'); ALL{1} = alltrigs;
load('APP.mat'); ALL{2} = alltrigs;
load('ADV.mat'); ALL{3} = alltrigs;
%%
boot = 10;
SAMP = cell(1,3);
Ks = cell(1,3);

for aa = 1:3
temp = ALL{aa};
Mm = zeros(1,boot);

tt = size(temp,2);
nsamp = size(temp,1);
subsamp = floor(nsamp/boot)*boot;
temp = temp(1:subsamp,:);
temp = temp(randperm(subsamp),:);
temp2 = reshape(temp,boot,subsamp/boot,tt);
allK = [];
for bt = 1:boot
    tempK = squeeze(mean(temp2(bt,:,:),2));
    allK = [allK; tempK'];
Mm(bt) = max(tempK)-min(tempK);
end

SAMP{aa} = Mm;
Ks{aa} = allK;

end

xdata = [ones(boot,1),2*ones(boot,1),3*ones(boot,1)];
ydata = [SAMP{1}' SAMP{2}' SAMP{3}'];
scatter(xdata(:), ydata(:),100, 'r.', 'jitter','on', 'jitterAmount', 0.05);
hold on
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(mean(ydata, 1), 2, 1), 'k-')
hold off

%%
waves = cat(1, Ks{1}, Ks{2});
waves = cat(1,waves,Ks{3});
waves = waves-mean(waves);
C = cov(waves);
[U,S,V] = svd(C);

PC1 = waves*U(:,1);
PC2 = waves*U(:,2);


