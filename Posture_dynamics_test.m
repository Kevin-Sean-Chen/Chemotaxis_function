%Posture_dynamics_test
%% load eigen-weights and centerlines
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
fields_to_load = {'Path','Time','SmoothSpeed','Behaviors','LEDPower','Centerlines','AngSpeed','ProjectedEigenValues'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);%(1:3)

allt = 0; for ii = 1:length(Tracks); temp = Tracks(ii).Time; allt = allt+ temp(end)-temp(1); end  %animal-hours

%% load eigen worms
load('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Centerline\EigenVectors.mat')

%% find a1 and a2
alla12s = [];
for ii = 1:length(Tracks)
    alla12s = [alla12s Tracks(ii).ProjectedEigenValues(1:3,:)];  %second and third modes are sine and cosine
end
H = histogram2(alla12s(2,:),alla12s(3,:));
imagesc(H.Values)
%% compute phi
filt = 1;
phi = zeros(1,size(alla12s,2));
a1 = smooth(1*alla12s(2,:)/var(alla12s(2,:)),filt);
a2 = smooth(1*alla12s(3,:)/var(alla12s(3,:)),filt);
for tt = 1:length(phi)
    phi(tt) = atan(a1(tt)/a2(tt));  %phase angle in the a1-a2 phase plane
end

ww = diff(phi);
phi = phi(1:end-1);
H2 = histogram2(ww/(2),phi,[100,100]);
imagesc(H2.Values)

%% speed versus agnular
allvs = [];
for ii = 1:length(Tracks)
    allvs = [allvs Tracks(ii).SmoothSpeed];
end
allvs = smooth(allvs(1:end-1),1)';  %smoothing out jittering speed
pos = find(allvs>0.4);  %remove fast reversal or artifact
allvs(pos) = [];
ww = diff(phi);
ww(pos(1:end-1)) = [];
H = histogram2(abs(allvs),abs(ww));
imagesc(H.Values,'XData',[0 max(allvs)],'YData',[0 max(ww)]);
set(gca,'YDir','normal');

%% stimulus correlation
allSs = [];
for ii = 1:length(Tracks)
    allSs = [allSs Tracks(ii).LEDPower];
end
plot(xcorr(allSs,alla12s(1,:)))

%% stimulus modulated covariance
win = 100;
len = 10000;
xx = hankel(1:len-win, len-win:len);  %design matrix in a time window
stim = allSs(1:len);%-mean(allSs(1:len));
sinp = sin(phi(1,1:len));%-mean(sin(phi(1,1:len)));
a3_ = alla12s(1,1:len);%-mean(alla12s(1,1:len));
X = stim(xx) - repmat(mean(stim(xx),2),1,size(xx,2));
ss = sinp(xx) - repmat(mean(sinp(xx),2),1,size(xx,2));
aa = a3_(xx) - repmat(mean(a3_(xx),2),1,size(xx,2));

imagesc(ss'*aa)  %%%triggered covariance of a3 turns
