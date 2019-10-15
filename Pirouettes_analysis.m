%Pirouettes_analysis
%101419
clear
clc
%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Frames','SmoothX','SmoothY'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% pre-processing
figure;
M = imread('Z:\Kevin\20190817\Data20190817_165031\Frame_000000.jpg');
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);

pix2mm = 1/31.5;  %pixel to mm (camera position before the flow chanber setup)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
nn = length(Tracks); %number of worms selected
mint = 60*1; %minimum time in seconds
minx = 100;  %minimum displacement (in terms of pixels)
target = [2517,975];  %position of target/sourse of odorant (approximated from images)
endingt = 60*15;  %only taking the first few minutes

hold on
cand = [];  %index of tracks as candidates for analysis
alldists = [];
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint  %time cutoff
        displace = mean((Tracks(i).Path(:,1)-mean(Tracks(i).Path(:,1))).^2 + (Tracks(i).Path(:,2)-mean(Tracks(i).Path(:,2))).^2); %pixel displacement
        alldists = [alldists displace*pix2mm^2];  %all dispacements in mm
        if displace > minx^2  %space cutoff
            pos = find(Tracks(i).Time>endingt);  %time window cutoff (the later time points are less correct...)
            if isempty(pos)~=1
                 plot(Tracks(i).Path(pos,1)*pix2mm,Tracks(i).Path(pos,2)*pix2mm,'k'); hold on;
                 cand = [cand i];
            end
        end
    end
end

%% turning rate
bin = 14;  %down sampling
filt = 14;  %smoothing the path
allTR = {};   %accumulating angle rate

for c = 1:length(cand)
    
    %pre-processing for each track
    id = cand(c);
    temp = [Tracks(id).SmoothX; Tracks(id).SmoothY]';  %Tracks(id).Path;
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
    
    %calculate angles
    dCs = zeros(1,size(vecs,1));  
    for dd = 2:length(angs)
        tempa = angles(vecs(dd-1,:),vecs(dd,:));
        %tempa = wrapToPi(tempa);
        angs(dd) = tempa;
        %dCs(dd) = Est_con(subs(dd-1,1),subs(dd-1,2),target(1),target(2),50);
    end
    
    allTR{c} = angs(2:end)/(fr*bin);  %all truning rates (angles (-180 to 180) per second)
    
end

%% cutoff with defined turns
turn_threshold = 120;   %<--visually identified threshold in the original Pirouette paper
allSP = {};  %all binary marker for turn or not
allRD = {};  %all run durations
for c = 1:length(cand)
    id = cand(c);
    pos_turn = find(abs(allTR{c})>=turn_threshold);
    tempSP = zeros(1,length(allTR{c}));
    tempSP(pos_turn) = 1;
    allSP{c} = tempSP;
    
    runs = diff(pos_turn)*(fr*bin);  %duration between turns in seconds
    if isempty(runs)==1
        allRD{c} = [];
    else
        allRD{c} = runs;
    end
end

%% plotting
histogram(cell2mat(allTR),100,'Normalization','probability')
xlabel('turning rate (deg/s)')
ylabel('pdf')

%% run durations
[cnt,dur] = hist(cell2mat(allRD),5000);
pos_z = find(cnt==0);
dur(pos_z) = [];
cnt(pos_z) = [];
% cutoff = 150;
pos_c = find(dur>=cutoff);
dur(pos_c) = [];
cnt(pos_c) = [];

%fit and find critical time
%test with doulbe exponent
Fexponent = fittype('a*exp(b*x)+c*exp(d*x)','dependent',{'y'},'independent',...
{'x'},'coefficients',{'a', 'b', 'c', 'd'});  %'e'
xVals = dur;
normf = max(cnt);
rVals = cnt/normf;
%all the fitting options required
minWindow = 1;
fitOptions = fitoptions(Fexponent);
fitOptions.Lower = [0,-.2,0];
fitOptions.Upper = [1000,0,10000];
fitOptions.StartPoint=[range(rVals(rVals~=0)),-.0001,range(rVals(rVals~=0)),-.0001];
fitOptions.Weights=zeros(size(rVals));
fitOptions.Weights(minWindow:end-minWindow)=1;
        
%do exponential fitting
[f,fout] = fit(xVals',rVals',Fexponent,fitOptions);
ts = max(f.b,f.d);
tl = min(f.b,f.d);
As = max(f.a,f.c);
Al = min(f.a,f.c);
tc = 1/(ts-tl)*log(As/Al)
plot(dur,cnt,'-o')
hold on
recon = f.a*exp(f.b*dur)+f.c*exp(f.d*dur);  %reconstruction from the exponential fitting
recon1 = f.a*exp(f.b*dur);
recon2 = f.c*exp(f.d*dur);
plot(dur,normf*recon,'r--'); hold on; plot(dur,normf*recon1,'g--'); hold on; plot(dur,normf*recon2,'b--')
xlabel('run duration (s)')
ylabel('counts')
%% bearing angle
allBe = {};  %all bearing angles
delta_B = 10;  %time steps measuring after the turn
for c = 1:length(cand)
    %find turns
    pos_turn = find(allSP{c}==1);
    
    if isempty(pos_turn)==0        
        id = cand(c);
        temp = [Tracks(id).SmoothX; Tracks(id).SmoothY]';  %Tracks(id).Path;
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
    
        %calculate angles
        temp_bear = [];
        for dd = 1:length(pos_turn)
            dd_ = pos_turn(dd);
            
            if dd_+delta_B<=length(temp2)  %conditioning on not out of the vector size
            tempa = angles(vecs(dd+1,:),[target-subs(dd,:)]);  %angle between next heading and the direction to target
            temp_bear = [temp_bear tempa];
            end
        end
        allBe{c} = temp_bear;
    else
        allBe{c} = [];
    end
    
end

%% plotting
histogram(cell2mat(allBe),100,'Normalization','probability')
xlabel('bearing angle (deg)')
ylabel('pdf')

%% conditioning on odor


