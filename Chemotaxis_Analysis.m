%Chemotaxis_Analysis
%02272018_test
clear
clc
%%
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%load data from one folder
% Tracks = [];
% load('Path.mat')
% paths = values;
% load('Time.mat')
% times = values;
% for ii = 1:length(paths); Tracks(ii).Path = paths{ii}; end
% for ii = 1:length(pat0hs); Tracks(ii).Time = times{ii}; end

%or batch analysis
fields_to_load = {'Path','Time','Frames'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%%% criteria %%%
nn = length(Tracks); %number of worms selected
mint = 60;%60*1; %minimum time in seconds
disth = 500;  %radius of pixels from target
target = [2517,975];%[950,1100];%  %position of target/sourse of odorant (approximated from images)

%visualize all paths and check criteria
cand = [];
figure;
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        plot(Tracks(i).Path(:,1),Tracks(i).Path(:,2));% pause()
        hold on
        cand = [cand i];
    end
end

%% plot by time slots (with Mochi's code)
tt = zeros(1,length(cand));
for c = 1:length(cand)
    id = cand(c);
    temp2 = Tracks(id).Time;
    tt(c) = temp2(end);
end
Mt = max(tt);
%%%time windows
partt = 5;
bin_times = [0:Mt/(partt):Mt];

for tt = 1:length(bin_times)-1
    sub_tracks = FilterTracksByTime(Tracks,bin_times(tt)*14,bin_times(tt+1)*14);  %Mochi's code to select tracks within a time window
    figure;
    for i = 1:length(sub_tracks)
        plot(sub_tracks(i).Path(:,1),sub_tracks(i).Path(:,2))
        hold on
    end
end

%% plot for distance and angle
%test 
figure;
filt = 10;
bin = 3;
fract = [];
for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
    temp1 = zeros(size(temp));
    temp1(:,1) = smooth(temp(:,1),filt);
    temp1(:,2) = smooth(temp(:,2),filt);
    %temp1 = reshape(temp1,size(Tracks(id).Path,1),2);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [diff(subs); [0,0]];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    p1 = subs(end,:); p2 = target;
%     for dd = 1:length(dists)
%         dists(dd) = distance(temp1(dd,:),target);
%     end
    %if distance(p1,p2) < disth  &&  distance(subs(1,:),p2) > disth
%     if isempty(find(dists<disth)) ~= 1
    %if sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) < 300 && sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) > 300%dist
    %if 1==1
        
    %%%for distance
    for dd = 1:length(dists)
        dists(dd) = distance(subs(dd,:),target);
    end
    infrac = zeros(1,size(subs,1));
    infrac(dists<disth) = 1;
    fract = [fract sum(infrac)/length(infrac)];

    %%%for angle
    for dd = 2:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        angs(dd) = angles(vecs(dd,:),target-subs(dd,:));%ThetaInDegrees;%angles(vecs(dd-1,:),vecs(dd,:));  %
    end
    
%     subplot(4,1,1); plot(subs(:,1),subs(:,2)); hold on; plot(target(1),target(2),'r*'); hold off
%     subplot(4,1,2); plot(newtime,dists)
%     ylabel('distance')
%     title(num2str(id))
%     subplot(4,1,3); 
%     plot(newtime,infrac)
%     subplot(4,1,4); plot(newtime,angs)
%     ylabel('angle')
%     xlabel('time')
    %hold on
 
    angs(isnan(angs)==1) = 0;
    subplot(2,1,1);
    autocorr(temp1(:,1),size(temp1,1)-1);%( angs,length(angs)-1 ); %
    subplot(2,1,2);
    plot(angs)
    pause()
    
    %end
    
end

%% turning analysis
figure;
bin = 1;
filt = 10;
allas = [];
allps = [];
allrs = [];
allcs = [];
p_thr = 60;
ct = 0;

timelim = 60*10;

for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1),filt);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2),filt);
    %temp1 = reshape(temp1,size(Tracks(id).Path,1),2);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [[0,0]; diff(subs)];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    p1 = subs(end,:); p2 = target;
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
        
%     %%%for distance
%     for dd = 1:length(dists)
%         dists(dd) = distance(subs(dd,:),target);
%     end
%     infrac = zeros(1,size(subs,1));
%     infrac(dists<disth) = 1;
%     fract = [fract sum(infrac)/length(infrac)];

    %%%for angle
    for dd = 2:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        %angs(dd) = angles(vecs(dd,:),target-subs(dd,:)); %ThetaInDegrees;%
        angs(dd) = angles(vecs(dd-1,:),vecs(dd,:));
%         if isreal(angs(dd)) ~=1
%             c
%             dd
%             break;
%         end
    end
    
    %%%for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(real(angs),p_thr,vecs);
%     plot(timestamps); pause()
    
    allps = [allps length(timestamps)/(newtime(end)-newtime(1))];  %turn rate!
    allas = [allas angs];
    allrs = [allrs runs];
    
    [acf,lags,bounds] = autocorr( temp1(:,1),size(temp1,1)-1 );
    acf(acf<0) = 0;
    [val,pos] = min(acf-0.5);
    allcs = [allcs pos*bin/14];
    
    end
    
    
end

%hist(real(allas),100)
%hist(allas,100)
%hist(allrs,100)

% %functions in the end~
% function dd = distance(p1,p2)%p2=target
%     dd = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
% end
% function aa = angles(v,u)%u=target
%     CosTheta = dot(u,v)/(norm(u)*norm(v));
%     ThetaInDegrees = acosd(CosTheta);
%     aa = ThetaInDegrees;
% end
% function nt = normtime(time)
%     nt = time-time(1);
% end


%% trial-based turn analysis!!
figure;
bin = 1;  %use for down-sampling
filt = 10;  %filter trackes
allas = [];  %angles
allps = [];  %peurettes (sharp turns)
allrs = [];  %runs
allds = [];  %test variables
p_thr = 90;  %angle threshold
ct = 0;  %count incident concidered


for c = 1:length(cand)
    
    id = cand(c);%cand(randi([1,length(cand)]));%
    temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
    temp1 = zeros(round(size(temp,1)/1),2);
    temp1(:,1) = smooth(temp(1:round(length(temp)/1),1),filt);
    temp1(:,2) = smooth(temp(1:round(length(temp)/1),2),filt);
    %temp1 = reshape(temp1,size(Tracks(id).Path,1),2);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [[0,0]; diff(subs)];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));
    angT = zeros(1,size(vecs,1));
    
    %%% select criteria
    p1 = subs(end,:); p2 = target;
    for dd = 1:length(dists)
        dists(dd) = distance(temp1(dd,:),target);
    end
    %if distance(p1,p2) < disth  &&  distance(subs(1,:),p2) > disth 
    %if isempty(find(dists<disth)) ~= 1  %&&   min(newtime) < 600
    %if sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) < 300 && sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2) > 300%dist
    %if newtime(end)<1500  &&  newtime(1)<1500
    if 1==1
        ct = ct +1;
        
    %%% for distance
%     for dd = 1:length(dists)
%         dists(dd) = distance(subs(dd,:),target);
%     end
%     infrac = zeros(1,size(subs,1));
%     infrac(dists<disth) = 1;
%     fract = [fract sum(infrac)/length(infrac)];

    %%% for angle
    for dd = 2:length(angs)
        %CosTheta = dot(vecs(dd,:),(target-subs(dd,:)))/(norm(vecs(dd,:))*norm((target-subs(dd,:))));
        %ThetaInDegrees = acosd(CosTheta);
        angT(dd) = angles(vecs(dd-1,:),target-subs(dd-1,:)); %ThetaInDegrees;%
        angs(dd) = angles(vecs(dd-1,:),vecs(dd,:));
    end
    
    %%% for "Pirouttes" frequnecy
    [timestamps,runs] = def_turns(angs,p_thr,subs);
    %plot(timestamps); pause()
    
    prate = length(timestamps)/(newtime(end)-newtime(1));
    allps = [allps prate];  %turn rate!
    angs(2) = 0;
    R = corrcoef(angs,angT);
    dist_turn = [];
    if isempty(timestamps)~=1  %&& length(timestamps)<300
        for d = 1:length(timestamps)
            dist_turn = [dist_turn distance(temp1(timestamps+1,:),target)];
        end
        time_turn = newtime(timestamps);
    [aa,bb] = hist(angs(timestamps),[0:2:180]);
    %aa = aa/sum(aa);
    %plot(bb,aa);pause()
    allds = [allds  angT(timestamps)];%aa];%newtime];%time_turn];%R(2)];%(length(find(angT>90)))];  %%%turn time, angle, correlation...
    
    allas = [allas angs];  %%%all angles
    allrs = [allrs runs];  %%% all run lengths
    
    %%%visualization
%     subplot(1,2,1); hist(angT(timestamps),50); 
%     subplot(1,2,2); plot(temp1(:,1),temp1(:,2)); 
%     hold on; plot(temp1(timestamps,1),temp1(timestamps,2),'r*'); hold off;
%     pause()
    end
    
    end
    
    %%% visualize turns
%     plot(temp1(:,1),temp1(:,2),'b')
%     hold on
%     plot(temp1(timestamps+1,1),temp1(timestamps+1,2),'r*')

    
end

% hist(real(allds),100)
% plot(bb,sum(allds))

[yy,xx] = hist(allrs,100);
semilogy(xx,yy,'r-.')

%% Time course
figure;
bin = 1;  %use for down-sampling
filt = 10;  %filter trackes
p_thr = 90;  %angle threshold
ct = 0;  %count incident concidered
meanp = [];  %statistics for turn rate

tt = zeros(1,length(cand));
for c = 1:length(cand)
    id = cand(c);
    temp2 = Tracks(id).Time;
    tt(c) = temp2(end);
end
Mt = max(tt);

%%%time windows
partt = 1;
bin_times = [0:Mt/(partt):Mt];
trknum = zeros(1,partt);
ccc = hsv(partt);

for bt = 1%:length(bin_times)-1
    
    allas = [];  %angles
    allps = [];  %peurettes (sharp turns)
    allrs = [];  %runs
    allds = [];  %test variables
    allcs = [];  %autocorrelation (half time)

    for c = 1:length(cand)
        id = cand(c);
        temp2 = Tracks(id).Time;
%         plot(temp2); hold on
        overlap = find(temp2>bin_times(bt) & temp2<bin_times(bt+1));
        %%%only analyze overlaping time points
        if length(overlap)>2%isempty(overlap)~=1
            trknum(bt) = trknum(bt)+1;
            
            
            %%%pre-processing
            temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
            temp1 = zeros(length(overlap),2);
            temp1(:,1) = smooth(temp(overlap,1),filt);
            temp1(:,2) = smooth(temp(overlap,2),filt);
            subs = temp1(1:bin:end,:);
            newtime = temp2(1:bin:end);
            vecs = [[0,0]; diff(subs)];
            newtime = newtime(1:size(vecs,1));
            dists = zeros(1,size(subs,1));
            angs = zeros(1,size(vecs,1));
            angT = zeros(1,size(vecs,1));
            
            %%% for distance to target
            for dd = 1:length(dists)
                dists(dd) = distance(subs(dd,:),target);
            end
            %dists(dists>2000) = [];
            
            %%% for angle & angle to target
            for dd = 2:length(angs)
                angT(dd) = angles(vecs(dd-1,:),target-subs(dd-1,:)); %ThetaInDegrees;%
                angs(dd) = angles(vecs(dd-1,:),vecs(dd,:));
            end
            angs(2) = [];
            angT(2) = [];
            

            %%% for "Pirouttes" frequnecy
            [timestamps,runs] = def_turns(angs,p_thr, subs);

            prate = length(timestamps)/(newtime(end)-newtime(1));
            %angs(2) = 0;
            R = corrcoef(angs,angT);
            dist_turn = [];
            
            %%%turning conditions
            if isempty(timestamps)~=1  %&& length(timestamps)<300
                for d = 1:length(timestamps)
                    dist_turn = [dist_turn distance(temp1(timestamps+1,:),target)];
                end
            time_turn = newtime(timestamps);
            [aa,bb] = hist(angT(timestamps),[0:10:180]);
            aa = aa/sum(aa);
%             plot(bb,aa);pause()
            
            tempT = angT(max(timestamps(2:end),1));
            ini_runs = find(tempT<=45);  %??--> picking the angle to target before runs
            
            allds = [allds angT(timestamps(1:end-1))];%angs(ini_runs)];%; aa];newtime];%time_turn];%R(2)];%(length(find(angT>90)))];  %%%turn time, angle, correlation...
            
            allrs = [allrs runs(ini_runs)];  %%%all runs
            allps = [allps prate];  %turn rate!
            allas = [allas angT];  %%%all angles
            
            end
            
            
            %temp = dists(timestamps(1:end-1));
            %ini_runs = find(temp>1000);  %?? try with distance
            
            
%             if length(angT)>140
%             angT(isnan(angT)==1) = 0;
%             [acf,lags,bounds] = autocorr( temp1(:,1),size(temp1,1)-1 );%( angT,length(angT)-1 ); %
%             acf(acf<0) = 0;
%             [val,pos] = min(abs(acf-bounds(1)));%(abs(acf-0.5));
%             allcs = [allcs pos*bin/14];
%             end
            %%%plot(angT);pause();  %visualization
            
            
        end
    end
    
%     figure;
%     hist(real(allds),100)
%     meanp = [meanp; [mean(allds) std(allds) median(allds)]];
%     hist(real(allas),100)
%     hist(real(allrs),100)

%     meanp = [meanp; [mean(allps) std(allps)]];
%     plot(bb,sum(allds)/max(sum(allds)),'color',ccc(bt,:));
%     hold on

    [yy,xx] = hist(allrs,100);
    hist(allrs,100)
    semilogy(xx,yy,'r-.')
    %xlim([0,2])
    
end

%% run vs. angle
Rs = allrs*0.0367;
Ts = allds;
 
cri = 1.5;
pos = find(Rs<cri);
Rs(pos) = [];
Ts(pos) = [];
ang_run = [];

binsz = 18;
bb = 0:binsz:180;
ang_runs = zeros(size(bb));
for bs = 1:length(bb)-1
    poss = find(bb(bs)<Ts & bb(bs+1)>Ts);
    temp = Rs(poss)*0.0367;
    ang_run(bs) = mean(temp);
end

yyaxis left
% plot(allds,allrs*0.0367,'o')
plot(Ts,Rs,'o')
yyaxis right
plot(bb(1:end-1),ang_run,'r-o')

%%
cri = 90;
pos = find(Ts>cri);
op = Rs(pos);
pos = find(Ts<=cri);
sa = Rs(pos);
[yy,xx] = hist(op,100);
semilogy(xx,yy,'r-.')
hold on
[yy,xx] = hist(sa,100);
semilogy(xx,yy,'b-.')
set(gca,'Fontsize',30)


