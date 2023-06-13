%%% Path_triggered
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothSpeed','AngSpeed','SmoothX','SmoothY','Behaviors', 'LEDVoltages'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% simple triggering analysis
poly_degree = 3;  %polynomial fitting for moving window
filt = 7;  %window the path (has to be odd because it is +/- points around the center)
fr = 1/14;  %1/14 seconds between each frame  (~0.0714 second for each frame)
bin = 5;  %down-sampling
nn = length(Tracks); %number of worms selected
windt = 14 * (1/(bin*fr)); %time window in seconds
acst = 4 * (1/(bin*fr));  % acausal window
trig_angs = [];

figure
for ii = 1:nn
    
    %%% calculate angles
    temp = Tracks(ii);  % for loading saved tracks
    temp1 = zeros(round(size(temp.Path,1)/1),2);
    temp1(:,1) = smooth(temp.SmoothX, filt,'sgolay',poly_degree); %smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp.SmoothY, filt,'sgolay',poly_degree); %smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    vecs = diff(subs);
    vecs = [vecs; vecs(end,:)];   % compensating for the size change/time shift after taking difference 
    angs = zeros(1,size(vecs,1));    
    
    %%% iterate through worms
    for dd = 2:length(angs)
        %%% angle function
        angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
    end
    
    %%% LED signal
    stim = Tracks(ii).LEDVoltages;
    stim = stim(1:bin:end);
    pos = find(diff(stim)>1);% & diff(stim)>1);  %0.01 1 3 thresholding to find impulse
    pos = randi([1,length(stim)],1,length(pos)*5);%randi(length(stim));
    if isempty(pos)~=1
        for pp = 1:length(pos)
            if (pos+windt)<length(stim) & pos-acst>1
                wind = floor((pos(pp)-acst):(pos(pp)+windt));
                trig_angs = [trig_angs; angs(wind)];
             
%                 plot(subs(wind,1),subs(wind,2)); hold on
            end
        end
    end
    
end

figure;
plot([-acst:windt]*((bin*fr)), mean(abs(trig_angs)) / (fr*bin))

%%
figure
patch([0 5 5 0], [15 15, 55 55], [0.8 0.8 0.8])
hold on
plot([-acst:windt]*((bin*fr)), mean(abs(trig_angs)) / (fr*bin))
xlabel('time (s)')
ylabel('angular speed (|deg|/s)')