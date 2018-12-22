%Track_analysis_plot
clear
clc
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%for batch analysis
fields_to_load = {'Path','Time'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%%% criteria %%%
nn = length(Tracks); %number of worms selected
mint = 60*1; %minimum time in seconds
disth = 300;  %radius of pixels from target
target = [2517,975];%[950,1100];%  %position of target/sourse of odorant (approximated from images)
center = [1353 975];  %poition of the starting point in the middle of plate (Ali's parameter)

%visualize all paths and check criteria
cand = [];
for i = 1:nn
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        plot(Tracks(i).Path(:,1),Tracks(i).Path(:,2))
        hold on
        cand = [cand i];
    end
end

%% tracks and distribution
cand = [];
alltr = [];
t_window = [2 30]*60;
mint = 60;

for i = 1:nn%/1.5
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        
        temp2 = Tracks(i).Time;
        overlap = find(temp2>t_window(1) & temp2<t_window(2));
        %%%only analyze overlaping time points
        if length(overlap) > 2
            temp = Tracks(i).Path;%smooth(Tracks(id).Path,6);
            temp1 = zeros(length(overlap),2);
            temp1(:,1) = temp(overlap,1)-center(1);
            temp1(:,2) = temp(overlap,2)-center(2);
            
            plot(temp1(:,1),temp1(:,2))
            hold on
            cand = [cand i];
            alltr = [alltr temp1(:,1)'*0.0367];
        end
        
    end
end

%% alignment
samp = 8000;
cand = [];
alltr = [];
dirs = [];
t_window = [2 20]*60;
mint = 120;

for i = 1:samp
    if Tracks(i).Time(end)-Tracks(i).Time(1) > mint
        
        %i = randi([1,nn]);
        temp2 = Tracks(i).Time;
        overlap = find(temp2>t_window(1) & temp2<t_window(2));
        %%%only analyze overlaping time points
        if length(overlap) > 2
            temp = Tracks(i).Path;%smooth(Tracks(id).Path,6);
            temp1 = zeros(length(overlap),2);
            temp1(:,1) = temp(overlap,1)-temp(overlap(1),1);
            temp1(:,2) = temp(overlap,2)-temp(overlap(1),2);
            
            plot(temp1(:,1),temp1(:,2))
            hold on
            cand = [cand i];
            alltr = [alltr temp1(:,1)'];
            if temp1(end,2)>0; tempa = angles(temp1(end,:),[1,0]); else; tempa = 360-angles(temp1(end,:),[1,0]); end
%             tempa = angles(temp1(end,:),[1,0]);
            dirs = [ dirs  tempa];
            
        end
        
    end
end


xlim([-1500 1500]);
ylim([-1000 1000]);

figure;
[x,y] = hist(dirs,50);
x = [x x(1)];
x = x/sum(x);
polarplot(0:2*pi/50:2*pi,x)

