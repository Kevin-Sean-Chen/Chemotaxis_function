% index_salt
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

%% landscape
rows=2500; cols=3000;
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
gradient_x = x_ * 1;
M0 = (y_*0+1) .* gradient_x;
[x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M100 = (y_*0+1) .* gradient_x;

[ci_, brw_index, wv_index] = compute_index(Tracks(1:end), M0, 15)

%% looping files and conditions
rng(123)
% Tracks = loadtracks(folder_names{1},fields_to_load);
temp = load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/salt_data_folder.mat');
folder_all = {temp.salt_data_50_0, temp.salt_data_50_100};
track_learn_all = cell(1,2);
BWCs = cell(1,2);
cols = ['b','k','r'];
rep_samp = 10;
n_samps = 500;
figure
for cc = 1:2
    folder_names = folder_all{cc};
    BWC = zeros(length(folder_names), 3);
    sub_learn_tracks = cell(1,length(folder_names));
    Tracks = loadtracks(folder_names,fields_to_load);
    for ff = 1:rep_samp%length(folder_names)
%         Tracks = loadtracks(folder_names{ff},fields_to_load);
        samp_id = randperm(length(Tracks));
        Tracks_i = Tracks(samp_id(1:n_samps));
        if cc == 1
            [ci_, brw_index, wv_index] = compute_index(Tracks_i, M0, 20); %20
            BWC(ff,:) = [ci_, brw_index, wv_index];
        else
            [ci_, brw_index, wv_index] = compute_index(Tracks_i, M100, 15);
            BWC(ff,:) = [ci_, brw_index, -wv_index];
        end
        
        sub_learn_tracks{ff} = Tracks;
    end
    
    plot(BWC',cols(cc)); hold on
    BWCs{cc} = BWC;
    track_learn_all{cc} = sub_learn_tracks;
end

%% bar plot
figure
m_0 = mean(BWCs{1},1);
m_100 = mean(BWCs{2},1);
s_0 = std(BWCs{1},[],1)/sqrt(rep_samp);
s_100 = std(BWCs{2},[],1)/sqrt(rep_samp);

hBar = bar([m_0; m_100]');

for k1 = 1:2
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');     
    ydt(k1,:) = hBar(k1).YData;  hold on
%     plot(ctr(k1,:), data_cell{k1},'ko')
end
% set(hBar, {'DisplayName'}, {'App','Naive','Ave'}')
hold on
% errorbar(ctr, ydt, [s_ap;s_na;s_av]', '.r')  
errorbar(ctr, ydt, [m_0;m_100]*0,[s_0;s_100], '.k')         
hold off
ylabel('Index')
% set(gca,'linewidth',2,'FontSize',20)
names = {'CI'; 'BRW'; 'WV'};
set(gca,'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
legend([hBar(1), hBar(2)], '0-50','50-100')
ylim([-1,1])
