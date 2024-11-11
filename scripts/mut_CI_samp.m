% mut_CI_samp
%%% given the simulated tracks from "mute_param_sim" and the training data,
%%% we bootstrap to copute the distribution of these CIs

%% load data and simulation
%%% fitted variables
load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/mutant_gof_vars_n2.mat')
% load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/mute_sim_vars.mat')

%%% mutant Data
data_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/Data_files';
fileList = dir(fullfile(data_path, '*.mat'));

mut_sim_path = '/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/mutant_sim/';
mut_List = dir(fullfile(mut_sim_path, '*.mat'));

%%% N2 Data
data_n2 = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat'};

rng(37) %37

%% setup
reps = 50;  % repeated samples
% n_samps = 30;  % number of samples  % limited by min tracsk... may want to sampled adaptively
frac_samp = 1/2.5;  % now as fraction

ci_samp_data = zeros(length(fileList), reps);
ci_samp_sims = zeros(length(fileList), reps);

ci_samp_data_n2 = zeros(3, reps);
ci_samp_sims_n2 = zeros(3, reps);

%% loop for data and sim tracks
for li = 1:numel(fileList)
    
    li
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    ids = all_ids{li+3};  % retrieve training tracks
    
    %%% record data
    data_i = Data(ids==1);
    n_samps = floor(frac_samp*length(data_i));
    for rr = 1:reps
        temp = randperm(length(data_i));
        data_ir = data_i(temp(1:n_samps));
    	ci_samp_data(li,rr) = Data2CI(data_ir);
    end
    
    %%% load sim tracks
%     fileName = mut_List(li).name;
    fileName = sprintf('tracks_%d.mat', li);
    filePath = fullfile(mut_sim_path, fileName);
    load(filePath);
    
    %%% record data 
    n_samps = floor(frac_samp*length(data_i));
    for rr = 1:reps
        temp = randperm(length(tracks));
        track_ir = tracks(temp(1:n_samps));
    	ci_samp_sims(li,rr) = Data2CI(track_ir);
    end
end

%% for N2
rng(0) %42

for li = 1:3
    
    li
    %%% load Data structure
    load(data_n2{li});
    ids = all_ids{li};  % retrieve training tracks
    
    %%% record data
    data_i = Data(ids==1);
    if li==1
        data_i = Data(347:end);
    end
    n_samps = floor(frac_samp*length(data_i));
    for rr = 1:reps
        temp = randperm(length(data_i));
        data_ir = data_i(temp(1:n_samps));
    	ci_samp_data_n2(li,rr) = Data2CI(data_ir);
    end
    
    %%% load sim tracks
    fileName = sprintf('N2tracks_%d.mat', li);
    filePath = fullfile(mut_sim_path, fileName);
    load(filePath);
    
    %%% record data 
    n_samps = floor(frac_samp*length(data_i)/1);
    for rr = 1:reps
        temp = randperm(length(N2tracks));
        track_ir = N2tracks(temp(1:n_samps));
    	ci_samp_sims_n2(li,rr) = Data2CI(track_ir);
    end
end

%%
% fileName = sprintf('N2tracks_%d.mat', 2);
% filePath = fullfile(mut_sim_path, fileName);
% load(filePath);
% Data2CI(N2tracks)

%% plottng
x = mean(ci_samp_data,2);
y = mean(ci_samp_sims,2);
xerror = std(ci_samp_data,1,2);
yerror = std(ci_samp_sims,1,2);
figure;
scatter(x, y, 'b', 'filled');
hold on
errorbar(x, y, yerror, 'bo', 'LineStyle', 'none');
errorbar(x, y, xerror, 'horizontal', 'bo', 'LineStyle', 'none');
plot([-.1,1.1],[-.1,1.1],'r--')

for i = 1:length(fileList)
    temp = fileList(i).name;
    pattern = 'Data_(.*)\.mat';
    tokens = regexp(temp, pattern, 'tokens');
    modifiedString = strrep(tokens{1}, '_', '-');
    text(x(i) + 0.02, y(i) - 0.02, modifiedString, 'FontSize', 12, 'Color', 'k');
end
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('CI data'); ylabel('CI dPAW')

correlationMatrix = corr(x(:), y(:))

%%% hold on N2
hold on
x = mean(ci_samp_data_n2,2);
y = mean(ci_samp_sims_n2,2);
xerror = std(ci_samp_data_n2,1,2);
yerror = std(ci_samp_sims_n2,1,2);

scatter(x, y, 'k', 'filled');
hold on
errorbar(x, y, yerror, 'ko', 'LineStyle', 'none');
errorbar(x, y, xerror, 'horizontal', 'ko', 'LineStyle', 'none');

y

%% function
function [CI] = Data2CI(Data)
    %%% CI
    cii = 0;
    for ii = 1:length(Data)
        temp = Data(ii).dc;
        if temp(end)>temp(1)+0
            cii = cii + 1;
        end
    end
    CI = 2*cii/length(Data) - 1;
end