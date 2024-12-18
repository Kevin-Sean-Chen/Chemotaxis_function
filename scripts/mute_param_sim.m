% mute_param_sim
%%%% this is a script to load in the dPAW fit to all strains and simulate
%%%% trajectories to comapre chemotaxis index to the measurements

%% load
%%% fitted variables
load('/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/mutant_gof_vars_n2.mat')

%%% mutant Data
data_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/Data_files';
fileList = dir(fullfile(data_path, '*.mat'));

save2path = '/projects/LEIFER/Kevin/Publications/Chen_learning_2023/revision/mutant_sim_test/';

%%% N2 Data
data_n2 = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat'};
     
%% prep
rng(42)
clear specs
specs = struct();

Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
M = Cmap.vq1;
M = fliplr(flipud(M));
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);

specs.M = M;
specs.fr = 14/5;
specs.cosBasis = cosBasis;
specs.T = floor(30*60*14/5);
specs.dt = 1;
specs.REP = 100;

data_ci = zeros(1,length(fileList));
sim_ci = zeros(1,length(fileList));

%% looping to compute simulation and data
for li = 1:numel(fileList)
    
    li
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    ids = all_ids{li+3};  % retrieve training tracks
    
    %%% record data
    data_i = Data(ids==1);
    cii = Data2CI(data_i);
    data_ci(li) = cii;
    
    %%% record simulation
    xi = mle_mut_params(li,:);
    v_dist = post_compute_dis(data_i);
    [tracks, CI] = param2tracks_mute(xi, specs, v_dist, data_i);
    sim_ci(li) = CI;
    
    %%% save tracks
    filename4track =  sprintf('%stracks_%d.mat',save2path,li);
    save(filename4track, 'tracks'); %% important to comment out to prevent replacement! use it to bootstrap later
end

%% plotting
figure
plot(data_ci(:), sim_ci(:),'o')
xlim([0 1]); ylim([0 1])

figure;
% cols = ['b','r','k'];
subplot(211); 
b1=bar(reshape(data_ci, 3, 5)', 'grouped');
subplot(212); 
% hold on
b2=bar(reshape(sim_ci, 3, 5)', 'grouped');

%% for N2!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(13);

data_ci_n2 = zeros(1,3);
sim_ci_n2 = zeros(1,3);

specs.REP = 100;
%% looping to compute simulation and data
for li = 1:3
    
    li
    %%% load Data structure
    load(data_n2{li});
    ids = all_ids{li};  % retrieve training tracks
    
    %%% record data
    data_i = Data;%(ids==1);
    cii = Data2CI(data_i);
    data_ci_n2(li) = cii;
    
    %%% record simulation
    xi = mle_params(li,:);
    v_dist = post_compute_dis(data_i);
    [N2tracks, CI] = param2tracks_mute(xi, specs, v_dist, data_i);
    sim_ci_n2(li) = CI;
    
    %%% save tracks
    filename4track =  sprintf('%sN2tracks_%d.mat',save2path,li);
    save(filename4track, 'N2tracks'); % important to comment out to prevent replacement! use it to bootstrap later
end

%% plotting
figure
plot(data_ci_n2(:), sim_ci_n2(:),'o')
xlim([0 1]); ylim([0 1])

figure;
subplot(211); 
b1=bar(reshape(data_ci_n2, 3, 1)', 'grouped');
subplot(212); 
b2=bar(reshape(sim_ci_n2, 3, 1)', 'grouped');

%% functions
function [CI] = Data2CI(Data)
    %%% CI
    cii = 0;
    for ii = 1:length(Data)
        temp = Data(ii).dc;
        if temp(end)>temp(1)
            cii = cii + 1;
        end
    end
    CI = 2*cii/length(Data) - 1;
    
end

function alldis = post_compute_dis(Data)
    allxys = [];
    for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end
    [xxf, yyf, alltrials, time] = data2xy(Data);
    alldis = vecnorm(diff(allxys,1,2));
    alldis(alldis>10) = [];
end