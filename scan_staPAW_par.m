% scan_staPAW
% this script is used to scan through test-sets, n-states, and repeats
% should run on the cluster given the heavy loading
% make sure that the loaded train/test data and functional fits are correct

%% load path and data
% set path to reach the subfolder functions
addpath(genpath(pwd))   %%%% use this for cluster computation!  %%%

filepath = '/projects/LEIFER/Kevin/Data_salt/data_analysis/';
load([filepath,'Data_salt100_50.mat']);     %%%%%%%% load Data structure here %%%%%%%%

%% scanning parameters
num_cv = 3;  % half is enough data
num_states = 4;  % number of states scanned
num_repeats = 5;  % repeat for EM selection
num_subsamp = 20;  % testing subsamples
len_data = 10000;  % test-data length for the subsamples... can turn to track-based too

data_id = 1:length(Data);  % original track ID
data_sv = crossvalind('Kfold',data_id, num_cv);

% important structure that carries all info
all_record = struct();
all_record(num_cv, num_states, num_repeats).params = [];  % trained parameter
all_record(num_cv, num_states, num_repeats).logp = [];  % train ll trace
all_record(num_cv, num_states, num_repeats).test_ll = []; % test ll

%% parallel computing
% Open a parallel pool
parpool;%('local', 16); % Opens a pool with 16 workers

%%
% loop through train-test sets
for c = 1:num_cv
    % loading train-test Data tracks
%     test_set = (data_sv==c);
%     train_set = ~test_set;
    train_set = (data_sv==c);  % tweek here to speed up...
    test_set = ~train_set;
    
    Data_test = Data(test_set);
    Data_train = Data(train_set);
    c
    
    % scan through states
    for k = 1:num_states
%         mmhat0 = init_staPAW(k);      %%%%%%%% initialize parameters here %%%%%%%%... include this in EM function for now
    k
    
        % run repeats
        parfor r = 1:num_repeats
            [mmhat, logp] = runEM_staPAW(Data_train, k);      %%%%%%%% EM training here %%%%%%%% % check K=1 and run null model separately
            testLL = staPAW_test(Data_test, mmhat, num_subsamp, len_data);      %%%%%%%% testing here %%%%%%%%
            
            % putting results in
            all_record(c, k, r).params = mmhat;  % trained parameter
            all_record(c, k, r).logp = logp;  % train ll trace
            all_record(c, k, r).test_ll = testLL; % test ll
        end
    end
end

%% close paralllel computing
delete(gcp)

%% saving
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');  % in case of overwriting the large data...
timestampString = datestr(currentDateTime, 'yyyymmdd_HHMMSS');
filename = [filepath, timestampString, '_cv_staPWA.mat']; % save all variable... with caution~
save(filename);
