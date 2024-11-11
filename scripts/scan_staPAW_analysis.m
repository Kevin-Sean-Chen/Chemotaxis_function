%% scan_staPAW_analysis

%% load data
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240224_132726_cv_staPWA.mat') %%% repeat 0_50
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat') %%% Data_salt0_50 %%%
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240324_114402_cv_staPWA.mat') %%% Data_salt100_50
% load('/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/20240402_110348_cv_staPWA.mat') %%% Data_nai and Data_ave
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240520_023653_cv_staPWA.mat')  %%% Data_salt0_50_0513
% %%% Data_salt100_50_0511
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240618_233016_cv_staPWA.mat')

%% mutants!
% load('/projects/LEIFER/Kevin/Data_learn/Mutants/ssm_analysis/20240603_164117_cv_staPWA.mat')  %% AIB
% load('/projects/LEIFER/Kevin/Data_learn/Mutants/ssm_analysis/20240604_031244_cv_staPWA.mat')  %% AIY
% %%% AIZ...

%% ablation for null model control...
neg_control = struct();
neg_control(num_cv, 1, num_repeats).params = [];  % trained parameter
neg_control(num_cv, 1, num_repeats).logp = [];  % train ll trace
neg_control(num_cv, 1, num_repeats).test_ll = []; % test ll
hmm_control = neg_control;  % for state control

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
    for k = 1:1
%         mmhat0 = init_staPAW(k);      %%%%%%%% initialize parameters here %%%%%%%%... include this in EM function for now
    k
    
        % run repeats
        for r = 1:num_repeats
%             [mmhat, logp] = runEM_staPAW(Data_train, k);      %%%%%%%% EM training here %%%%%%%% % check K=1 and run null model separately
            mmhat = all_record(c,1,r).params;
            temp_ws = mmhat.wts;
            temp_ws([2:5]) = zeros(1,4);
            temp_ws([14:17]) = ones(1,4);
%             temp_ws([12:13]) = [mean(temp_ws([12:13])) mean(temp_ws([12:13]))];
            mmhat.wts = temp_ws;  % ad-hoc ablation!
            testLL = staPAW_test(Data_test, mmhat, num_subsamp, len_data);      %%%%%%%% testing here %%%%%%%%
            
            % putting results in
            neg_control(c, 1, r).params = mmhat;  % trained parameter
%             neg_control(c, 1, r).logp = logp;  % train ll trace
            neg_control(c, 1, r).test_ll = testLL; % test ll
            
            % for state control
            mmhat = all_record(c,1,r).params;
            temp_wt_state = mmhat.wts_state;
            temp_wt_state = temp_wt_state*0;
            mmhat.wts_state = temp_wt_state;
            testLL = staPAW_test(Data_test, mmhat, num_subsamp, len_data);      %%%%%%%% testing here %%%%%%%%
            hmm_control(c, 1, r).params = mmhat;  % trained parameter
%             hmm_control(c, 1, r).logp = logp;  % train ll trace
            hmm_control(c, 1, r).test_ll = testLL; % test ll
        end
    end
end

%% state analysis
cvi = 2;  %3
null_ll = [];
state_null = [];
for cc = 1:cvi
    for rr = 1:5
        temp_null = neg_control(cc,1,rr).test_ll;
        null_ll = [null_ll temp_null];
        temp_null = hmm_control(cc,1,rr).test_ll;
        state_null = [state_null temp_null];
    end
end
base_ll = mean(null_ll);
figure;
for cc = 1:cvi
for rr = 1:5
    for kk = 1:4
        temp = all_record(cc,kk,rr).test_ll;
        plot(zeros(1,length(temp))+kk + randn(1, length(temp))*0.07, (temp - base_ll)/log(2), 'k.')
        hold on
    end
    

end
end
plot(randn(1, length(null_ll))*0.07, (null_ll - base_ll)/log(2), 'k.')
set(gcf,'color','w'); set(gca,'Fontsize',20);
ylabel('test LL (bits/s)'); xlabel('models')

plot(randn(1, length(null_ll))*0.07 + 5, (state_null - base_ll)/log(2), 'k.')
set(gcf,'color','w'); set(gca,'Fontsize',20);
ylabel('test LL (bits/s)'); xlabel('models')

%% test LL
figure;
for cc = 3
    for rr = 1:5
        for kk = 3
            temp = all_record(cc,kk,rr).logp;
            plot(temp, 'k-o'); hold on
        end
    end
end
