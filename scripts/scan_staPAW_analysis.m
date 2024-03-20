%% scan_staPAW_analysis

%% load data
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240224_132726_cv_staPWA.mat')
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat')
%% ablation for null model control...
neg_control = struct();
neg_control(num_cv, 1, num_repeats).params = [];  % trained parameter
neg_control(num_cv, 1, num_repeats).logp = [];  % train ll trace
neg_control(num_cv, 1, num_repeats).test_ll = []; % test ll

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
            neg_control(c, 1, r).logp = logp;  % train ll trace
            neg_control(c, 1, r).test_ll = testLL; % test ll
        end
    end
end

%% state analysis
cvi = 3;
null_ll = [];
for cc = 1:3
    for rr = 1:5
        temp_null = neg_control(cc,1,rr).test_ll;
        null_ll = [null_ll temp_null];
    end
end
base_ll = mean(null_ll);
figure;
for cc = 1:3
for rr = 1:5
    for kk = 1:4
        temp = all_record(cc,kk,rr).test_ll;
        plot(zeros(1,length(temp))+kk + randn(1, length(temp))*0.07, temp - base_ll, 'k.')
        hold on
    end
    

end
end
plot(randn(1, length(null_ll))*0.07, null_ll - base_ll, 'k.')
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
