% Figure 4
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4
% extracted from script 'classify_learn3.m'

%% load data, with fully trained parameters
load(fullfile(datadir,'data4plots', 'classify_learn3_vars7.mat'));

%% plotting predictions
rtime = squeeze(mean(mean(mean(data_len,1),2),4))*5/14/60/60;
figure;
plot(rtime, mean(mean(cv_perf,3)),'k-o')
hold on
errorbar(rtime, mean(mean(cv_perf,3)), std(mean(cv_perf,3))/sqrt(K), '.k')

plot(rtime, mean(mean(cv_perf_bin,3)),'g--')
errorbar(rtime, mean(mean(cv_perf_bin,3)), std(mean(cv_perf_bin,3))./sqrt(K.*scal_vec), '.g')

plot(rtime, mean((cv_nbc),2),'r--')
errorbar(rtime, mean((cv_nbc),2), std((cv_nbc),0,2)'./sqrt(K.*scal_vec), '.r')

plot(rtime, mean(mean(cv_perf_null,3)),'b-*')
errorbar(rtime, mean(mean(cv_perf_null,3)), std(mean(cv_perf_null,3))./sqrt(K.*scal_vec), '*b')

xlabel('mean data length (hr)')
ylabel('cross-validation (ratio)')
legend({'mean prediction','','binomial prediction', '','\Delta C', '','behavior',''})
set(gcf,'color','w'); set(gca,'Fontsize',20);

figure
for ii = 1:3
    temp_cv = squeeze(cv_perf(ii,:,:));
    plot(rtime, mean(temp_cv,2),'-o')
    hold on
    errorbar(rtime, mean(temp_cv,2), std(temp_cv,[],2)/sqrt(K), '.k')
end
legend({'appetitive','','naive','','aversive',''})
xlabel('mean data length (hr)')
ylabel('cross-validation (ratio)')
set(gcf,'color','w'); set(gca,'Fontsize',20);

%% some post analysis for variability
ttl = {'appetitive','naive','aversive'};
figure
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
for cc = 1:3
    subplot(1,3,cc);
for ii = 1:K
    temp = squeeze(mle_params(ii,cc,3:6))'*cosBasis';
    plot(tt,temp)
    hold on
end; title(ttl{cc})
end

figure; 
xx = 0:length(cosBasis)-1;
for cc = 1:3
    subplot(1,3,cc)
for ii = 1:K
    amp = mle_params(ii,cc,8); tau = mle_params(ii,cc,9); temp = (-amp*exp(-xx/tau));
    plot(tt,temp)
    hold on
end; title(ttl{cc})
end

%% information analysis
figure
col = {'b','k','r'};
tils = {'full','K_c','K_{c^\perp}'};
for nn = 1:3
    subplot(1,3,nn)
    temp_m = squeeze(mean(testLL(:,:,nn),1));
    temp_s = squeeze(std(testLL(:,:,nn),0,1))/sqrt(K);
    for cc = 1:3
        bar(cc, temp_m(cc), 'FaceColor',col{cc})
        hold on
        errorbar(cc, temp_m(cc), temp_s(cc), 'k.');
    end
    title(tils{nn})
    xticklabels([]);
end
