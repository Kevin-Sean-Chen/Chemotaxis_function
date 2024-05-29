% model_comparison
% given the cross-validated models with different states,
% compare dPAW with staPAW with control models
% this is done by plotting a time series of everything but headings, then
% look at the heading prediction, and further quantify the performance with
% some metric...

%% load trained models
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat')
rng(42) %37

%% asign models
rr = 1;
dPAW_fit = all_record(rr,1,1).params;  % single state
staPAW_fit = all_record(rr,2,1).params;  % two-state model
temp_ws = dPAW_fit.wts;
temp_ws([2:5]) = zeros(1,4);
temp_ws([14:17]) = ones(1,4);
null_fit = dPAW_fit;
null_fit.wts = temp_ws;  % ad-hoc ablation!

%% configure data targets
[xxf, yyf, alltrials, alltime] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1100:1600];%[3000:3500];
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);
maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(staPAW_fit,xx,yy,mask);

%% show time series 
basis = dPAW_fit.basis;
wind = length(basis);
plt_states = gams_(:,wind:end)';  % states
plt_dtheta = yy(1,wind:end);  % headings
plt_dr = yy(2,wind:end);  % speed

%% model prediction
[dth_pred_null, dr_pred_null] = model_predit_PAW(null_fit, xx, yy, gams_);
[dth_pred_dPAW, dr_pred_dPAW] = model_predit_PAW(dPAW_fit, xx, yy, gams_);
[dth_pred_staPAW, dr_pred_staPAW] = model_predit_PAW(staPAW_fit, xx, yy, gams_);

%% plot for dtheta
time = [1:length(dth_pred_null)]*5/14;
figure()
subplot(5,1,1)
plot(time, plt_states); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); ylabel('P(Z)'); ylim([-0.05, 1.05]); xlim([0,max(time)])
title('d\theta')
subplot(5,1,2)
plot(time, plt_dtheta); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); xlim([0,max(time)]); ylabel('data')
subplot(5,1,3)
plot(time, dth_pred_staPAW); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); xlim([0,max(time)]); ylabel('staPAW')
subplot(5,1,4)
plot(time, dth_pred_dPAW); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); xlim([0,max(time)]); ylabel('dPAW')
subplot(5,1,5)
plot(time, dth_pred_null); set(gcf,'color','w'); set(gca,'Fontsize',15); xlabel('time (s)'); xlim([0,max(time)]); ylabel('control')

%%  plot for dr
dr_factor = 1/30;
figure()
subplot(5,1,1)
plot(time, plt_states); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); ylabel('P(Z)'); ylim([-0.05, 1.05]); xlim([0,max(time)])
title('dr')
subplot(5,1,2)
plot(time, plt_dr*dr_factor); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); xlim([0,max(time)]); ylabel('data')
subplot(5,1,3)
plot(time, dr_pred_staPAW*dr_factor); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); xlim([0,max(time)]); ylabel('staPAW')
subplot(5,1,4)
plot(time, dr_pred_dPAW*dr_factor); xticks([]);set(gcf,'color','w'); set(gca,'Fontsize',15); xlim([0,max(time)]); ylabel('dPAW')
subplot(5,1,5)
dr_pred_null = plt_dr(randperm(length(plt_dr)))*dr_factor;
plot(time, dr_pred_null); set(gcf,'color','w'); set(gca,'Fontsize',15); xlabel('time (s)'); xlim([0,max(time)]); ylabel('control')

%% repeat for data statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(43)

%% statsitics
samp_len = 1000;
n_samps = 70;
threshold = 70;  % angle threshold
pred_power = zeros(3, n_samps);
pred_dr = zeros(3, n_samps);

for ss = 1:n_samps
    %%% sampling
    randomInteger = randi([1, length(xxf)-samp_len]);
    wind_i = randomInteger:randomInteger+samp_len;
    xxi = xxf(:,wind_i);
    yyi = yyf(:,wind_i);
    maski = maskf(wind_i);
    [logp_test,gams_i,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(staPAW_fit,xxi,yyi,maski);
    
    %%% dth_vectors
    ith_data = yyi(1,wind:end);
    ith_run = yyi(2,wind:end);
    [ith_pred_null, ith_dr_null] = model_predit_PAW(null_fit, xxi, yyi, gams_i);
    [ith_pred_dPAW, ith_dr_dPAW] = model_predit_PAW(dPAW_fit, xxi, yyi, gams_i);
    [ith_pred_staPAW, ith_dr_staPAW] = model_predit_PAW(staPAW_fit, xxi, yyi, gams_i);
    
    %%% measure predictions
    pred_power(1,ss) = predict_turns(ith_data, ith_pred_staPAW, threshold);
    pred_power(2,ss) = predict_turns(ith_data, ith_pred_dPAW, threshold);
    pred_power(3,ss) = predict_turns(ith_data, ith_pred_null, threshold);
    pred_dr(1,ss) = predict_runs(ith_run, ith_dr_staPAW);
    pred_dr(2,ss) = predict_runs(ith_run, ith_dr_dPAW);
    pred_dr(3,ss) = predict_runs(ith_run, ith_dr_null);
end

%% bar for dtheta
figure;
bar([1,2,3], mean(pred_power')); hold on; 
errorbar([1,2,3],mean(pred_power'),std(pred_power')/n_samps^0.5,'ko')
custom_labels = {'staPAW', 'dPAW', 'control'};
xticks([1,2,3]);
xticklabels(custom_labels); set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('turn prediction')
% set(gca, 'YScale', 'log');

%% bar for dr
figure;
bar([1,2,3], mean(pred_dr')); hold on; 
errorbar([1,2,3],mean(pred_dr'),std(pred_dr')/n_samps^0.5,'ko')
custom_labels = {'staPAW', 'dPAW', 'control'};
xticks([1,2,3]);
xticklabels(custom_labels); set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('speed prediction')
% set(gca, 'YScale', 'log');

%%
function [pp] = predict_turns(v1,v2,threshold)
    v1b = zeros(1,length(v1));
    v2b = zeros(1,length(v2));
    v1b(find(abs(v1)>threshold)) = 1;
    v2b(find(abs(v2)>threshold)) = 1;
    
    %%% performance options
    % direct correlation
%     pp = sum(v1b.*v2b)/(sum(v1b)*sum(v2b));  %%% some metric here...
    % Jarrad similarity
%     windw = ones(1,5);
%     v1b = conv(v1b, windw, 'same');
%     v2b = conv(v2b, windw, 'same');
    intersection = sum(v1b & v2b);
    union = sum(v1b | v2b);
    pp = intersection / union;
    % convolved correlation
%     window = ones(1,5);
%     correlation_coefficient = corrcoef(conv(v1b, window, 'same'), conv(v2b, window, 'same'));
%     pp = correlation_coefficient(1, 2);
end
function [pp] = predict_runs(v1,v2)
    windw = ones(1,28); %14
    correlation_coefficient = corrcoef(conv(v1, windw, 'same'), conv(v2, windw, 'same'));%corrcoef(v1, v2);
    pp = correlation_coefficient(1, 2);
end