% model_compare_corr
% given the cross-validated models with different states,
% compare dPAW with staPAW for the autocorrelation sturcture

%% load trained models
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat')
% rng(42) %37

%% asign models
rr = 1;  %1,3
dPAW_fit = all_record(rr,1,1).params;  % single state
staPAW_fit = all_record(rr,2,1).params;  % two-state model

%% configure data targets
[xxf, yyf, alltrials, alltime] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:150000];
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);
maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(staPAW_fit,xx,yy,mask);

%% show time series 
basis = dPAW_fit.basis;
wind = length(basis);
dtheta_data = yy(1,:);  % headings
dr_data = yy(2,:);  % speed

%% model prediction
[dth_pred_dPAW, dr_pred_dPAW] = model_predit_PAW(dPAW_fit, xx, yy, gams_, 1);
[dth_pred_staPAW, dr_pred_staPAW] = model_predit_PAW(staPAW_fit, xx, yy, gams_, 1);

%% analyze <dr,dr'>
nlag = 60;
% dv = dr_data(find(dr_data<100));
[As_data,lgs] = autocorr(dr_data, nlag);
[As_dpaw,lgs] = autocorr(dr_pred_dPAW, nlag);
[As_stapaw,lgs] = autocorr(dr_pred_staPAW, nlag);

figure;
plot(lgs*5/14, As_data, 'k-', 'LineWidth',3); hold on
plot(lgs*5/14, As_stapaw, 'LineWidth',3);
plot(lgs*5/14, As_dpaw, '--', 'LineWidth',3);
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('<drdr''>')

%% further analysis of ITI!
pos_data = find(abs(dtheta_data)>50); %.7
pos_dpaw = find(abs(dth_pred_dPAW)>50);
pos_stapaw = find(abs(dth_pred_staPAW)>50);

%%% model prediction
figure; 
[iti_cnt, iti_bin] = iti2hist(pos_data);
plot(iti_bin, iti_cnt, 'LineWidth',3); hold on;
[iti_cnt, iti_bin] = iti2hist(pos_dpaw);
plot(iti_bin, iti_cnt, '--','LineWidth',3);
[iti_cnt, iti_bin] = iti2hist(pos_stapaw);
plot(iti_bin, iti_cnt, 'LineWidth',3);
% plot(iti_bin, iti_cnt, 'k-','MarkerSize',15); 
set(gcf,'color','w'); set(gca,'Fontsize',20);  ylabel('probability');  xlabel('turn interval (s)'); xlim([0, 70]); set(gca, 'YScale', 'log')

%% 
function [iti_cnt, iti_bin] = iti2hist(pos)
    bins = 0:0.7:70;
    iti = diff(pos)*5/14;
    H1 = histogram(iti, bins, 'Normalization', 'pdf'); H1.Visible = 'off'; %,
    iti_cnt = H1.Values;
    iti_bin = H1.BinEdges(1:end-1) + diff(H1.BinEdges)/2;
end