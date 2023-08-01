% mVM_analysis
%%% load Data and fitted parameters for analysis
%%% this includes plotting kernels, variances, decision function, and
%%% marginal turning rates

%% load Data and params
datas = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test2.mat'};

params_mat = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param4.mat');

%% assign mle parameters
time_scale = 5/14;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
t_vec = 1:length(cosBasis)-1;
% param_mle = params_mat.mle_params;
K1_mle = median(squeeze(param_mle(:,:,1)), 1);
K2_mle = median(squeeze(param_mle(:,:,12)), 1);
A_mle = median(squeeze(param_mle(:,:,2)), 1);
C_mle = median(squeeze(param_mle(:,:,7)), 1);
Ah_mle = median(squeeze(param_mle(:,:,10)), 1);
th_mle = median(squeeze(param_mle(:,:,11)), 1);
alpha_mle = squeeze(median((param_mle(:,:,3:6)), 1))';
b_mle = median(squeeze(param_mle(:,:,13)), 1);
Kc_mle = zeros(3, size(cosBasis,1));
for ii = 1:3
    Kc_mle(ii,:) = alpha_mle(:,ii)'* cosBasis';
end

%% simple plots
% learn_bar_plot(K1_mle, '\kappa_{wv}')
% learn_bar_plot(K2_mle, '\kappa_{pr}')
learn_bar_plot(C_mle*time_scale, 'baseline turn rate (per s)')

%% plotting turn rates
total_turns = zeros(1,3);

for ii = 1:3
    %%% load data
    load(datas{ii});
    [xx_train, yy_train, mask_train] = data2xy(Data(1:400));
    ddc_fit = xx_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
    ang_fit = yy_train;
    
    %%% calculation of turn
    filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
    K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
    dc_dth = filt_ddc + filt_dth*1;
    Pturns = (A_mle(ii)-C_mle(ii))./ (1 + exp( -(dc_dth) )) + C_mle(ii);
    Pturns = Pturns*time_scale;
    total_turns(ii) = nansum(Pturns) / (length(Pturns)-sum(isnan(Pturns)));
end

learn_bar_plot(total_turns, 'turn rate (per s)')

%% multi-plots
figure
subplot(121)
learn_bar_plot(K1_mle, '\kappa_{wv}',0)
subplot(122)
learn_bar_plot(C_mle*time_scale, 'min turn rate (per s)',0)
% subplot(141)
% learn_bar_plot(K1_mle, '\kappa_{wv}',0)
% subplot(142)
% learn_bar_plot(K2_mle, '\kappa_{brw}',0)
% subplot(143)
% learn_bar_plot(A_mle*time_scale, 'max turn rate (per s)',0)
% subplot(144)
% learn_bar_plot(C_mle*time_scale, 'min turn rate (per s)',0)

%% histogram and turning curves
total_turns = zeros(1,3);
numBins = 50;
cols = ['b', 'k', 'r'];
figure
for ii = 1:3
    %%% load data
    load(datas{ii});
    [xx_train, yy_train, mask_train] = data2xy(Data(1:400));
    ddc_fit = xx_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
    ang_fit = yy_train;
    
    %%% calculation of turn
    filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
    K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
%     dc_dth = (filt_ddc + b_mle(ii)*0 + filt_dth)/norm(Kc_mle(ii, :));
    dc_dth = filt_ddc/norm(Kc_mle(ii, :)) + filt_dth;
%     if ii ==1 
%         dc_dth = -( filt_ddc/norm(Kc_mle(ii, :)) + 1*(filt_dth) );
%     else
%         dc_dth = ( filt_ddc/norm(Kc_mle(ii, :)) + 1*(filt_dth) );%/norm(K_h_rec)*1;
%     end

    subplot(212)
%     histogram(dc_dth + mean(filt_dth)*0, 'FaceColor', cols(ii),'Normalization', 'probability'); hold on
    histogram(dc_dth, numBins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', cols(ii), 'FaceAlpha', 0.7); hold on
    subplot(211)
    xx_cond = linspace(-180,80,100); %linspace(min(dc_dth),max(dc_dth),100);
    plot(xx_cond, ( (A_mle(ii)-C_mle(ii))./(1+exp(-(xx_cond*norm(Kc_mle(ii, :)) + 1*median(filt_dth))))+C_mle(ii) )*time_scale,'Color',cols(ii)); hold on
%     scatter(dc_dth + 0*(filt_dth), ( (A_mle(ii)-C_mle(ii))./(1+exp(-(dc_dth*norm(Kc_mle(ii, :)) + 0*(filt_dth) )))+C_mle(ii) )*time_scale,'Marker', '.','MarkerFaceColor',cols(ii),'MarkerFaceAlpha', 0.5); hold on
%     dc_dth_scaled = filt_ddc + filt_dth;
%     plot(dc_dth, (A_mle(ii)-C_mle(ii))./(1+exp(dc_dth_scaled))+C_mle(ii),'.'); hold on
    
end

%% emperical vs. model turn probability
model_Pturn = zeros(1,3);
exp_Pturn = zeros(1,3);
turn_thr = 130;
gamm = 0.25;
d2r = pi/180;

figure
for ii = 1:3
    %%% load data
    load(datas{ii});
    [xx_train, yy_train, mask_train] = data2xy(Data(1:400));
    ddc_fit = xx_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
    ang_fit = yy_train;
    
    % emperical
    exp_Pturn(ii) = length(find(abs(ang_fit)>turn_thr)) / (length(ang_fit)-sum(isnan(trials_fit)));
    
    % model
    filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
    K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
    dc_dth = (filt_ddc + filt_dth + b_mle(ii)*0);
    Pturns = (A_mle(ii)-C_mle(ii))./ (1 + exp( -(dc_dth) )) + C_mle(ii);
    Pturns = Pturns*1;
%     p_thr_beta = 2*pi*besseli(0,K2_mle(ii)^1) * exp(K2_mle(ii)^1*cos( ang_fit*d2r-pi ))*(0.1) + (1-0.1)/(2*pi);
    model_Pturn(ii) = nansum(Pturns) / (length(Pturns)-sum(isnan(Pturns))) * ((180-turn_thr)/180*(1-gamm)+gamm); %%%%% estimate the P(>thr|beta=1) density!!!

end

figure
% bar([exp_Pturn; model_Pturn]')
learn_bar_plot([exp_Pturn; model_Pturn], 'turn/s')

%% angular density plots
ii = 1;
gamm = 0.2;
load(datas{ii});
figure
[xx_train, yy_train, mask_train] = data2xy(Data);
%%% turn P
filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
dc_dth = (filt_ddc + filt_dth + b_mle(ii)*0);
Pturns = (A_mle(ii)-C_mle(ii))./ (1 + exp( -(dc_dth) )) + C_mle(ii);
P_beta = nansum(Pturns) / (length(Pturns)-sum(isnan(Pturns)));
% dth histogram
hh = histogram(yy_train, 100, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1)*pi/180;
scal = sum(1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)  + ... 
       ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta );
plot( bb*180/pi, 1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)*1/scal, 'b'); hold on
plot( bb*180/pi, ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta*1/scal,'r');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot function
function learn_bar_plot(ap_na_av, y_name, hh)
    if nargin==2
        figure;
    else
        disp('hold on')
    end
    
    cols = ['b','k','r'];
    for ii = 1:3
        bar(ii, ap_na_av(:,ii), 'FaceColor',cols(ii));
        hold on
    end
    xticks(1:3);
    xticklabels( {'appetitive' 'naive', 'aversive'});
    % set(gca, 'XTickLabel', {'appetitive' 'naive', 'aversive'})
    ylabel(y_name)
    set(gcf,'color','w'); set(gca,'Fontsize',20);
end
