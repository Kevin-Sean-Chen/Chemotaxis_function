% Figure 4
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4
% extracted from script 'mVM_analysis.m'

%% load all data, with fully trained parameters
load(fullfile(datadir,'data4plots', 'classify_learn3_vars7.mat'));
datas = {fullfile(datadir,'data4plots', 'Data_app_test2.mat'),...
         fullfile(datadir,'data4plots', 'Data_nai_test2.mat'),...
         fullfile(datadir,'data4plots', 'Data_ave_test2.mat')};

param_mle = mle_params;  % directly import
% ids = params_mat.ids;  % load id for training

%% assign mle parameters
time_scale = 5/14;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
t_vec = 1:length(cosBasis)-1;

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

%% histogram and turning curves; with normalization
numBins = 100;
desiredPointsPerBin = 1000;
xx_cond = linspace(-30,30, numBins);
cols = ['b', 'k', 'r'];
figure
for ii = 1:3
    %%% load data
    load(datas{ii});
    [xx_train, yy_train, mask_train] = data2xy(Data(1:end));
    ddc_fit = xx_train(1,:);
    ang_fit = yy_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
    ddc_fit = ddc_fit(2:end);
    ang_fit = ang_fit(1:end-1);
    trials_fit = trials_fit(2:end);
    
    %%% calculation of turn
    filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
    K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
    filt_dth = conv_kernel(abs(ang_fit).*trials_fit, K_h_rec);
    dc_dth_raw = filt_ddc + 1*(filt_dth);
    dc_dth = dc_dth_raw; %

    %%% plots
    subplot(223)
    h = histogram(dc_dth, xx_cond, 'Normalization', 'probability','Visible', 'off');
    bin_edges_input = h.BinEdges;
    probabilities_input = h.Values/sum(h.Values);
    plot(bin_edges_input(1:end-1), probabilities_input, 'Color', cols(ii),'LineWidth',1.5); hold on
    xlim([min(xx_cond), max(xx_cond)])
    xlabel('filtered signal'); ylabel('probability'); set(gcf,'color','w'); set(gca,'Fontsize',20);
    
    subplot(221)
    P_pir_F = ( (A_mle(ii)-C_mle(ii))./(1+exp(-bin_edges_input(1:end-1)))+C_mle(ii) )*time_scale;
    plot(bin_edges_input(1:end-1), P_pir_F,'Color',cols(ii),'LineWidth',1.5); hold on
    xlim([min(xx_cond), max(xx_cond)])
    xlabel('filtered signal'); ylabel('P(\beta=1|C,d\theta)'); set(gcf,'color','w'); set(gca,'Fontsize',20);
    
    subplot(222)
    sig_output = ((A_mle(ii)-C_mle(ii))./(1+exp(-dc_dth_raw*1))+C_mle(ii))*time_scale;
    yy_cond = linspace(min(sig_output),max(sig_output),numBins);
    yy_cond = linspace(0, 0.1,numBins);
    h = histogram(sig_output, yy_cond, 'Normalization', 'probability','Visible', 'off');
    bin_edges_output = h.BinEdges;
    probabilities_output = h.Values/sum(h.Values);
    P_pir = probabilities_output;       %.*probabilities_input; %.*(P_pir_F); %
    P_pir = P_pir/sum(P_pir);
    plot(P_pir, bin_edges_output(1:end-1), 'Color', cols(ii),'LineWidth',1.5); hold on
    ylim([0, 0.1])
    xlabel('probability'); ylabel('P(\beta=1|C,d\theta)'); set(gcf,'color','w'); set(gca,'Fontsize',20);
    
    subplot(224)
    semilogx(P_pir, bin_edges_output(1:end-1), 'Color', cols(ii),'LineWidth',1.5); hold on
    ylim([min(sig_output),max(sig_output)])
    xlabel('probability'); ylabel('P(\beta=1|C,d\theta)'); set(gcf,'color','w'); set(gca,'Fontsize',20);
%     h = histogram(filt_ddc/1 + 1*(filt_dth), numBins, 'Normalization', 'probability','Visible', 'off');
%     bin_edges = h.BinEdges;
%     probabilities = h.Values;
%     plot(bin_edges(1:end-1), probabilities, 'Color', cols(ii),'LineWidth',1.5); hold on
%     xlabel('raw filtered signal'); ylabel('probability'); set(gcf,'color','w'); set(gca,'Fontsize',20);
end

%% emperical vs. model turn probability (histogram method!)
model_Pturn = zeros(1,3);
exp_Pturn = zeros(1,3);
turn_thr = 150;
gamm = 0.20;
nbins = 50;
d2r = pi/180;

figure
for ii = 1:3
    %%% load data
    load(datas{ii});
%     idvec = randperm(length(Data));

    %%% test with the same CV testing data
%     Data = Data(data_ids(ii,:)); % control the number of tracks 
    indices = ids{ii};  % indices for CV (pre-assigned!)
    test_set = (indices==kth);   %%% this matters!!
    train_set = ~test_set;
    Data_train = Data(train_set);
    %%%
    
    [xx_train, yy_train, mask_train] = data2xy(Data_train);%(Data(1:400));%(Data(idvec(1:400)));
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
    P_beta = nansum(Pturns) / (length(Pturns)-sum(isnan(Pturns)));
    [n, edges] = histcounts(yy_train, nbins, 'Normalization', 'probability');
    bb = edges(1:end-1)*pi/180;
    pos = find(abs(bb) > turn_thr*pi/180);
    scal = sum(1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)  + ... 
       ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta );
    p_thr_beta = (1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) * P_beta / scal;
    model_Pturn(ii) = sum(p_thr_beta(pos));
    exp_Pturn(ii) = sum(n(pos));

end

% figure
% bar([exp_Pturn; model_Pturn]')
learn_bar_plot([exp_Pturn; model_Pturn], 'turn/s')

%% angular density plots
nbins = 100;
ii = 2;
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
hh = histogram(yy_train, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1)*pi/180;
scal = sum(1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)  + ... 
       ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta );
plot( bb*180/pi, 1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)*1/scal, 'b'); hold on
plot( bb*180/pi, ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta*1/scal,'r');

%% direct concentration analysis
numBins = 30;
cols = ['b', 'k', 'r'];
figure
for ii = 1:3
    %%% load data
    load(datas{ii});
    dcs = zeros(1,length(Data));
    for jj = 1:length(Data)
        dcs(jj) = Data(jj).dc(end) - Data(jj).dc(1);
    end
    hh = histogram(dcs, numBins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
end

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
