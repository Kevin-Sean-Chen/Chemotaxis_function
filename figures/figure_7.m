%% Figure 7 for mutants
clear
clc
% USER INSTRUCTION: download and unzip the demo data pack Chen_states_2024.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_states_2024');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_states_2024');

%% load data
%%%
% after running inference, with the data set allXs (as,dis,dC,dcp) and the
% trained mmhat structure, we can analyze the distributions and parameters
% correding label of the specific figure panel are shown in the heading of
% each block of code
%%%
% load example tracks (comment out others to plot the given mutant)
load(fullfile(datadir,'data4plot/mutant', 'Data_AIB_nai_staPAW_vars.mat'));
% load(fullfile(datadir,'data4plot/mutant', 'Data_AIY_nai_staPAW_vars.mat'));
% load(fullfile(datadir,'data4plot/mutant', 'Data_AIZ_ave_staPAW_vars.mat'));

%% vectorize data content from structure
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:length(alldis)];
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end

%% emission analysis
%%% load state condition, change "stateK" to inspect
stateK = 2;
x = squeeze(mmhat.wts(:,:,stateK));
[aa,bb] = max( gams_ ,[], 1 );
pos = find(bb==stateK)+offset;
% pos = wind_test;
dcp_K = xxf(2,pos); ddc_K = xxf(1,pos);  ang_K = yyf(1,pos)*1; dis_K = yyf(2,pos);

%%% compute summary statistics
K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
ks_ = x(14:15);  thetas_ = x(16:17);

xv = 1:length(cosBasis);
tt = [0:length(xv)-1]*5/14;
K_h_rec = Amp_h*exp(-xv/tau_h);
K_dc_rec = B_*cosBasis';
K_dcp_rec = Amp*exp(-xv/tau);

%% Figure 3C
%%% transitional kernels
tt = [0:size(cosBasis,1)-1]*5/14;
figure
subi = 1;
for ii = 1:nStates
    for jj = 1:nStates
        if ii ~= jj
%         subplot(nStates, nStates, subi)
        alpha_tran_ = squeeze(mmhat.wts_state(ii,jj,:));
        K_trans_ = alpha_tran_' * cosBasis';
        plot(tt,(K_trans_)); hold on
%         alpha_tran_ = squeeze(mmhat.wts_state(ii,jj,:));
%         K_trans_ = alpha_tran_' * cosBasis';
%         plot(tt,(K_trans_)); 
        end
        subi = subi + 1;
    end
end
set(gcf,'color','w'); set(gca,'Fontsize',20); xlabel('time (s)'); ylabel('weights')

%% Figure 4 panels
%%% plots of densities, fits, conditioned across states
figure()
CM = ['r','k'];%['k','r','g','b'];
dis_cutoff = 10;  % displacement cutoff
fac_mms = 1/30;  % convert to mm displacement
nbins = 100;  %77
for sk = 1:nStates %nStates:-1:1 %
    stateK = sk;
    x = squeeze(mmhat.wts(:,:,stateK));
    %%% armax method
    [aa,bb] = max( gams_ ,[], 1 );
    pos = find(bb==stateK)+offset;
%     %%% probablistic sampling method
%     rand_temp = rand(1,length(gams_));
%     pos = find(gams_(stateK,:)>rand_temp) + offset;
    
    %%% asigning state and parameters
    dcp_K = xxf(2,pos); ddc_K = xxf(1,pos);  ang_K = yyf(1,pos)*1; dis_K = yyf(2,pos);
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
    ks_ = x(14:15);  thetas_ = x(16:17); phis = [0,0];
    K_h_rec = Amp_h*exp(-xv/tau_h);
    K_dc_rec = B_*cosBasis';
    K_dcp_rec = Amp*exp(-xv/tau);
    
    Pk_frac = length(pos)/length(gams_);  % fraction in state K
    subplot(132)
    [counts, edges] = histcounts(ang_K, nbins);
    logCounts = (counts)/sum(counts) * Pk_frac;
    bar((edges(2:end)+edges(1:end-1))/2*pi/180, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); hold on
    bb = (edges(2:end)+edges(1:end-1))/2*pi/180; 
    filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
    filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
    dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
    Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_;
    Pturn_fac = sum(Pturns)/length(Pturns);
    wv_dense = 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb ))*(1-Pturn_fac);
    pr_dense = (1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi))*Pturn_fac;
    scal_fac = sum(wv_dense+pr_dense);
%     rescale_factor = max(logCounts)/max(wv_dense/1);
    plot(bb, ( wv_dense * 1/scal_fac *Pk_frac), CM(sk), 'LineWidth',2); hold on
    plot(bb, ( pr_dense * 1/scal_fac *Pk_frac), '--', 'Color',CM(sk), 'LineWidth',2)
    title('d\theta'); set(gca,'Fontsize',20); set(gcf,'color','w');
    
    subplot(133)
    bar((edges(2:end)+edges(1:end-1))/2*pi/180, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); hold on
    plot(bb, ( wv_dense * 1/scal_fac *Pk_frac), CM(sk), 'LineWidth',2); hold on
    plot(bb, ( pr_dense * 1/scal_fac *Pk_frac), '--', 'Color',CM(sk), 'LineWidth',2)
    ylim([0,0.02]);set(gca,'Fontsize',20); set(gcf,'color','w');
    
    subplot(131)
    dv_k = dis_K(find(dis_K<dis_cutoff));
    Pturns = [1 Pturns];
    p_turn_k = Pturns(find(dis_K<dis_cutoff));
    [counts, edges] = histcounts(dv_k, nbins);
    bb = edges(1:end-1);
    logCounts = (counts)/sum(counts) * Pk_frac;
    bar(edges(1:end-1)*fac_mms, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); hold on
    bb = edges+0.1; 
    gam_p = gampdf(bb, ks_(1), thetas_(1)) * Pturn_fac;
    gam_w = gampdf(bb, ks_(2), thetas_(2)) * (1-Pturn_fac);
    scal_fac = sum(gam_p + gam_w);
    plot(bb*fac_mms, gam_p* 1/scal_fac * Pk_frac, '--', 'Color',CM(sk), 'LineWidth',2)
    plot(bb*fac_mms, gam_w* 1/scal_fac * Pk_frac, 'Color',CM(sk), 'LineWidth',2)
    title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w'); ylabel('probability')
    plot(bb*fac_mms, (gam_p+gam_w)/scal_fac*Pk_frac,'k--')
    
end