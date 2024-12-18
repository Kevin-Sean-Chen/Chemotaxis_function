%% dr_dth_density
%%%
% this code is specialize to produce densities and fits for staPAW
%%%
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_staPAW.mat')

%%
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:250000];  %200000,  40000
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end

%% plot densities, for dr
figure()
CM = ['r','k'];%['k','r','g','b'];
dis_cutoff = 20;  % displacement cutoff
fac_mms = 1/30;  % convert to mm displacement
nbins = 100;  %77
marge_sum = zeros(1,nbins+0);
bbv = linspace(0.05, max(yy(2,:)), nbins);
for sk = nStates:-1:1 %1:nStates %
    stateK = sk;
    x = squeeze(mmhat.wts(:,:,stateK));
    %%% armax method
    [aa,bb] = max( gams_ ,[], 1 );
    pos = find(bb==stateK)+offset;
%     %%% probablistic sampling method
%     rand_temp = rand(1,length(gams_));
%     pos = find(gams_(stateK,:)>rand_temp) + offset;
    
    %%% asigning state and parameters
    dcp_K = xx(2,pos); ddc_K = xx(1,pos);  ang_K = yy(1,pos)*1; dis_K = yy(2,pos);
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
    ks_ = x(14:15);  thetas_ = x(16:17); phis = [0,0]; %phis_ = x(18:19);
    K_h_rec = Amp_h*exp(-xv/tau_h);
    K_dc_rec = B_*cosBasis';
    K_dcp_rec = Amp*exp(-xv/tau);
    
    Pk_frac = length(pos)/length(gams_);  % fraction in state K
    filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
    filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
    dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
    Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_; %+sb
    Pturn_fac = sum(Pturns)/length(Pturns);
    
    subplot(1,3,sk+1)
    dv_k = dis_K(find(dis_K<dis_cutoff));
    Pturns = [1 Pturns];
    p_turn_k = Pturns(find(dis_K<dis_cutoff));
    [counts, edges] = histcounts(dv_k, bbv);%nbins);
    logCounts = (counts)/sum(counts) * 1; %Pk_frac
    nbbv = (bbv(1:end-1)+bbv(2:end))/2;
    bar(nbbv*fac_mms, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5, 'EdgeColor','none'); hold on
    gam_p = gampdf(bbv, ks_(1), thetas_(1)) * Pturn_fac;
    gam_w = gampdf(bbv, ks_(2)*1.0, thetas_(2)) * (1-Pturn_fac);
    scal_fac = sum(gam_p + gam_w);
    plot(bbv*fac_mms, gam_p* 1/scal_fac * 1, '--', 'Color',CM(sk), 'LineWidth',2)  %Pk_frac
    plot(bbv*fac_mms, gam_w* 1/scal_fac * 1, 'Color',CM(sk), 'LineWidth',2)  %Pk_frac
    title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w'); ylabel('probability'); xlim([0 0.4])
    marge_sum = marge_sum + (gam_p+gam_w)/scal_fac*Pk_frac;
    
end

subplot(131);
plot(bbv*fac_mms, marge_sum/sum(marge_sum),'b--'); hold on
dv_f = yy(2, find(yy(2,:)<dis_cutoff));
[counts, edges] = histcounts(dv_f,   bbv);%nbins);
bb = edges(1:end-1);
logCounts = (counts)/sum(counts) * 1;
bar(bb*fac_mms, logCounts,'FaceColor', 'k', 'FaceAlpha',0.5, 'EdgeColor','none');
title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w'); ylabel('probability'); xlim([0 0.4])

%% for dth
figure()
CM = ['r','k'];%['k','r','g','b'];
dis_cutoff = 10;  % displacement cutoff
fac_mms = 1/30;  % convert to mm displacement
nbins = 77;  %77
marge_sum = zeros(1,nbins);
for sk = nStates:-1:1 %1:nStates %
    stateK = sk;
    x = squeeze(mmhat.wts(:,:,stateK));
    %%% armax method
    rand_temp = rand(1,length(gams_));
    pos = find(gams_(stateK,:)>rand_temp) + offset;
    
    %%% asigning state and parameters
    dcp_K = xxf(2,pos); ddc_K = xxf(1,pos);  ang_K = yyf(1,pos)*1; dis_K = yyf(2,pos);
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
    ks_ = x(14:15);  thetas_ = x(16:17); phis = [0,0]; %phis_ = x(18:19);
    K_h_rec = Amp_h*exp(-xv/tau_h);
    K_dc_rec = B_*cosBasis';
    K_dcp_rec = Amp*exp(-xv/tau);
    
    Pk_frac = length(pos)/length(gams_);  % fraction in state K
    subplot(1,3,sk+1)
    [counts, edges] = histcounts(ang_K, nbins);
    logCounts = (counts)/sum(counts) * 1;  %Pk_frac
    bar((edges(2:end)+edges(1:end-1))/2*pi/180, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5, 'EdgeColor','none'); hold on
    bb = (edges(2:end)+edges(1:end-1))/2*pi/180;  %edges*pi/180;
    filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
    filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
    dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
    Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_; %+sb
    Pturn_fac = sum(Pturns)/length(Pturns);
    wv_dense = 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb ))*(1-Pturn_fac);
    pr_dense = (1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi))*Pturn_fac;
    scal_fac = sum(wv_dense+pr_dense);
    plot(bb, ( wv_dense * 1/scal_fac *1), CM(sk), 'LineWidth',2); hold on  %Pk_frac
    plot(bb, ( pr_dense * 1/scal_fac *1), '--', 'Color',CM(sk), 'LineWidth',2)  %Pk_frac
    title('d\theta'); set(gca,'Fontsize',20); set(gcf,'color','w');
    
    marge_sum = marge_sum + (wv_dense+pr_dense) * 1/scal_fac *Pk_frac;
    
end

subplot(131);
plot(bb, marge_sum/sum(marge_sum),'b--','LineWidth',2); hold on
dth_f = yy(1, :);
[counts, edges] = histcounts(dth_f, nbins);
bb = edges(1:end-1);
logCounts = (counts)/sum(counts) * 1;
bar((edges(2:end)+edges(1:end-1))/2*pi/180, logCounts,'FaceColor', 'k', 'FaceAlpha',0.5, 'EdgeColor','none');
title('d\theta'); set(gca,'Fontsize',20); set(gcf,'color','w');
