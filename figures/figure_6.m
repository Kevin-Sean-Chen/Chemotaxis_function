%% Figure 6 kernels
clear
clc
% USER INSTRUCTION: download and unzip the demo data pack Chen_states_2024.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_states_2024');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_states_2024');

%% load data
load(fullfile(datadir,'data4plot/optogenetics', 'odor_opto_AWC_vars_linear2.mat'));  %%% for aversive worm in odor
% load(fullfile(datadir,'data4plot/optogenetics', 'odor_opto_RIM_vars_linear2.mat'));

%% format struct to vectors
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
allopto = extractfield(Data, 'opto');  % Data
xxf = [xxf; allopto];
yyf = [yyf; alldis];

% testing
wind_test = [1:100000];
offset = min(wind_test)-1;
yy = yyf(:,wind_test);
xx = xxf(:,wind_test);
mask = maskf(wind_test);
% mask = ones(1,length(wind_test));
% mask(find(mask==false)) = NaN;
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data);
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end
% allxys = reshape(extractfield(Data, 'xy'),[],2)';

%% transitional kernels
tt = [0:size(cosBasis,1)-1]*5/14;
figure
subplot(121)
alpha_tran_ = squeeze(mmhat.wts_state(1,2,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(tt,(K_trans_)); hold on
alpha_tran_ = squeeze(mmhat.wts_state(2,1,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(tt,(K_trans_)); 

subplot(122)
alpha_tran_ = squeeze(mmhat.w_state_opto(1,2,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(tt,(K_trans_)); hold on
alpha_tran_ = squeeze(mmhat.w_state_opto(2,1,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(tt,(K_trans_)); 

%% emission analysis
%%% load state condition
stateK = 2;
x = squeeze(mmhat.wts(:,:,stateK));
[aa,bb] = max( gams_ ,[], 1 );
pos = find(bb==stateK);%+offset;
% pos = wind_test;
mask_K = mask(pos);
dcp_K = xx(2,pos).*mask_K; ddc_K = xx(1,pos).*mask_K;  ang_K = yy(1,pos).*mask_K; dis_K = yy(2,pos).*mask_K; 

%%% compute summary statistics
K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(16); b_dcp=x(17);
k_ = x(14);  theta_ = x(15);
Opto_ = x(18:21);

tv = 1:length(cosBasis);
tt = [0:length(tv)-1]*5/14;
K_h_rec = Amp_h*exp(-tv/tau_h);
K_dc_rec = B_*cosBasis';
K_dcp_rec = Amp*exp(-tv/tau);
K_opto_rec = Opto_*cosBasis';

% kernels
figure
subplot(1,4,1); plot(tt,K_dc_rec); xlabel('time (s)'); title('K_C'); set(gca,'Fontsize',20);
subplot(1,4,2); plot(tt,K_h_rec); xlabel('time (s)'); title('K_h'); set(gca,'Fontsize',20);
subplot(1,4,3); plot(tt,K_dcp_rec); xlabel('time (s)'); title('K_{dC^{\perp}}'); set(gca,'Fontsize',20);
subplot(1,4,4); plot(tt,K_opto_rec); xlabel('time (s)'); title('K_{opto}'); set(gca,'Fontsize',20);

% densities
filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_; %+sb
Pturn_fac = sum(Pturns)/length(Pturns);

figure
nbins = 100;
hh = histogram(ang_K*pi/180, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1);
wv_dense = 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb ))*(1-Pturn_fac);
pr_dense = (1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi))*Pturn_fac;
scal_fac = sum(wv_dense+pr_dense);
plot(bb, wv_dense/scal_fac); hold on
plot(bb, pr_dense/scal_fac)
set(gcf,'color','w'); set(gca,'Fontsize',20);

figure
plot(dc_dth/length(K_dcp_rec) , Pturns,'o')
title('P(\beta=1)')

figure
hh = histogram(dis_K, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1);
gam_p = gampdf(bb, k_, theta_);
plot(bb, gam_p/sum(gam_p))
    
