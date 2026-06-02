%% Figure panels for fig.1,3,4
clear
clc
% USER INSTRUCTION: download and unzip the demo data pack Chen_states_2024.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_states_2024');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_states_2024');

%% load data
%Data_salt100_50_staPAW2.mat;  Data_nai.mat;  Data_ave.mat
load(fullfile(datadir,'data4plot/salt_chemotaxis', 'Data_salt100_50_staPAW2.mat'));  %%% for 50-100 salt
load(fullfile(datadir,'data4plot/odor_chemotaxis', 'Data_nai_staPAW.mat'));  %%% for naive worm in odor
load(fullfile(datadir,'data4plot/odor_chemotaxis', 'Data_ave_staPAW.mat'));  %%% for aversive worm in odor

%%
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:100000];  %200000,  40000, 150000
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end

%% transitional kernels, shown in the panels
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

%% emission analysis
%%% load state condition
stateK = 2;
x = squeeze(mmhat.wts(:,:,stateK));
[aa,bb] = max( gams_ ,[], 1 );
pos = find(bb==stateK)+offset;
% pos = wind_test;
dcp_K = xxf(2,pos); ddc_K = xxf(1,pos);  ang_K = yyf(1,pos)*1; dis_K = yyf(2,pos);

%%% compute summary statistics
K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
ks_ = x(14:15);  thetas_ = x(16:17); %phis_ = x(18:19);

xv = 1:length(cosBasis);
tt = [0:length(xv)-1]*5/14;
K_h_rec = Amp_h*exp(-xv/tau_h);
K_dc_rec = B_*cosBasis';
K_dcp_rec = Amp*exp(-xv/tau);

% kernels
figure
subplot(1,3,1); plot(tt,K_dc_rec); xlabel('time (s)'); title('K_C'); set(gca,'Fontsize',20);
subplot(1,3,2); plot(tt,K_h_rec); xlabel('time (s)'); title('K_h'); set(gca,'Fontsize',20);
subplot(1,3,3); plot(tt,K_dcp_rec); xlabel('time (s)'); title('K_{dC^{\perp}}'); set(gca,'Fontsize',20);

% densities
filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)/std(dc_dth)) + 0) + C_; %+sb
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
bb = hh.BinEdges(1:end-1)+0.1;
gam_p = gampdf(bb, ks_(1), thetas_(1));
gam_w = gampdf(bb, ks_(2), thetas_(2));
plot(bb, gam_p/sum(gam_p))
plot(bb, gam_w/sum(gam_w))

