%%VMG_analysis3
%%%
% after running inference, with the data set allXs (as,dis,dC,dcp, opto) and the
% trained mmhat structure, we can analyze the distributions and parameters
%%%

% Data stucture imported to construct testing data
% If this is right after 

% testing
wind_test = [1:300000]; %[100000:148982];%500000:length(allas)];%max(wind):length(allas);
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

%% emission analysis
%%% load state condition
stateK = 2;
x = squeeze(mmhat.wts(:,:,stateK));
[aa,bb] = max( gams_ ,[], 1 );
pos = find(bb==stateK)+offset;
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

%% spatial analysis
pix2mm = 1/31.5;
[aa,bb] = max( gams_ ,[], 1 );
CM = ['k','r','w','g'];%jet(stateK);  % See the help for COLORMAP to see other choices.
figure;
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on
for kk = nStates:-1:1 %1:nStates %
    pos = find(bb==kk)+offset;
%     plot(allxys(1,pos)*pix2mm, allxys(2,pos)*pix2mm,'.')%,'color',CM(kk))
    plot(allxys(1,pos)*pix2mm, allxys(2,pos)*pix2mm,'.','color',CM(kk))
    hold on
end
xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');

%% check emission densities
fac_mms = 1/33*14/5;
figure
[aa,bb] = max( gams_ ,[], 1 );
for kk = nStates:-1:1 %1:nStates %
    pos = find(bb==kk)+offset;
    subplot(121)
    histogram(yy(2,wind_test(pos))*fac_mms,100, 'FaceColor', CM(kk)); hold on
    title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w');
    subplot(122)
%     histogram(allas(wind_test(pos)),100, 'FaceColor', CM(kk)); hold on
    [counts, edges] = histcounts(yy(1,wind_test(pos)), 100);
    logCounts = log(counts);
    bar(edges(1:end-1), logCounts,'FaceColor', CM(kk)); hold on
    title('d\theta'); set(gca,'Fontsize',20); set(gcf,'color','w');
end

%% example
wind_ex = offset+(1:1800);  % 800:1200
figure;
subplot(311)
yyaxis left 
plot([1:length(wind_ex)]*5/14, yy(1,wind_ex)); ylabel('d\theta')
yyaxis right
plot([1:length(wind_ex)]*5/14, xx(1,wind_ex)); ylabel('ppm');set(gca,'Fontsize',20);
subplot(312)
plot([1:length(wind_ex)]*5/14, yy(2,wind_ex)*fac_mms, 'k'); ylabel('mm/s'); set(gca,'Fontsize',20);
subplot(313)
% plot([1:length(wind)]*5/14, reshape(smooth(gams(:,wind),10),2,length(wind)))
plot([1:length(wind_ex)]*5/14, gams_(:,wind_ex-offset))
ylim([-0.05,1.05])
xlabel('time (s)'); ylabel('P(Z)')
set(gca,'Fontsize',20); set(gcf,'color','w');

%% time evolution analysis
bins = 11;
[aa,bb_] = max( gams_ ,[], 1 );
cnt_i = zeros(1,bins-1);
cnt_b = cnt_i*1;
epsv = zeros(2,bins);

tvec = linspace(0,max(alltime),bins); zz = alltime(wind_test);
% tvec = linspace(min(alldC), max(alldC), bins); zz = alldC(wind_test); 
figure
for kk = 1:nStates
for bi = 2:bins%+2
    pos = find(zz>tvec(bi-1) & zz<tvec(bi));
    cnt_i(bi-1) = length(find(bb_(pos)==kk));  % for state occupency
%     cnt_i(bi-1) = length(find(abs(yy(pos))>100));  % for simple turning
    cnt_b(bi-1) = length(pos); 
end

EE = ( ((cnt_i-1)./(cnt_b.^2)) + (cnt_i.^2.*(cnt_b-1)./(cnt_b.^4)) ).^0.5 *1/1;
pstatet = cnt_i./cnt_b;
errorbar(tvec(1:end-1)+mean(diff(tvec))/2, pstatet, EE,'-o','LineWidth',2,'color',CM(kk)); hold on
% plot(tvec(1:end-1)+mean(diff(tvec))/2, pstatet, '-o')
ylim([0.,1.])
end

%% compute dwell time statsitics
bins = 500;
rescale_t = 5/14;
[aa,bb] = max( gams_ ,[], 1 );
[dwell_times] = compute_dwell_time(bb, []);

max_dt = max(max([dwell_times{1} , dwell_times{2}]));
bin_edges = linspace(0,1,bins).*max_dt*rescale_t;
figure;
H2 = histogram(dwell_times{1}*rescale_t, bin_edges,'Normalization', 'pdf');hold on
H1 = histogram(dwell_times{2}*rescale_t, bin_edges,'Normalization', 'pdf'); 
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('pdf'); title('dwell time')
