%% staPAW_analysis
%%%
% after running inference, with the data set allXs (as,dis,dC,dcp) and the
% trained mmhat structure, we can analyze the distributions and parameters
%%%
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_staPAW.mat')

%%
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:150000];  %200000,  40000
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

%%
figure;
for ss = 1:nStates
    stateK = ss;
    x = squeeze(mmhat.wts(:,:,stateK));
    [aa,bb] = max( gams_ ,[], 1 );
    pos = find(bb==stateK)+offset;
    dcp_K = xxf(2,pos); ddc_K = xxf(1,pos);  ang_K = yyf(1,pos)*1; dis_K = yyf(2,pos);
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
    ks_ = x(14:15);  thetas_ = x(16:17); %phis_ = x(18:19);
    K_h_rec = Amp_h*exp(-xv/tau_h);
    K_dc_rec = B_*cosBasis';
    K_dcp_rec = Amp*exp(-xv/tau);
    filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
    filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
    dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
    subplot(122); 
    Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)/std(dc_dth)) + 0) + C_;
%     plot(dc_dth(1:100:end)/1/std(dc_dth) , Pturns(1:100:end),'o'); ylabel('P(q=1|state)'); xlabel('projected signal'); set(gca,'Fontsize',20); hold on
    temp_dc_dth = dc_dth(1:50:end);
    temp_pt = Pturns(1:50:end);
    [~, temp_order] = sort(temp_pt);
    plot(temp_dc_dth(temp_order)/1/std(dc_dth) , temp_pt(temp_order),'-o'); ylabel('P(q=1|state)'); xlabel('projected signal'); set(gca,'Fontsize',20); hold on
    subplot(121);
    plot(tt,K_dc_rec); xlabel('time (s)'); ylabel('weights'); set(gca,'Fontsize',20); set(gcf,'color','w'); hold on
end

%% plot densities, fits, across states
figure()
CM = ['r','k'];%['k','r','g','b'];
dis_cutoff = 10;  % displacement cutoff
fac_mms = 1/30;  % convert to mm displacement
nbins = 77;  %77
for sk = 1:nStates %nStates:-1:1 %
    stateK = sk;
    x = squeeze(mmhat.wts(:,:,stateK));
    %%% armax method
%     [aa,bb] = max( gams_ ,[], 1 );
%     pos = find(bb==stateK)+offset;
%     %%% probablistic sampling method
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
    subplot(132)
    [counts, edges] = histcounts(ang_K, nbins);
    logCounts = (counts)/sum(counts) * Pk_frac;
    bar((edges(2:end)+edges(1:end-1))/2*pi/180, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); hold on
    bb = (edges(2:end)+edges(1:end-1))/2*pi/180;  %edges*pi/180;
    filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
    filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
    dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
    Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_; %+sb
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
%     hh = histogram(dis_K, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
    [counts, edges] = histcounts(dv_k, nbins);
    bb = edges(1:end-1);
    logCounts = (counts)/sum(counts) * Pk_frac;
    bar(edges(1:end-1)*fac_mms, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); hold on
    bb = edges+0.1;  %hh.BinEdges(1:end-1);
    gam_p = gampdf(bb, ks_(1), thetas_(1)) * Pturn_fac;
    gam_w = gampdf(bb, ks_(2), thetas_(2)) * (1-Pturn_fac);
    scal_fac = sum(gam_p + gam_w);
    plot(bb*fac_mms, gam_p* 1/scal_fac * Pk_frac, '--', 'Color',CM(sk), 'LineWidth',2)
    plot(bb*fac_mms, gam_w* 1/scal_fac * Pk_frac, 'Color',CM(sk), 'LineWidth',2)
    title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w'); ylabel('probability')
    plot(bb*fac_mms, (gam_p+gam_w)/scal_fac*Pk_frac,'k--')
    

%     bb = (edges(2:end)+edges(1:end-1))/2;
%     logCounts = (counts)/sum(counts) *Pk_frac;
%     bar(edges(1:end-1)*fac_mms, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); hold on
%     binIndices = discretize(dv_k, bb);
%     prob_dr_wv = zeros(1,length(bb));
%     prob_dr_pr = zeros(1,length(bb));
%     for bin = 1:numel(bb)-1
%         pointsInBin = find(binIndices == (bin));
%         pointsInBin = pointsInBin(~isnan(pointsInBin));
%         prob_dr_pr(bin) = sum(gampdf(zeros(1,length(pointsInBin)-1)+bb(bin), ks_(1)+1*phis_(1)*mean(dv_k(pointsInBin(1:end-1))), thetas_(1))  .*mean(p_turn_k(pointsInBin(2:end))));
%         prob_dr_wv(bin) = sum(gampdf(zeros(1,length(pointsInBin)-1)+bb(bin), ks_(2)+1*phis_(2)*mean(dv_k(pointsInBin(1:end-1))), thetas_(2))  .*(1-mean(p_turn_k(pointsInBin(2:end)))));
%         %%%% TEST joint !?.. fix angle part!!
% %         prob_dr_pr(bin) = sum(lognpdf(zeros(1,length(pointsInBin)-1)+bb(bin), ks_(1)+1*phis_(1)*dv_k(pointsInBin(1:end-1)), thetas_(1)).*...
% %                             (1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( zeros(1,length(pointsInBin)-1)+bb(bin)-pi ))*(gamma) + (1-gamma)/(2*pi))*Pturn_fac.*p_turn_k(pointsInBin(2:end)));
% %         prob_dr_wv(bin) = sum(lognpdf(zeros(1,length(pointsInBin)-1)+bb(bin), ks_(2)+1*phis_(2)*dv_k(pointsInBin(1:end-1)), thetas_(2)).*...
% %                             (1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( zeros(1,length(pointsInBin)-1)+bb(bin) ))*(1-Pturn_fac).*(1-p_turn_k(pointsInBin(2:end)))));
%         
%     end
%     norm_prwv = nansum(prob_dr_pr+prob_dr_wv);
%     plot(bb*fac_mms, prob_dr_pr/norm_prwv*Pk_frac,'--')
%     plot(bb*fac_mms, prob_dr_wv/norm_prwv*Pk_frac,'-')
%     plot(bb*fac_mms, (prob_dr_wv+prob_dr_pr)/norm_prwv*Pk_frac,'k--')
    
end

%% proper probability conditions
% binIndices = discretize(dv_k, bb);
% prob_dr_wv = zeros(1,length(bb));
% prob_dr_pr = zeros(1,length(bb));
% for bin = 1:numel(bb)-1
%     pointsInBin = find(binIndices == (bin));
%     pointsInBin = pointsInBin(~isnan(pointsInBin));
%     prob_dr_pr(bin) = mean(gampdf(zeros(1,length(pointsInBin)-1)+bb(bin), gam_shapes_(1)+gam_phis_(1)*dv_k(pointsInBin(1:end-1)), gam_scales_(1)).*p_turn_k(pointsInBin(2:end)));
%     prob_dr_wv(bin) = mean(gampdf(zeros(1,length(pointsInBin)-1)+bb(bin), gam_shapes_(2)+gam_phis_(2)*dv_k(pointsInBin(1:end-1)), gam_scales_(2)).*(1-p_turn_k(pointsInBin(2:end))));
% end
% plot(bb, prob_dr_pr)
% plot(bb, prob_dr_wv)

%% transitional kernels
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

%% spatial analysis
pix2mm = 1/31.5;
[aa,bb] = max( gams_ ,[], 1 );
CM = ['k','r','w','g'];%jet(stateK);  % See the help for COLORMAP to see other choices.
figure;
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on
for kk = 1:nStates %nStates:-1:1 %
    pos = find(bb==kk)+offset;
%     plot(allxys(1,pos)*pix2mm, allxys(2,pos)*pix2mm,'.')%,'color',CM(kk))
    plot(allxys(1,pos)*pix2mm, allxys(2,pos)*pix2mm,'.','color',CM(kk))
    hold on
end
xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');

%% single track
figure;
imagesc(M); hold on
CM = ['k','r','g','b'];
exp_wind = (2200:3150)+5000;%[1800:2800]+0; %25000:31000;
[aa,cc] = max( gams_(:,exp_wind) ,[], 1 );
for kk = 1:nStates %nStates:-1:1 %
    pos = find(cc==kk)+min(exp_wind)+offset;
    plot(allxys(1,pos), allxys(2,pos),'.','color',CM(kk),'MarkerSize',15)
    hold on
end

%% check emission densities
fac_mms = 1/33*14/5;
figure
[aa,bb] = max( gams_ ,[], 1 );
for kk = 1:nStates %nStates:-1:1 %
    pos = find(bb==kk)+offset*0;
    subplot(121)
    histogram(yyf(2,wind_test(pos))*fac_mms,100, 'FaceColor', CM(kk)); hold on
    title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w');
    subplot(122)
%     histogram(allas(wind_test(pos)),100, 'FaceColor', CM(kk)); hold on
    [counts, edges] = histcounts(yyf(1,wind_test(pos)), 100);
    logCounts = log(counts);
    bar(edges(1:end-1), logCounts,'FaceColor', CM(kk)); hold on
    title('d\theta'); set(gca,'Fontsize',20); set(gcf,'color','w');
end

%% example
wind_ex = offset+(2200:3150)+0;%2800;  %1200:2200   %800:1200
figure;
subplot(311)
yyaxis left 
plot([1:length(wind_ex)]*5/14, yyf(1,wind_ex)); ylabel('d\theta')
yyaxis right
plot([1:length(wind_ex)]*5/14, xxf(1,wind_ex)); ylabel('ppm');set(gca,'Fontsize',20);
% subplot(412)
% plot([1:length(wind_ex)]*5/14, xxf(2,wind_ex)*30, 'k'); ylabel('ppm/mm'); set(gca,'Fontsize',20);
subplot(312)
plot([1:length(wind_ex)]*5/14, yyf(2,wind_ex)*fac_mms, 'k'); ylabel('mm/s'); set(gca,'Fontsize',20);
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

% tvec = linspace(0,max(alltime),bins); zz = alltime(wind_test);
tvec = linspace(0, 20*60, bins); zz = alltime(wind_test); 
figure
for kk = 1:nStates
for bi = 2:bins%+2
    pos = find(zz>tvec(bi-1) & zz<tvec(bi));
    cnt_i(bi-1) = length(find(bb_(pos)==kk));  % for state occupency
%     cnt_i(bi-1) = length(find(abs(yy(pos))>100));  % for simple turning
    cnt_b(bi-1) = length(pos); 
end

EE = ( ((cnt_i-1)./(cnt_b.^2)) + (cnt_i.^2.*(cnt_b-1)./(cnt_b.^4)) ).^.5;
pstatet = cnt_i./cnt_b;
errorbar(tvec(1:end-1)+mean(diff(tvec))/2, pstatet, EE,'-o','LineWidth',2,'color',CM(kk)); hold on
% plot(tvec(1:end-1)+mean(diff(tvec))/2, pstatet, '-o')
ylim([0.,1.])
end

%% compute dwell time statsitics
bins = 100;
rescale_t = 5/14;
[aa,bb] = max( gams_ ,[], 1 );
[dwell_times] = compute_dwell_time(bb, []);

max_dt = max(max([dwell_times{1} , dwell_times{2}]));
bin_edges = linspace(0,1,bins).*max_dt*rescale_t;
figure;
H2 = histogram(dwell_times{1}*rescale_t, bin_edges,'Normalization', 'pdf'); H2.Visible = 'off';%hold on
counts2 = H2.Values;
bin_edges2 = H2.BinEdges(1:end-1) + diff(H2.BinEdges)/2;
H1 = histogram(dwell_times{2}*rescale_t, bin_edges,'Normalization', 'pdf'); H1.Visible = 'off';
counts1 = H1.Values;
bin_edges1 = H1.BinEdges(1:end-1) + diff(H1.BinEdges)/2;
semilogy(bin_edges1, counts1); hold on
semilogy(bin_edges2, counts2)
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('pdf'); title('dwell time'); xlim([0,50])

%% further analysis for autocorrelation
% given state posterior, simulate state-dependent AR process to check autocorrelations for dr
[aa,bb] = max( gams_ ,[], 1 );
sim_dr = zeros(1,length(bb));
sim_dth = sim_dr*1;
drt = yyf(2,1);
xv = 1:length(cosBasis); wind = length(xv);
for tt = wind+1:length(bb)
    %%% choose state!
%     state_t = bb(tt);
    prob = gams_(:,tt);
    if rand()<prob(1)
        state_t = 1;
    else
        state_t = 2;
    end
    
%     state_t
    %%% unpack parameters given state
    x = squeeze(mmhat.wts(:,:,state_t));
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
    ks_ = x(14:15);  thetas_ = x(16:17); phis_ = [0,0]; %phis_ = x(18:19);
    
    %%% for dr
    if state_t==1
        drt = gamrnd(ks_(2) + 1*drt*phis_(2)/1, thetas_(1)/1);
    else
        if rand(1) < mean(Pturns)
            drt = gamrnd(ks_(2) + 1*drt*phis_(2)/1, thetas_(1)/1);
        else
            drt = gamrnd(ks_(1) + 1*drt*(phis_(1))/1, thetas_(2)/1);
        end
    end
    
    %%% for dth
    K_dcp = Amp*exp(-xv/tau);
    K_h = Amp_h*exp(-xv/tau_h);
    K_dc = B_ * cosBasis';
    dCv = xx(1,tt-wind+1:tt); dthv = yy(1,tt-wind+1:tt); dCpv = xx(2,tt-wind+1:tt);
    wv = (1*sum(K_dcp.*dCpv) + b_dcp*1) + (vmrand(0,K_))*180/pi;
    P_event = (A_-C_) / (1+exp( -(sum(K_dc.*dCv)*1. + 1*(sum(K_h.*abs(dthv)*180/pi)) *1 + 1*b_dc)+0) ) + C_;%Pturn_base;%length(wind)
    if rand < P_event
        beta = 1;
%         vv = gamrnd(k_(1), theta_(1))/1;
    else
        beta = 0;
%         vv = gamrnd(k_(2), theta_(2))/1;
    end
    
    if rand < gamma
        rt = beta*(vmrand(pi,1*K2_)*180/pi);
    else
        rt = beta*(vmrand(0,0)*180/pi);
    end
    
    dth = wv+rt; dth = wrapToPi(dth*pi/180)*180/pi;
    sim_dr(tt) = drt;
    sim_dth(tt) = dth;
end

%% analyze <dr,dr'>
model = @(params, x) params(1)*exp(-1/params(2)*x) + params(3)*exp(-1/params(4)*x);% + params(5);

nlag = 60;
dr_data = yyf(2,1:150000);
dv = dr_data(find(dr_data<dis_cutoff));
[As_data,lgs] = autocorr(dv, nlag);
[As_sim,lgs] = autocorr(sim_dr, nlag);

%%% data fitting
initialGuess = [.5, 1, .1, 10];
paramsFit = lsqcurvefit(model, initialGuess, lgs*5/14, (As_data));
yFit = model(paramsFit, lgs*5/14);
figure; plot(lgs*5/14, yFit, 'LineWidth',3); hold on;
plot(lgs*5/14, As_data, 'k.', 'MarkerSize',15);
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('<drdr''>')

%%
%%% model prediction
figure;
plot(lgs*5/14, As_data, '--'); hold on
plot(lgs(1:1:end)*5/14, As_sim(1:1:end),'k')
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('<drdr''>')

%% further analysis of ITI!
pos = find(abs(yy(1,:))>50); %.7
% pos = find(abs(sim_dth)>50); %.6
iti = diff(pos)*5/14;  % duration in seconds
% figure; hist(iti,100)
bins = 0:.7:70;
H1 = histogram(iti, bins, 'Normalization', 'pdf'); H1.Visible = 'off'; %,
iti_cnt = H1.Values;
iti_bin = H1.BinEdges(1:end-1) + diff(H1.BinEdges)/2;

% Double exponential model function
initialGuess = [.5, 1, .1, 100];
paramsFit = lsqcurvefit(model, initialGuess, iti_bin, (iti_cnt));
yFit = model(paramsFit, iti_bin);
figure; 
semilogy(iti_bin, yFit, 'LineWidth',3); hold on;
semilogy(iti_bin, iti_cnt, 'k.','MarkerSize',15); xlim([0, 70]);
a1=paramsFit(1); tau1=paramsFit(2); a2=paramsFit(3); tau2=paramsFit(4);
tau_c = 1/(1/tau1-1/tau2)*log(a1/a2)
set(gcf,'color','w'); set(gca,'Fontsize',20);  ylabel('probability');  xlabel('turn interval (s)')

%%
%%% model prediction
pos = find(abs(sim_dth)>50); %.6
iti_model = diff(pos)*5/14;  % duration in seconds
bins_m = 0:.7:70;
H1 = histogram(iti_model, bins_m, 'Normalization', 'pdf'); H1.Visible = 'off'; %,
iti_cnt_m = H1.Values;
iti_bin_m = H1.BinEdges(1:end-1) + diff(H1.BinEdges)/2;
figure; 
plot(iti_bin_m, iti_cnt_m, 'LineWidth',3); hold on;
plot(iti_bin, iti_cnt, 'k-','MarkerSize',15); xlim([0, 70]);
set(gcf,'color','w'); set(gca,'Fontsize',20);  ylabel('probability');  xlabel('turn interval (s)')

%% descriptive stats plots
figure;
full_data = yy(1,1:150000);
% temp_y = full_data(find(full_data<dis_cutoff));
temp_y = full_data;
[counts, edges] = histcounts(temp_y, 100); %nbins 60, 100
bb = edges(1:end-1);
logCounts = (counts)/sum(counts) * Pk_frac;
% bar(edges(1:end-1)*fac_mms, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); 
bar(edges(1:end-1)*1, logCounts,'FaceColor', CM(sk), 'FaceAlpha',0.5); 
title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w'); ylabel('probability'); %xlim([0, 0.3])

%% for demo
wind_ex = offset+(2350:2800)+0;%2800;  %1200:2200   %800:1200
figure;
subplot(312) 
plot([1:length(wind_ex)]*5/14, yyf(1,wind_ex)); ylabel('d\theta');set(gca,'Fontsize',20); xlim([0,160]); hold on; xticks([])
pos = find(abs(yyf(1,wind_ex))>50);
plot(pos*5/14, ones(1, length(pos)),'*'); plot([0,160],[50,50],'k--'); plot([0,160],[-50,-50],'k--')
subplot(311)
plot([1:length(wind_ex)]*5/14, xxf(1,wind_ex)); ylabel('ppm');set(gca,'Fontsize',20); xlim([0,160]); xticks([])
subplot(313)
plot([1:length(wind_ex)]*5/14, yyf(2,wind_ex)*fac_mms, 'k'); ylabel('mm/s'); set(gca,'Fontsize',20);set(gcf,'color','w'); xlim([0,160])
