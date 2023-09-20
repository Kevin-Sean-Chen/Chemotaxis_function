%%EM_analysis

% testing

wind_test = [1:150000]; %[100000:148982];%500000:length(allas)];%max(wind):length(allas);
offset = min(wind_test)-1;
yy = [allas(wind_test)*1;
      alldis(wind_test)];
xx = [alldC(wind_test); 
      alldcp(wind_test)];
mask = true(1,length(yy));
mask(isnan(alltrials(wind_test))) = false;%
[logp_test,gams_,xisum_test] = runFB_GLMHMM(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data);

%% emission analysis
%%% load state condition
stateK = 2;
x = squeeze(mmhat.wts(:,:,stateK));
[aa,bb] = max( gams_ ,[], 1 );
pos = find(bb==stateK)+offset;
% pos = wind_test;
dcp_K = alldcp(pos); ddc_K = alldC(pos);  ang_K = allas(pos)*1;

%%% compute summary statistics
K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(16); b_dcp=x(17);

figure
subplot(2,2,1)
xx = 1:length(cosBasis);
yyaxis left; plot(Amp*exp(-xx/tau)); hold on
yyaxis right; plot(Amp_h*exp(-xx/tau_h))
title('\delta C^{\perp}, \delta \theta kernel');set(gca,'Fontsize',20);
subplot(2,2,2)
plot(B_ * cosBasis')
title('\delta C kernel');set(gca,'Fontsize',20);

subplot(2,2,3)
K_dcp_rec = Amp*exp(-xx/tau);
filt_dcp = conv_kernel(dcp_K, K_dcp_rec);%conv(dcp_fit, fliplr(Amp*exp(-xx/tau)), 'same');
[aa,bb] = hist(-(ang_K - filt_dcp - b_dcp)*pi/180 , 500);
bar( bb, 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb )) , 5); hold on
bar( bb, 1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi) , 5,'r');
title('von Mises mixture')%for \delta C^{\perp}')
xlim([-pi, pi]);set(gca,'Fontsize',20);
subplot(2,2,4)
K_dc_rec = B_*cosBasis';
filt_ddc = conv_kernel(ddc_K(2:end), K_dc_rec);
xx_h = 1:length(xx)*1;
K_h_rec = Amp_h*exp(-xx_h/tau_h);
filt_dth = conv_kernel(abs(ang_K(1:end-1)), K_h_rec);
dc_dth = filt_ddc*1 + 1*filt_dth + b_dc;
Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_; %+sb
plot(dc_dth/length(K_dcp_rec) , Pturns,'o')
title('P(\beta=1)')%('Logistic for \delta C')
set(gcf,'color','w'); set(gca,'Fontsize',20);

%% transitional kernels
tt = [0:size(cosBasis,1)-1]*5/14;
figure
alpha_tran_ = squeeze(mmhat.wts_state(1,2,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(tt,(K_trans_)); hold on
alpha_tran_ = squeeze(mmhat.wts_state(2,1,:));
K_trans_ = alpha_tran_' * cosBasis';
plot(tt,(K_trans_)); 

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
% set ( gca, 'xdir', 'reverse' )

%% check emission densities
fac_mms = 1/33*14/5;
figure
[aa,bb] = max( gams_ ,[], 1 );
for kk = nStates:-1:1 %1:nStates %
    pos = find(bb==kk)+offset;
    subplot(121)
    histogram(alldis(wind_test(pos))*fac_mms,100, 'FaceColor', CM(kk)); hold on
    title('dr'); set(gca,'Fontsize',20); set(gcf,'color','w');
    subplot(122)
%     histogram(allas(wind_test(pos)),100, 'FaceColor', CM(kk)); hold on
    [counts, edges] = histcounts(allas(wind_test(pos)), 100);
    logCounts = log(counts);
    bar(edges(1:end-1), logCounts,'FaceColor', CM(kk)); hold on
    title('d\theta'); set(gca,'Fontsize',20); set(gcf,'color','w');
end

%% time series
% wind = 1:5000;
% figure;
% plot(allas(wind)); hold on
% plot(smooth((bb(wind)-1)*100,10))

%% example
wind_ex = offset+(1200:1800);  % 800:1200
figure;
subplot(311)
yyaxis left 
plot([1:length(wind_ex)]*5/14, allas(wind_ex)); ylabel('d\theta')
yyaxis right
plot([1:length(wind_ex)]*5/14, alldC(wind_ex)); ylabel('ppm');set(gca,'Fontsize',20);
subplot(312)
plot([1:length(wind_ex)]*5/14, alldis(wind_ex)*fac_mms, 'k'); ylabel('mm/s'); set(gca,'Fontsize',20);
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

%% CV ll
