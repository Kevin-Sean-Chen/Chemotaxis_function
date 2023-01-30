%%EM_analysis

%% emission analysis
%%% load state condition
stateK = 1;
x = squeeze(mmhat.wts(:,:,stateK));
[aa,bb] = max( gams ,[], 1 );
pos = find(bb==stateK);
dcp_K = alldcp(pos); ddc_K = alldC(pos);  ang_K = allas(pos);

%%% compute summary statistics
K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11);

figure
subplot(2,2,1)
xx = 1:length(cosBasis);
yyaxis left; plot(Amp*exp(-xx/tau)); hold on
yyaxis right; plot(Amp_h*exp(-xx/tau_h))
title('\delta C^{\perp}, \delta \theta kernel')
subplot(2,2,2)
plot(B_ * cosBasis')
title('\delta C kernel')

subplot(2,2,3)
K_dcp_rec = Amp*exp(-xx/tau);
filt_dcp = conv_kernel(dcp_K, K_dcp_rec);%conv(dcp_fit, fliplr(Amp*exp(-xx/tau)), 'same');
[aa,bb] = hist(-(ang_K - filt_dcp)*pi/180 , 500);
bar( bb, 1/(2*pi*besseli(0,K_^2)) * exp(K_^2*cos( bb )) , 100); hold on
bar( bb, 1/(2*pi*besseli(0,K2_^2)) * exp(K2_^2*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi) , 100,'r');
title('von Mises for \delta C^{\perp}')
xlim([-pi, pi])
subplot(2,2,4)
K_dc_rec = B_*cosBasis';
filt_ddc = conv_kernel(ddc_K, K_dc_rec);
xx_h = 1:length(xx)*1;
K_h_rec = Amp_h*exp(-xx_h/tau_h);
filt_dth = conv_kernel(abs(ang_K), K_h_rec);
dc_dth = filt_ddc + 0*filt_dth;
Pturns = 1 ./ (1 + exp( (dc_dth)) + 0) + 0; %+sb
plot(dc_dth/length(K_dcp_rec) , Pturns,'o')
title('Logistic for \delta C')

%% spatial analysis
pix2mm = 1/31.5;
[aa,bb] = max( gams ,[], 1 );
CM = ['k','w','y'];%jet(stateK);  % See the help for COLORMAP to see other choices.
figure;
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on
for kk = 1:nStates %nStates:-1:1 %
    pos = find(bb==kk);
%     plot(allxys(1,pos)*pix2mm, allxys(2,pos)*pix2mm,'.')%,'color',CM(kk))
    plot(allxys(1,pos)*pix2mm, allxys(2,pos)*pix2mm,'.','color',CM(kk))
    hold on
end
xlabel('x (mm)'); ylabel('y (mm)'); h = colorbar();  ylabel(h, 'ppm');
set(gca,'Fontsize',20); set(gcf,'color','w');
% set ( gca, 'xdir', 'reverse' )

%% time series
% wind = 1:5000;
figure;
plot(allas(wind)); hold on
plot(smooth((bb(wind)-1)*100,10))

%% example
wind = 800:1200;
figure;
subplot(211)
yyaxis left 
plot([1:length(wind)]*5/14, allas(wind)); ylabel('d\theta')
yyaxis right
plot([1:length(wind)]*5/14, alldC(wind)); ylabel('ppm');set(gca,'Fontsize',20);
subplot(212)
% plot([1:length(wind)]*5/14, reshape(smooth(gams(:,wind),10),2,length(wind)))
plot([1:length(wind)]*5/14, gams(:,wind))
ylim([-0.05,1.05])
xlabel('time (s)'); ylabel('P(Z)')
set(gca,'Fontsize',20); set(gcf,'color','w');

%% CV ll
figure()
LLs = [-29618 -5959.7  -6287.8  -13303];
ks = 1:4;
b = bar(ks, LLs);
xlabel('# states')
ylabel('test LL')
b(1).BaseValue = min(LLs)*0.9;
set(gca,'Fontsize',20);
set(gcf,'color','w');
%%
figure
app = [-17078, 6395, 6651];
nai = [-23962, -2833, -2714];
ave = [-27120, -1807, 850];

hBar = bar([app;nai;ave]);
ctr = zeros(3,3);
for k1 = 1:3
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');     
    ydt(k1,:) = hBar(k1).YData;                    
end
names = {'app'; 'nai'; 'ave'};
ylabel('test log-likelood')
set(gca,'xticklabel',names,'FontSize',20)
set(gcf,'color','w');
legend([hBar(1), hBar(2),hBar(3)], 'K=1','K=2','K=3')

%%
nai0 = -32575;
app0 = -61319;
ave0 = -36852;
allLL = [app-app0;
         nai-nai0;
         ave-ave0]';
figure;
plot([1,2,3],allLL,'-o')
xlim([0.9,3.1]);xticks([1 2 3])
ylabel('test log-likelood'); xlabel('# states')
set(gca,'FontSize',20);set(gcf,'color','w');
