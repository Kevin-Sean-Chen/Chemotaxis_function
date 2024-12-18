%%% compare kernels
% simple script to load the fitted parameters to compare kernels across conditions
%% load inferred parameters
x = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param4.mat');
x_app = squeeze(mean(x.mle_params(:,1,:),1))';
x_nai = squeeze(mean(x.mle_params(:,2,:),1))';
x_ave = squeeze(mean(x.mle_params(:,3,:),1))';

%% iterations
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
xx = [1:1:length(cosBasis)]*5/14;

figure();
[Kdc,Kdcp] = x2k(x_app,cosBasis)
plot(xx,Kdc); hold on
[Kdc,Kdcp] = x2k(x_nai,cosBasis)
plot(xx(1:end-2),Kdc(3:end)); hold on
[Kdc,Kdcp] = x2k(x_ave,cosBasis)
plot(xx,Kdc); hold on

%%
figure()
[Kdc,Kdcp] = x2k(x_app,cosBasis)
plot(xx,-Kdcp); hold on
[Kdc,Kdcp] = x2k(x_nai,cosBasis)
plot(xx,-Kdcp); hold on
[Kdc,Kdcp] = x2k(x_ave,cosBasis)
plot(xx,-Kdcp); hold on
xlim([0,5])
set(gcf,'color','w'); set(gca,'Fontsize',20);
xlabel('time (s)'); ylabel('weight'); title('K_{dC^{\perp}}')

%% function
function [Kdc,Kdcp] = x2k(x,cosBasis)
    
    K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7); Amp = x(8); tau = x(9); 
    Amp_h = x(10); tau_h = x(11); K2_ = x(12);  gamma = x(13);
    
    xx = 1:length(cosBasis);
    Kdcp = Amp*exp(-xx/tau);
    Kdc = (B_ * cosBasis');
end
