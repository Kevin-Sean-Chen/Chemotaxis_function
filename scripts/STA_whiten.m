%STA_whiten
%%%spike triggered analysis with de-correlation opterations
%%%investigate effects of different nonlinearity and color-noise for STC

%% Orthogonal basis
t = 0:1:30;
ll = length(t);
a = 0.8;
k11 = sin(t/2);  %two basis for covariance analysis
k22 = exp(t/10);
k33 = a*exp(-a*t).*((a*t).^5/factorial(4) - (a*t).^5/factorial(6));%;
k11 = k11./norm(k11);  k11 = k11-mean(k11); k22 = k22./norm(k22); k22 = k22-mean(k22);  k33 = k33./norm(k33); k33 = k33-mean(k33);
                                                                                                                                                                                                                                                                          
ks = [k11; k22; k33];
wall = eye(3);  wall(end) = 1;
Call = ks'* wall * ks;
Q = orth(Call);  %decompose to three orthonormal basis!!

k2 = Q(:,1)'; k1 = Q(:,2)'; k3 = Q(:,3)';
%STA filter
K = k1;
%STC basis
ws = [1, 0; 0 -1];
C = [k2; k3]' * ws * [k2; k3];

%% Kernels
% t = 0:1:30;
% ll = length(t);
% a = 0.8;
% K = a*exp(-a*t).*((a*t).^5/factorial(5) - (a*t).^7/factorial(7));  %biphasic kernel
% K = K-mean(K);  K = K./norm(K);
% % plot(K)
% 
% %%%covariance matrix
% a = 0.8;
% k1 = a*exp(-a*t).*((a*t).^5/factorial(4) - (a*t).^5/factorial(6));%sin(t/4);  %two basis for covariance analysis
% k2 = exp(t/10);%sin(t/4);
% k1 = k1./norm(k1);  k1 = k1-mean(k1); k2 = k2./norm(k2); k2 = k2-mean(k2);
% % plot(k1); hold on; plot(k2);  hold off;
% 
% ws = [-1, 0; 0 1];
% ks = [k1;k2];
% C = ks' * ws * ks;
% %%%construct with orthonormal basis
% Q = orth(C);
% k1 = Q(:,1);  k2 = Q(:,2);  ks = [k1'; k2'];
% C = ks' * ws * ks;
% % % imagesc(C)

%% Stimuli
T = 50000;%20240;
% test = dsp.ColoredNoise('Color','pink');   % pink noise stimuli
% stim = test();
stim = randn(1,T)';   % Gaussian white noise stimuli
tau_cor = 15;
stim = conv(stim,exp(-(1:20)/tau_cor),'same');   % correlated noise with
% correlation time tau_cor

stim = stim-mean(stim); stim = stim./var(stim);

%% Spike simulation
spkr = zeros(1,T);
for i = 1:T-ll
    swin = stim(i:i+ll-1);  %stimulus window
    L = 0.5*swin'*C*swin + (K)*swin;  %linear part
    spkr(i) = L;
    spkr(i) = exp(L*0.2);   % exponential nonlinearity
%     spkr(i) = 1./(1+exp(L*1));   % logistic function
%     spkr(i) = max(0,L);   % ReLU
end
% % temp = conv(fliplr(K),stim)';
% % spkr = temp(ll:end);
% plot(spkr); hold on; plot(stim,'r'); hold off;

%% %%%%%% STA %%%%%%
%% Design matrix
X = zeros(T,ll);
for i = 1:T-ll
    X(i,:) = stim(i:i+ll-1)';
end
%% STA estimation
Kest = (X'*X)\X'*spkr'; %same as %inv(X'*X)*X'*spkr';
Kest = Kest-mean(Kest);  Kest = Kest./norm(Kest);
% plot(Kest); hold on; plot(K,'--r')

%% STC estimation
Ctemp = zeros(T-ll,ll,ll);
for i = 1:T-ll
    swin = stim(i:i+ll-1);
    Ctemp(i,:,:) = spkr(i)* (swin-Kest)*(swin-Kest)';   % stacking STC matrix
end
STC = squeeze(sum(Ctemp,1))/sum(spkr);
% imagesc(STC)

%%% dimention reduction to remove noise %%%
XX = cov(X);
% [u,s,v] = svd(STC);
% cutoff = 2;
% STCn = u(:,1:cutoff) * s(1:cutoff,1:cutoff) * v(:,1:cutoff)';
%%%%%%

% Cest = (X'*X)\STC;  %%%whitening for correlated stimuli
Cest = (XX)\eye(ll) - STC\eye(ll);  %(X'*X)\eye(ll) - STC\eye(ll);   %Baysien STC ML method
% Cest = (X'*X)\eye(ll) * (STC-(X'*X)) * (X'*X)\eye(ll);   %Weiner method

% imagesc(Cest)
[u,s,v] = svd(Cest);
% plot(diag(s),'-o')
k1est = u(:,1);  k1est = k1est./norm(k1est);
k2est = u(:,end);  k2est = k2est./norm(k2est);
% plot(k1,'b');  hold on; plot(k1est,'--b'); plot(k2,'r'); plot(k2est,'--r');  hold off;

%% Adding Prior and minimize MSE
% lamb = 15;
% Lap = 2*eye(ll); Lap(1) = 1; Lap(end) = 1;
% Lap = Lap + diag(-1*ones(1,ll-1),1) +  diag(-1*ones(1,ll-1),-1);
% fun = @(x) sum(sum((Cest-x'*ws*x).^2)) + lamb*( x(1,:)*Lap*x(1,:)' + x(2,:)*Lap*x(2,:)' );  %lamb*sum(sum(x*Lap*x'));  %lamb*x(1,:)*Lap*x(1,:)';
% k0 = [k1est k2est]';%randn(2,ll);
% k = fminsearch(fun,k0);
% k(1,:) = k(1,:)./norm(k(1,:));  k(2,:) = k(2,:)./norm(k(2,:));
% 
% plot(k1,'b');  hold on; plot(k(1,:),'--b'); plot(k2,'r'); plot(k(2,:),'--r');  hold off;

%% %% Adding Prior and minimize MSE---2
% lamb = 100;
Lap = 2*eye(ll); Lap(1) = 1; Lap(end) = 1;
Lap = Lap + diag(-1*ones(1,ll-1),1) +  diag(-1*ones(1,ll-1),-1);
fun = @(x) -(0.5*trace((x'*ws*x) * Cest) + 0.5*K*(x'*ws*x)*K' - length(spkr)/sum(spkr) * det(eye(ll) - XX*(x'*ws*x))^(-0.5) * exp(0.5*K*((XX\eye(ll) ...
    - (x'*ws*x))\eye(ll))*K')  - 0.5*lamb*( x(1,:)*Lap*Lap'*x(1,:)' + x(2,:)*Lap*Lap'*x(2,:)' ) );%- 0.5*lamb*sum(sum(x*Lap*Lap'*x')) );
k0 = [k1est k2est]';%randn(2,ll);
options = optimset('MaxFunEvals',10000);
k = fminsearch(fun,k0,options);
k(1,:) = k(1,:)./norm(k(1,:));  k(2,:) = k(2,:)./norm(k(2,:));

% plot(k2,'b','LineWidth',2.5);  hold on; plot(k(1,:),'--b','LineWidth',5); plot(k3,'r','LineWidth',2.5); plot(k(2,:),'--r','LineWidth',5);  plot(k1est,'.-b'); plot(k2est,'.-r');  hold off;
%% MSE calculation
MSE_K = sum((Kest'-K).^2);
mink1 = min([sum((k1'-k1est).^2), sum((k1'+k1est).^2),sum((k1'-k2est).^2), sum((k1'+k2est).^2)]);   %might be different signs so try all combinations
mink2 = min([sum((k2'-k2est).^2), sum((k2'+k2est).^2),sum((k2'-k1est).^2), sum((k2'+k1est).^2)]);
MSE_C = mink1 + mink2;

%% %%%%% parameter screening %%%%%
%%%for smoothness
% smos = 1:2:30;
% MSE_STA = zeros(1,length(smos));  MSE_STC = zeros(1,length(smos));
% for rep = 1:10
% for itr = 1:length(smos); smn = smos(itr); STA_whiten; MSE_STA(itr) = MSE_K; MSE_STC(itr) = MSE_C; end
% plot(smos,MSE_STC); hold on
% end

%%%for data size
% Ts = logspace(2,5,5);
% MSE_STA = zeros(1,length(Ts));  MSE_STC = zeros(1,length(Ts));
% for rep = 1:10
% for itr = 1:length(Ts); T = round(Ts(itr)); STA_whiten; MSE_STA(itr) = MSE_K; MSE_STC(itr) = MSE_C; end
% semilogx(Ts,MSE_STC); hold on
% end

%%%for nonlinearity
% rr = 30;
% MSE_w1 = zeros(1,rr);
% for rep = 1:10
%     STA_whiten; MSE_w1(rep) = mink1; MES_w2(rep) = mink2;
% end

%%%for regularization
% for kkk = 1:10
% ls = [0 1 10 50 100 250 500 1000];
% MSE_STA = zeros(1,length(ls));  MSE_STC = zeros(1,length(ls));
% for itr = 1:length(ls)
%     lamb = ls(itr);
%     STA_whiten;
%     MSE_STA(itr) = MSE_K; 
%     MSE_STC(itr) = MSE_C;
% end
% plot(ls,MSE_STC); hold on
% end
