%%% STA_analysis 
load('NAI.mat')
%load('NAI_IO.mat')
%% pre-processing
pos = [];
for i = 1:size(alltrigs,1)
    if sum(abs(alltrigs(i,:)))==0
       pos = [pos i];
    end
end

%% covariance
X = alltrigs;
X(pos,:) = [];
mu = mean(X);
Xz = X;
for i = 1:size(X,1)
    Xz(i,:) = X(i,:) - mean(X(i,:));%mu;
end
%Xz = Xz - repmat(mu,size(Xz,1),1);  %%%
covX = cov(Xz);
 
%% dimension reduction
cutoff = 50;
[u,s,v] = svd(covX);
basisX = u(:,1:cutoff);
X_est = Xz*basisX;
covX_est = cov(X_est);

%% STA
STA_est = mean(X_est*basisX');
plot(STA_est)

%% STC
[u,s,v] = svd(covX_est);
plot(diag(s),'-o')

%% %%% stimuli
% stim = [];
% beha = [];
% for i = 1:30%length(Tracks)
%     stim = [stim Tracks(i).LEDPower];
%     beha = [beha Tracks(i).Behaviors];
% end
% beha = beha(8,:);  %focus on reversal
win = size(alltrigs,2);
stimX = zeros(length(stim)-win,win);
for i = 1:length(stim)-win
    stimX(i,:) = stim(i:i+win-1);
end

%% covariance
muS = mean(stimX);
stimXz = X;
for i = 1:size(X,1)
    stimXz(i,:) = stimX(i,:) - mean(stimX(i,:));%muS);
end
covS= cov(stimXz);
 
%% dimension reduction
[u,s,v] = svd(covS);
basisS = u(:,1:cutoff);
S_est = stimXz*basisS;
covS_est = cov(S_est);

%% %%% MLE results

%% STA
bML = inv(covS_est)*mean(covX_est)';%(covS_est\eye(size(covS_est,1)))*mean(covX_est)';%
bML = bML'*basisX';
plot(bML)

%% STC
CML = (covS_est\eye(size(covS_est,1))) * covX_est;%covS_est\covX_est;%
[u,s,v] = svd(CML);
plot(diag(s),'-o')
plot(u(:,1)'*basisX')
hold on
plot(u(:,end)'*basisX')

%% reconstruction
ww = [u(1,:); u(end,:)];
ss = [s(1) 0; s(end) 0];
red_S = stimX*basisS;
resp = zeros(1,length(beha));
for i = 1:length(beha)-win
    temp = red_S(i,:);
    resp(i) = temp*ww'*ss*ww*temp';% + bML
end

plot(resp/max(resp))
hold on
plot(beha)


