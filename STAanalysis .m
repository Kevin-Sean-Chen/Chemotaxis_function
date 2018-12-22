% STA_analysis 
load('NAI.mat')
%pre-processing
pos = [];
for i = 1:size(alltrigs,1)
    if sum(abs(alltrigs(i,:)))==0
       pos = [pos i];
    end
end

%covariance
X = alltrigs;
X(pos,:) = [];
mu = mean(X);
Xz = X;
for i = 1:size(X,1)
    Xz(i,:) = X(i,:) - mu;
end
covX = cov(Xz);
 
%dimension reduction
cutoff = 50;
[u,s,v] = svd(covX);
basisX = u(:,1:50);
X_est = Xz*basisX;
covX_est = cov(X_est);

%STA
STA_est = mean(X_est*basisX');
plot(STA_est)

%STC
[u,s,v] = svd(covX_est);
plot(s,'-o')

%%%%% stimuli
stim = [];
for i = 1:30%length(Tracks)
    stim = [stim Tracks(i).LEDPower];
end
win = size(alltrigs,2);
stimX = zeros(length(stim)-win,win);
for i = 1:length(stim)-win
    stimX(i,:) = stim(i:i+win-1);
end

%covariance
muS = mean(stimX);
stimXz = X;
for i = 1:size(X,1)
    stimXz(i,:) = stimX(i,:) - muS;
end
covS= cov(stimXz);
 
%dimension reduction
[u,s,v] = svd(covS);
basisS = u(:,1:50);
S_est = Xz*basisS;
covS_est = cov(S_est);

%%%%% MLE results
bML = inv(covS_eat)*
bML = bML*basis