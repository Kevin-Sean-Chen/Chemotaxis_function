%C elegans brain dynacmis _practice
clear
%load('Z:\PanNeuronal\20180221\BrainScanner20180221_150819\heatData.mat')
load('Z:\PanNeuronal\20180221\BrainScanner20180221_154553\heatData.mat')

%% PCA
dF_F = Ratio2(cgIdx,:);
pos = find(isnan(dF_F));
dF_F(pos) = 0;
%dF_F(:,1:400) = [];  %%%initial removal
C = cov(dF_F');
[U,S,V] = svd(C,'econ');

figure; plot(diag(S),'-o')

%% P-PCA
[coeff,score,pcvar] = ppca(dF_F',3);

%% dimension reduction
dM = zeros(size(dF_F));
for ii = 1:size(dF_F,1)
    dM(ii,:) = dF_F(ii,:) - mean(dF_F(ii,:));
end

fstPC = U(:,1);
sndPC = U(:,2);
trdPC = U(:,3);

all1 = dM'*fstPC;
all2 = dM'*sndPC;
all3 = dM'*trdPC;

figure;
subplot(4,3,1:9);scatter3(all1,all2,all3); xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
subplot(4,3,10); plot(all1,all2); xlabel('PC1'); ylabel('PC2');
subplot(4,3,11); plot(all1,all3); xlabel('PC1'); ylabel('PC3');
subplot(4,3,12); plot(all2,all3); xlabel('PC2'); ylabel('PC3');

%% state-space
figure;
p1t = zeros(1,size(dF_F,1));
p2t = zeros(1,size(dF_F,1));
p3t = zeros(1,size(dF_F,1));
for tt = 1:size(dF_F,2)-1000
    p1t(tt) = dM(:,tt)' * fstPC;
    p2t(tt) = dM(:,tt)' * sndPC;
    p3t(tt) = dM(:,tt)' * trdPC;
    
    %%%small animations
    %plot3(p1t,p2t,p3t); pause()
    %plot(p3t,p2t); pause()
end

s = repmat([240],numel(p1t),1);
c = linspace(1,tt,tt);
h = scatter3(p1t,p2t,p3t,s,c,'filled')

% plot3(p1t,p2t,p3t)

