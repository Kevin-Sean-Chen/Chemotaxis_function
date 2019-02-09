%Sparse-PCA test
%%%test with sparce PCA on neural dynamics
clc
clear
%% Data
load('Z:\PanNeuronal\20180221\BrainScanner20180221_154553\heatData.mat')

%% Sparse PCA
dF_F = Ratio2(cgIdx,:)';
pos = find(isnan(dF_F));
dF_F(pos) = 0;
%dF_F(:,1:400) = [];  %%%initial removal
dF_F = dF_F - mean(dF_F);
[n,p] = size(dF_F)
t = linspace(0, p, p);

K = 3;
delta = inf;
stop = -[10 5 3];%-[250 125 100];
maxiter = 3000;
convCriterion = 1e-9;
verbose = true;

[B SV L D PATHS] = spca(dF_F, [], K, delta, stop, maxiter, convCriterion, verbose);

%%
figure;
plot(t,sqrt(SV)*ones(1,p).*(B'));  %axis([0 1 -1.2 1.2]);
title('SPCA');
legend('pc1','pc2','pc3','Location','NorthWest');
%ylim([-0.2,0.8])

%%
uu = sqrt(SV)*ones(1,p).*(B');
all1 = dF_F*uu(1,:)';%*B(:,1);
all2 = dF_F*uu(2,:)';%*B(:,2);
all3 = dF_F*uu(3,:)';%*B(:,3);

figure;
subplot(4,3,1:9);scatter3(all1,all2,all3); xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
subplot(4,3,10); plot(all1,all2); xlabel('PC1'); ylabel('PC2');
subplot(4,3,11); plot(all1,all3); xlabel('PC1'); ylabel('PC3');
subplot(4,3,12); plot(all2,all3); xlabel('PC2'); ylabel('PC3');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normal/probablistic PCA
dF_F = Ratio2(cgIdx,:);
pos = find(isnan(dF_F));
dF_F(pos) = 0;
%dF_F(:,1:400) = [];  %%%initial removal
C = cov(dF_F');
[U,S,V] = svd(C,'econ');

figure; plot(diag(S),'-o')

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

