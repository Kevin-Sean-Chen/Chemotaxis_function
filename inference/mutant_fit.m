% mutant_fit
%%% looping through mutants and conditions to fit the same model

%% load Data folder
data_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/Data_files';
sav_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/';
fileList = dir(fullfile(data_path, '*.mat'));

%% parameter settings
rep = 3;  % repetition per fold to make sure that it is the better fit
n_strains = 5;  % number of stains in the folder
L = length(fileList);  % K-fold cross-validation
n_params = 13;  % number of parameters in our model for now
mle_mut_params = zeros(L, n_params); % L x N
trainLL = zeros(L,1);

%% looping
for li = 1:numel(fileList)
    
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    
    %%% training, with small repetition
    x_temp = zeros(rep,n_params);
    fv_temp = zeros(1,rep);
    for ri = 1:rep   % test with repreats to find the max MLE
        [x_MLE, fval] = MLE_mGLM(Data);
        fv_temp(ri) = fval;
        x_temp(ri,:) = x_MLE;
    end
        
    % record MLE
    [minv,pos] = min(fv_temp);
    x_MLE = x_temp(pos,:);
    mle_mut_params(li, :) = x_MLE;
    trainLL(li) = minv;
end

%% CI
CIs = zeros(1,3*n_strains);
for li = 1:numel(fileList)
    
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    
    cii = 0;
    c_down = 0;
    ci1 = 0;  ci2 = 0;
    for di = 1:length(Data)
        dC = Data(di).dc(end) - Data(di).dc(1);
        c0 = Data(di).dc(1);
        if  dC > 0
            cii = cii + dC/(100-c0);
            ci1 = ci1 + 1;
%             cii = cii + 1;
        else
            c_down = c_down + dC/(c0-10);
            ci2 = ci2 + 1;
        end
    end
%     CIs(li) = cii/di;
    CIs(li) = (cii/ci1 + c_down/ci2);%/ (cii + abs(c_down));
end
CIs = reshape(CIs,[3,n_strains]);
figure;
imagesc(CIs)

%% Analysis
load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/N2_params.mat')
mut_idex = zeros(2,3,n_strains);  % id x condition x strain
N2_param = squeeze(mle_params(5,:,:));
N2_param = N2_param([1,3,2],:);
mut_param = reshape(mle_mut_params,[3,n_strains,13]);
xtickLabels = {'AIA','AIB-','AIY-','AIZ-','RIA-'};
ytickLabels = {'app','ave','nai'};
ii = 1;
for ss = 1:n_strains  % strain
    for cc = 1:3  % condition
        temp_param = squeeze(mut_param(cc,ss,:))';
        temp_n2_param = N2_param(cc,:);
        ii = ii+1;
        mut_idex(1,cc,ss) = -(BRW_id(temp_n2_param)*1 - BRW_id(temp_param));% / BRW_id(temp_n2_param);  % BRW index comparison
        mut_idex(2,cc,ss) = (temp_n2_param(8)*1 - temp_param(8));% / temp_n2_param(8);  % WV index comparison
    end
end

figure;
subplot(131); imagesc(CIs); colorbar(); title('CI'); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]); yticks([1:3]);%colormap(jet(256))
subplot(132); imagesc(squeeze(mut_idex(1,:,:))); colorbar(); title('BRW'); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]);yticks([1:3]);%colormap(jet(256))
subplot(133); imagesc(squeeze(mut_idex(2,:,:))); colorbar(); title('WV'); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]);yticks([1:3]);%colormap(jet(256))

%% arrow plots
figure
N2_param = squeeze(mle_params(2,:,:));
N2_param = N2_param([1,3,2],:);
N2_ci = [0.42, 0.02, 0.21];
arrowHeadSize = 2;
ii = 9;
for cc = 1:3
    for ss = 1:5
        %%% plot arrow
        temp_param = squeeze(mut_param(cc,ss,:))';
        temp_n2_param = N2_param(cc,:);
        temp_n2_ci = N2_ci(cc);
        
        d_brw = (BRW_id(temp_param) - BRW_id(temp_n2_param));% / BRW_id(temp_n2_param);  % BRW index comparison
        d_wv =  WV_id(temp_param) -  WV_id(temp_n2_param);
%         d_wv = (-temp_param(8) - -temp_n2_param(8));% / temp_n2_param(8);  % WV index comparison
        d_ci = CIs(cc,ss) - temp_n2_ci;
        
        subplot(131)
        plot(N2_ci(cc), ii,'o'); hold on
        h0 = quiver(temp_n2_ci, ii, d_ci, 0,0,  'autoScale', 'off'); hold on
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('CI'); yticks([])
        subplot(132)
        semilogx(BRW_id(temp_n2_param), ii,'o'); hold on
        h1 = quiver((BRW_id(temp_n2_param)), ii, (d_brw), 0,0,  'autoScale', 'off'); hold on
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('BRW'); yticks([])
%         h1.MaxHeadSize = 10000000;
%         set(h, 'UData', arrowLength * d_brw, 'VData', arrowLength * 0);
%         h.AutoScaleFactor = arrow_size*d_brw;
%         h.MaxHeadSize  = 5;
        subplot(133)
%         plot(-temp_n2_param(8), ii,'o'); hold on
        plot(WV_id(temp_n2_param), ii,'o'); hold on
        h2 = quiver(WV_id(temp_n2_param), ii, d_wv, 0,0,  'autoScale', 'off'); hold on
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('WV'); yticks([])
%         h2.MaxHeadSize = arrowHeadSize/abs(d_wv);
%         set(h, 'UData', arrowLength * d_wv, 'VData', arrowLength * 0);
%         h.AutoScaleFactor = arrow_size*d_wv;
%         h.MaxHeadSize  = 0.5;
        ii = ii-1;
        
    end    
end

%% test parameter cluster
% Sample data (replace this with your actual data)
data = [mle_mut_params(:,:)];% reshape(CIs,12,1)]; %randn(100, 3);  % 100 observations with 3 features each

[coeff, score, ~, ~, explained] = pca(data);
firstTwoPCs = coeff(:, 1:2);
projectedData = data * firstTwoPCs;

figure
for ii = 1:12
    plot(projectedData(ii, 1), projectedData(ii, 2),'o'); hold on
end
% scatter(projectedData(:, 1), projectedData(:, 2));
xlabel('PC 1');
ylabel('PC 2');
title('Data Projected onto First Two Principal Components');
disp(['Variance explained by PC1: ', num2str(explained(1)), '%']);
disp(['Variance explained by PC2: ', num2str(explained(2)), '%']);

%% linear contribution analysis
P = ones(n_strains, n_strains);  % neuron x neuron perturbation matrix
P(1:size(P,1)+1:end) = 0;
B = zeros(2, n_strains);  % behavior x neuron
Ws = zeros(3,2, n_strains);  % condition x behavior x neuron
xtickLabels = {'AIA','AIB','AIY','AIZ','RIA'};
ytickLabels = {'BRW','WV'};
ttls = {'appetitive', 'aversive', 'naive'};

%%% color scheme
numColors = 256;  % Number of colors in the colormap
r = [linspace(1, 1, numColors), linspace(1, 0.5, numColors)];  % Red component
g = [linspace(0, 1, numColors), linspace(0.5, 1, numColors)];  % Green component
b = [linspace(0, 1, numColors), linspace(1, 1, numColors)];  % Blue component
customMap = [r', g', b'];

figure
for cc = 1:3
    for ss = 1:5
        temp_param = squeeze(mut_param(cc,ss,:))';
        B(1,ss) = BRW_id(temp_param);
        B(2,ss) =  WV_id(temp_param);
    end
    B = B./vecnorm(B,1,2);  % normalize behavior readout?
    wtemp = B/P;%B*inv(P);
    Ws(cc,:,:) = wtemp;
    
    subplot(1,3,cc)
    imagesc(wtemp); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]); yticks([1:2]); title(ttls(cc)); caxis([-0.4,.3])
    set(gcf,'color','w'); set(gca,'Fontsize',20);
%     caxis([-max(abs(wtemp(:))), max(abs(wtemp(:)))]);
end

%%
% save('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/mle_mut_params7.mat','mle_mut_params')

%% functions
function [x_MLE, fval] = MLE_mGLM(Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    lfun = @(x)nLL_kernel_hist5(x, ang_fit, dcp_fit, ddc_fit, cosBasis, 1, trials_fit);
    opts = optimset('display','iter');
%     LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, 1e-0*10, 0.1];
%     UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 50, 1];
%     prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -10, 25, 10, 25, 5, 1.];
    LB = [1e-0, 1e-5, ones(1,nB)*-inf, 0.0 -inf, 1e-0, -inf, 1e-1, -inf, -inf];
    UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, inf, inf];
    prs0 = [50, 0.2, randn(1,nB)*10, 0.01, -1, 25, 1, 25, -5, -0.5];
    prs0 = prs0 + prs0.*randn(1,length(UB))*0.0;
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
    x_MLE = x;
end

function [wv] = WV_id(x);
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    tt = 0:length(cosBasis)-1;
    vec = x(8)*exp(-tt/x(9));
    wv = abs(x(8));% norm(vec) *1; %
end

function [brw] = BRW_id(x)
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    delta = x(2) - x(7);
    beta = norm(x(3:6)*cosBasis');
    brw = (1)*beta*1;
end
