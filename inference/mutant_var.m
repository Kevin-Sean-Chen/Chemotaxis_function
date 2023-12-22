% mutant_var
%%% looping through mutants and conditions to fit the same model;
%%% this script samples CI for std and computes Hessian for uncertainty of
%%% the parameter estimates
%%% This may be important to draw conclusion from finite/noisy data from
%%% mutants...

%% load Data folder
data_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/Data_files';
sav_path = '/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/';
fileList = dir(fullfile(data_path, '*.mat'));

%% for N2
%%% CI computed in the BRW_WV_CI.m script, with weighted CI and std across plates
temp = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/learn_folders4.mat');
n2_ci_std = load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/N2_CI_std.mat');
N2_ci = [n2_ci_std.m_ap(1),    n2_ci_std.m_av(1)   n2_ci_std.m_na(1)];
N2_ci_std = [n2_ci_std.s_ap(1),  n2_ci_std.s_av(1)   n2_ci_std.s_na(1)];
%%% Hess computed from Hess4MLE, loading the result here
% load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/N2_params.mat')
load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/N2_params7.mat') %7,3

%% CI
samps = 20;
CIsamp = zeros(1,3*n_strains, samps);
for li = 1:numel(fileList)
    
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    
    cii = 0;
    c_down = 0;
    ci1 = 0;  ci2 = 0;
    for ss = 1:samps
        data_vec = 1:length(Data);
        shuffledVector = data_vec(randperm(length(data_vec)));
        Data_i = Data(shuffledVector(1:round(length(data_vec)/2)));  % subsample / bootstrap...
    for di = 1:length(Data_i)
        dC = Data_i(di).dc(end) - Data_i(di).dc(1);
        c0 = Data_i(di).dc(1);
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
    CIsamp(1,li,ss) = (cii/ci1 + c_down/ci2);%/ (cii + abs(c_down));
    end
end

%%
CIm = reshape(mean(CIsamp,3),[3,n_strains]);
CIstd = reshape(std(CIsamp,[],3),[3,n_strains]);
figure;
subplot(121); imagesc(CIm)
subplot(122); imagesc(CIstd)

[h,p] = ttest2(n2_ci_std.m_na, squeeze(CIsamp(1,6,:)))

%% sample test for CI
rng(4)
nsamp = 100;
mut_ci_sig = zeros(3,5);  % condition x strains
for cc = 1:3
    for ss = 1:5
        %%% sampling
        n2_samps = N2_ci(cc) + N2_ci_std(cc)*randn(1, nsamp);
        mut_samps = CIm(cc,ss) + CIstd(cc,ss)*randn(1,nsamp);
        %%% tests
        [h, p] = ttest2(n2_samps, mut_samps);
        p
        %%% IMPORTANT: Bonferroni correction for multiple tests
        if p<0.05/(5)
            mut_ci_sig(cc,ss) = 1;
        end
    end
end
figure; imagesc(mut_ci_sig)

%% compute Hessian of parameters
% given this: mle_mut_params ( (learning x strain) x parameters)
mut_param_std = zeros(3*n_strains, n_params);
dx = 0.1;

for li = 1:numel(fileList)
    %%% load Data structure
    fileName = fileList(li).name;
    filePath = fullfile(data_path, fileName);
    load(filePath);
    
    %%% compute Hessian
    mlee = squeeze(mle_mut_params(li,:));  % cond x param x repeat
    [H, g] = compHess(@ll_mut, mlee', dx, Data)
    
    %%% compute std
    % Compute the inverse of the Hessian matrix
    H_inv = inv(H);
    % Extract the diagonal elements of the inverse Hessian matrix
    variances = diag(H_inv);
    % Compute the standard errors by taking the square root of the variances
    standard_errors = sqrt(variances);

    mut_param_std(li,:) = real(standard_errors)';
end

%%
mut_param_std2 = reshape(mut_param_std,[3,n_strains,n_params]);
% figure;
% imagesc(mut_param_std2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to-do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Given Data and Hessian
% Sample posterior of weights N(MLE, Hess) for variability and t-test
% Convolve with Data dc,dcp to get the effect strength


%%
N2_param = squeeze(mle_params(2,:,:));
% N2_param = squeeze(median(mle_params(:,:,:),1));
% N2_param = N2_param([1,3,2],:);
% N2_std = MLE_std([1,3,2],:);

%% make strain x condition mutant filenames
data_mut = cell(3,5);
for ii = 1:15
    fileName = fileList(ii).name;
    data_mut{ii} = fullfile(data_path, fileName);
end

%% sampling posterior of kernels for significance test!
rng(1)
nsamp = 100;  % sample from kernel posterior
mut_sig = zeros(2,3,5);  % 2 x condition x strains
Knorms_n2 = zeros(3, nsamp, 2);  % condition x samples x behavior (brw and wv)
for cc = 1:3
    %%% sampling with normalization
    Data_n2_i = load(data_n2{cc});
    temp_n2_param = N2_param(cc,:);
    mean_brw = BRW_Kc_conv(temp_n2_param, Data_n2_i.Data);
    std_brw = BRW_Kc_conv(N2_std(cc,:), Data_n2_i.Data);
    mean_wv = WV_Kcp_conv(temp_n2_param, Data_n2_i.Data);
    std_wv = WV_Kcp_conv(N2_std(cc,:), Data_n2_i.Data);
    
% for jj = 1:nsamp
    %%% without normalization
%     temp_samp = (N2_param(cc,:)) + (N2_std(cc,:)) .* randn(1,13);
%     Knorms_n2(cc,jj,1) = BRW_id(temp_samp);
%     Knorms_n2(cc,jj,2) = WV_id(temp_samp); 
% end

    Knorms_n2(cc,:,1) = mean_brw + std_brw*randn(1,nsamp);
    Knorms_n2(cc,:,2) = mean_wv + std_wv*randn(1,nsamp);
end
  
for cc = 1:3
    for ss = 1:5

%         load(data_mut{cc,ss});
%         Data_n2_i = load(data_n2{cc});
%         temp_param = squeeze(mut_param(cc,ss,:))';
%         temp_n2_param = N2_param(cc,:);
        %%%%
        %%%% try to sample and convolve here!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%
        Data_mut_i = load(data_mut{cc,ss});       
        
%         %%% samples
%         samp_param = zeros(nsamp, n_params);
%         Knorms = zeros(nsamp, 2);
%         for jj = 1:nsamp
%             temp_samp = squeeze(mut_param(cc,ss,:))' + squeeze(mut_param_std2(cc,ss,:))' .* randn(1,n_params)*1;
%             samp_param(jj,:) = temp_samp;
%             Knorms(jj,1) = BRW_id(temp_samp); %BRW_Kc_conv(temp_param, Data_mut_i.Data, 50000); %
%             Knorms(jj,2) = WV_id(temp_samp); %WV_Kcp_conv(temp_param, Data_mut_i.Data, 50000); %
%         end
        %%% sampling with Gaussian assumption        
        mean_brw = BRW_Kc_conv(squeeze(mut_param(cc,ss,:))', Data_mut_i.Data);
        std_brw = BRW_Kc_conv(squeeze(mut_param_std2(cc,ss,:))', Data_mut_i.Data);
        mean_wv = WV_Kcp_conv(squeeze(mut_param(cc,ss,:))', Data_mut_i.Data);
        std_wv = WV_Kcp_conv(squeeze(mut_param_std2(cc,ss,:))', Data_mut_i.Data);
        Knorms = zeros(nsamp, 2);
        Knorms(:,1) = mean_brw + std_brw*randn(1,nsamp);
        Knorms(:,2) = mean_wv + std_wv*randn(1,nsamp);
        
        %%% tests
        [h, p] = ttest2(Knorms(:,1)', squeeze(Knorms_n2(cc,:,1))');
        p
        %%% IMPORTANT: Bonferroni correction for multiple tests
        if p<0.05/(5*3)
            mut_sig(1,cc,ss) = 1;
        end
        [h, p] = ttest2(Knorms(:,2)', squeeze(Knorms_n2(cc,:,2))');
        if p<0.05/(5*3)
            mut_sig(2,cc,ss) = 1;
        end
    end
end

figure;
subplot(121); imagesc(squeeze(mut_sig(1,:,:))); colorbar(); title('BRW'); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]);yticks([1:3]);%colormap(jet(256))
subplot(122); imagesc(squeeze(mut_sig(2,:,:))); colorbar(); title('WV'); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]);yticks([1:3]);%colormap(jet(256))

%%
WV_Kcp_conv(N2_std(cc,:), Data_n2_i.Data)

%% arrow plots with error bars
figure
% N2_ci = [0.42, 0.02, 0.21];  %%%  compute this with the same way above!!!
arrowHeadSize = 2;
ii = 15;
jj = 1;
testv = zeros(3,5);
for cc = 1:3
    for ss = 1:5
        %%% plot arrow
        temp_param = squeeze(mut_param(cc,ss,:))';
        temp_std = squeeze(mut_param_std2(cc,ss,:))';
        temp_n2_param = N2_param(cc,:);
        temp_n2_ci = N2_ci(cc);
        
        %%% kernel norm setting
%         d_brw = (BRW_id(temp_param) - BRW_id(temp_n2_param));% / BRW_id(temp_n2_param);  % BRW index comparison
%         d_wv =  WV_id(temp_param) -  WV_id(temp_n2_param);
%         d_wv = (-temp_param(8) - -temp_n2_param(8));% / temp_n2_param(8);  % WV index comparison
        d_ci = CIm(cc,ss) - temp_n2_ci;
        
        %%% convolution setting
        load(data_mut{cc,ss});%{cc,ss});
%         fileName = fileList(jj).name;
%         filePath = fullfile(data_path, fileName);
%         load(filePath);
        Data_n2_i = load(data_n2{cc});
        d_brw = BRW_Kc_conv(temp_param, Data) - BRW_Kc_conv(temp_n2_param, Data_n2_i.Data); 
        d_wv =  WV_Kcp_conv(temp_param, Data) - WV_Kcp_conv(temp_n2_param, Data_n2_i.Data); 
%         d_brw = mut_idex2(1,cc,ss);
%         d_wv = mut_idex2(2,cc,ss);
        brw_n2 = BRW_Kc_conv(temp_n2_param, Data_n2_i.Data);
        wv_n2 = WV_Kcp_conv(temp_n2_param, Data_n2_i.Data);
        testv(cc,ss) = d_brw;
        
        subplot(131)
        plot(N2_ci(cc), ii,'o'); hold on
        errorbar(N2_ci(cc), ii, N2_ci_std(cc),'Horizontal',  'Color', 'r');
        h0 = quiver(temp_n2_ci, ii, d_ci, 0,0,  'autoScale', 'off'); hold on
        errorbar(CIm(cc,ss), ii, CIstd(cc,ss),'Horizontal',  'Color', 'k');
        
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('CI'); yticks([])
        subplot(132)
        semilogx(brw_n2, ii,'o'); hold on
        h1 = quiver((brw_n2), ii, (d_brw), 0,0,  'autoScale', 'off'); hold on
%         errorbar((BRW_id(temp_param)), ii, (BRW_std(temp_std)),'Horizontal',  'Color', 'k');
%         errorbar(BRW_id(temp_n2_param), ii, BRW_std(N2_std(cc,:)),'Horizontal',  'Color', 'r');
        errorbar(BRW_Kc_conv(temp_param, Data), ii, BRW_Kc_conv(temp_std, Data),'Horizontal',  'Color', 'k');
        errorbar(BRW_Kc_conv(temp_n2_param, Data_n2_i.Data), ii, BRW_Kc_conv(N2_std(cc,:), Data_n2_i.Data),'Horizontal',  'Color', 'r');
        
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('BRW'); yticks([])

        subplot(133)
%         plot(-temp_n2_param(8), ii,'o'); hold on
        plot(wv_n2, ii,'o'); hold on
        h2 = quiver(wv_n2, ii, d_wv, 0,0,  'autoScale', 'off'); hold on
%         errorbar((WV_id(temp_param)), ii, (WV_std(temp_std)),'Horizontal',  'Color', 'k');
%         errorbar(WV_id(temp_n2_param), ii, WV_std(N2_std(cc,:))/sqrt(7),'Horizontal',  'Color', 'r');
        errorbar(WV_Kcp_conv(temp_param, Data), ii, WV_Kcp_conv(temp_std, Data),'Horizontal',  'Color', 'k');
        errorbar(WV_Kcp_conv(temp_n2_param, Data_n2_i.Data), ii, WV_Kcp_conv(N2_std(cc,:), Data_n2_i.Data),'Horizontal',  'Color', 'r');
        
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('WV'); yticks([])
%         h2.MaxHeadSize = arrowHeadSize/abs(d_wv);
%         set(h, 'UData', arrowLength * d_wv, 'VData', arrowLength * 0);
%         h.AutoScaleFactor = arrow_size*d_wv;
%         h.MaxHeadSize  = 0.5;
        ii = ii-1;
        jj = jj+1;
        
    end    
end

%% functions
function [nll] = ll_mut(theta, Data_train)
    [xx_train, yy_train, mask_train] = data2xy(Data_train);
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
    ang_fit = yy_train;
    dcp_fit = xx_train(2,:);
    ddc_fit = xx_train(1,:);
    trials_fit = mask_train;
    nll = -nLL_kernel_hist5(theta, ang_fit, dcp_fit, ddc_fit, cosBasis, .1, trials_fit);
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

function [wv] = WV_std(x);
    wv = abs(x(8));
end

function [brw] = BRW_std(x)
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    delta = x(2) - x(7);
    beta = norm(x(3:6)*cosBasis');
    brw = (1)*beta*1;
end

function [varargout] = BRW_Kc_conv(x, Data,  nsamp);
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    tt = 0:length(cosBasis)-1;
    Kc = x(3:6)*cosBasis';
    [xx_train, yy_train, mask_train] = data2xy(Data);%(Data(1:400));%(Data(idvec(1:400)));
    ddc_fit = xx_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
%     brw = nanstd(conv_kernel(ddc_fit.*trials_fit, Kc));
    K_h_rec = x(10)*exp(-tt/x(11));
    filt_dth = conv_kernel(abs(yy_train), K_h_rec);
    filt_dc = conv_kernel(ddc_fit.*trials_fit, Kc/1);
    dc_dth = filt_dc+filt_dth*0;
    compute_nl = (x(2)-x(7))./(1+exp(-dc_dth - 0*(filt_dth))) + x(7);
    test_nl = (x(2)-x(7))./(1+exp(-ddc_fit - 0*(filt_dth))) + x(7);
    if nargout > 1
%         varargout{1} = 1/nanstd((ddc_fit.*trials_fit)) * nanstd(conv_kernel(ddc_fit.*trials_fit, Kc));
        varargout{1} = norm(Kc)/nanstd(ddc_fit.*1);%norm(Kc);   %(nanstd(compute_nl));% - nanstd(filt_dc);
        temp = (conv_kernel(ddc_fit.*trials_fit, Kc));
        varargout{2} = temp(~isnan(temp));
%         varargout{2} = dc_dth;
    else
%         varargout{1} = 1/nanstd((ddc_fit.*trials_fit)) * nanstd(conv_kernel(ddc_fit.*trials_fit, Kc));
        if nargin>2
            rand_t = randi(length(ddc_fit)-nsamp);
            varargout{1} = norm(Kc)/nanstd(ddc_fit(rand_t:rand_t+nsamp).*1);%norm(Kc);   %(nanstd(compute_nl));% - nanstd(filt_dc);
        else
            varargout{1} = norm(Kc)/nanstd(ddc_fit.*1);
        end
    end
end

function [varargout] = WV_Kcp_conv(x, Data,  nsamp)
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    tt = 0:length(cosBasis)-1;
    Kcp = x(8)*exp(-tt/x(9));
    [xx_train, yy_train, mask_train] = data2xy(Data);%(Data(1:400));%(Data(idvec(1:400)));
    dcp_fit = xx_train(2,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
%     wv = nanstd(conv_kernel(dcp_fit.*trials_fit, Kcp));
    if nargout > 1
        varargout{1} = nanstd(conv_kernel(dcp_fit.*trials_fit, Kcp));
        temp = (conv_kernel(dcp_fit.*trials_fit, Kcp));
        varargout{2} = temp(~isnan(temp));
    else
        if nargin>2
            rand_t = randi(length(dcp_fit)-nsamp);
            varargout{1} = nanstd(conv_kernel(dcp_fit(rand_t:rand_t+nsamp).*trials_fit(rand_t:rand_t+nsamp), Kcp));
        else
            varargout{1} = nanstd(conv_kernel(dcp_fit.*trials_fit, Kcp));
        end
    end
end

function [h_brw, h_wv] = significant_test(x1, Data1, x2, Data2);
    alpha = 0.001; % Significance level
    %%% for BRW
    [~,samp1_brw] = BRW_Kc_conv(x1, Data1);
    [~,samp2_brw] = BRW_Kc_conv(x2, Data2);
    [h_brw, p] = vartest2(samp1_brw(1:end), samp2_brw(1:end), alpha, 'Both');
%     [h_brw, p] = vartest(samp1_brw, var(samp2_brw));
    %%% for WV
    [~,samp1_wv] = WV_Kcp_conv(x1, Data1);
    [~,samp2_wv] = WV_Kcp_conv(x2, Data2);
    [h_wv, p] = vartest2(samp1_wv(1:end), samp2_wv(1:end), alpha, 'Both');
%     [h_wv, p] = vartest(samp1_wv, var(samp2_wv));
end