% Hess4MLE
% compute Hessian around the MLE fit for each Data conditions

%% load Data and MLE fits
% load data files
datas = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test2.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test2.mat'};
load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param4.mat');  %4
% load('/projects/LEIFER/Kevin/Data_learn/Mutants/data_analysis/N2_params.mat')

%%
MLE_std = zeros(3,13);  % three conditions and 13 parameters
dx = 0.1;

for ii = 1:3
    load(datas{ii})  % load Data
    mlee = squeeze(median(mle_params(:,ii,:)))';  % load the fitted MLE as x0
    [H, g] = compHess(@pop_LL, mlee', dx, Data(1:400))

% Compute the inverse of the Hessian matrix
H_inv = inv(H);
% Extract the diagonal elements of the inverse Hessian matrix
variances = diag(H_inv);
% Compute the standard errors by taking the square root of the variances
standard_errors = sqrt(variances);

MLE_std(ii,:) = real(standard_errors)';

end

%%
figure
plot(MLE_std','-o')

%%
%% some post analysis for variability!
% for Kc kernel
ttl = {'appetitive','naive','aversive'};
col = {'b','k','r'};
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
K = size(mle_params,1);
figure
for cc = 1:3
%     subplot(1,3,cc);
    mlee = squeeze(median(mle_params(:,cc,:)));
    y = mlee(3:6)'*cosBasis';
    mle_hess = MLE_std(cc, 3:6)/sqrt(1);%/sqrt(length(Data));
    standardError = mle_hess*cosBasis';
    plot(tt,y,col{cc},'LineWidth',3)
    hold on
    
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y + standardError, fliplr(y - standardError)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    title(ttl{cc})
end

%% 
% for Kc_perp kernel
ttl = {'appetitive','naive','aversive'};
col = {'b','k','r'};
tt = [1:length(cosBasis)]*5/14;
figure
for cc = 1:3
%     subplot(1,3,cc);
    mlee = squeeze(median(mle_params(:,cc,:)));
    y = -mlee(8).*exp(-tt./mlee(9));
    mle_hess = MLE_std(cc,8:9)/sqrt(8);%/sqrt(length(Data));
    standardError = -mle_hess(1).*exp(-tt./mle_hess(2));
    plot(tt,y,col{cc},'LineWidth',3)
    hold on
    
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y + standardError, fliplr(y - standardError)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    title(ttl{cc})
end

%%
figure
xx = [1,2,3];
param_id = 1;  % the element for plot
for cc = 1:3
    base = squeeze(mean(mle_params(:,cc,param_id)));
    bar(xx(cc), base)
    hold on
%     berr = squeeze(std(mle_params(:,cc,param_id)));
    berr = MLE_std(cc,param_id);%/sqrt(length(Data));
    errorbar(xx(cc),base,berr)
end