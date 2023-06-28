% Hess4MLE
% compute Hessian around the MLE fit for each Data conditions

%% load Data and MLE fits
% load data files
datas = {'/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_app_test.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_nai_test.mat',...
         '/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Data_ave_test.mat'};
load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Kfold_mle_param3.mat');

%%
MLE_std = zeros(3,13);  % three conditions and 13 parameters
dx = 0.1;

for ii = 1:3
    load(datas{ii})  % load Data
    mlee = squeeze(mean(mle_params(:,ii,:)))';  % load the fitted MLE as x0
    [H, g] = compHess(@pop_LL, mlee', dx, Data)

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
ttl = {'appetitive','naive','aversive'};
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
K = size(mle_params,1);
figure
for cc = 1:3
    subplot(1,3,cc);
    mlee = squeeze(mean(mle_params(:,cc,:)));
    y = mlee(3:6)'*cosBasis';
    mle_hess = MLE_std(3:6)/sqrt(length(Data));
    standardError = mle_hess*cosBasis';
    plot(tt,y,'k','LineWidth',3)
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
for cc = 1:3
    base = squeeze(mean(mle_params(:,cc,7)));
    bar(xx(cc), base)
    hold on
    berr = squeeze(std(mle_params(:,cc,7)));
    errorbar(xx(cc),base,berr)
end