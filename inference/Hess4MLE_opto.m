% Hess4MLE
% compute Hessian around the MLE fit for each Data conditions

%% load Data and MLE fits
% load data files
datas = {'/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_app.mat',...
         '/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_nai.mat',...
         '/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_ave.mat'};
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/mle_param_opto.mat');
% %%% load mle_params_opt here...  %%%%% this is saved in a condition x parameter (3 x 13 x rep) matrix...

%%
MLE_std_opto = zeros(3,13);  % three conditions and 13 parameters
dx = 0.1;

for ii = 1:3
    load(datas{ii})  % load Data
%     mlee = squeeze((mle_params_opto(ii,:)))';  % load the fitted MLE as x0
    mlee = squeeze(nanmedian(mle_params_opto(ii,:,:),3));  % cond x param x repeat
    [H, g] = compHess(@pop_nLL_opto, mlee', dx, Data)

% Compute the inverse of the Hessian matrix
H_inv = inv(H);
% Extract the diagonal elements of the inverse Hessian matrix
variances = diag(H_inv);
% Compute the standard errors by taking the square root of the variances
standard_errors = sqrt(variances);

MLE_std_opto(ii,:) = real(standard_errors)';

end

%%
figure
plot(MLE_std_opto','-o')

%%
%% some post analysis for variability!
% for Kc kernel
ttl = {'appetitive','naive','aversive'};
col = {'b','k','r'};
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
tt = [1:length(cosBasis)]*5/14;
K = size(mle_params_opto,1);
figure
for cc = 1:3
    subplot(1,3,cc);
    mlee = squeeze(nanmedian(mle_params_opto(cc,:,:),3));  %((mle_params_opto(cc,:,5))); %
    y_odor = mlee(3:6)*cosBasis';
    y_opto = mlee(7:10)*cosBasis';
    mle_hess_odor = MLE_std_opto(cc,3:6)*4;%/sqrt(length(Data));   % odor
    mle_hess_opto = MLE_std_opto(cc,7:10)*1;  % opto
    
    standardError_odor = mle_hess_odor*cosBasis';
    yyaxis left; 
    plot(tt,y_odor,col{cc},'LineWidth',3); ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
    ylabel('K_{c}');
    hold on
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y_odor + standardError_odor, fliplr(y_odor - standardError_odor)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    hold off
    
    yyaxis right
    standardError_opto = mle_hess_opto*cosBasis';
    plot(tt,y_opto,'Color',col{cc},'LineWidth',1); yliml = get(gca,'Ylim');% ,'Alpha', 0.3)
    ylabel('K_{opto}');
    hold on
    % Create the shaded area
    xArea = [tt, fliplr(tt)];
    yArea = [y_opto + standardError_opto, fliplr(y_opto - standardError_opto)];
    fill(xArea, yArea, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold off
    
    if yliml(2)*ratio<yliml(1)
        set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
    else
        set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
    end

    xlabel('time (s)');
    set(gca,'FontSize',20); set(gcf,'color','w'); title(ttl{cc})
end

%%
figure
xx = [1,2,3];
param_id = 1;  % the element for plot
for cc = 1:3
    base = squeeze(mean(mle_params_opto(cc,param_id,:)));
    bar(xx(cc), base)
    hold on
%     berr = squeeze(std(mle_params(:,cc,param_id)));
    berr = MLE_std_opto(cc,param_id);%/sqrt(length(Data));
    errorbar(xx(cc),base,berr)
end