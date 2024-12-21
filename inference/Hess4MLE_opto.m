% Hess4MLE_opto
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
    mle_hess_odor = MLE_std_opto(cc,3:6)*1;%4/sqrt(length(Data));   % odor
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

%% compare emperical and model-predictive turns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/mle_param_opto2.mat');
%% bar analysis
fr = 1/14;  % frame rate
bin = 7;  % filtering bins
turn_thre = 50;  % turning threshold for emperical analysis
acst = 4 * (1/(bin*fr)); % pre-stim
windt = 14 * (1/(bin*fr));  % post-stim
wind_analysis = 4*2:4*2+5*2;  % time window with impulse
wind_baseline = 1:4*2;  % baseline comparison
t_vec = [-acst:windt]*((bin*fr));
Impulse = zeros(1,length(t_vec));
Impulse(5:5+5*2) = 4.2780;  % opto strength
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/app_trigs.mat')
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/nai_trigs.mat')
load('/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/ave_trigs.mat')
all_ang_trigs = {[app_up_trigs; app_down_trigs],... 
                 [nai_up_trigs; nai_down_trigs],...
                 [ave_up_trigs; ave_down_trigs]};
data_pr_mean = zeros(1,3);
data_pr_std = zeros(1,3);
model_pr_mean = zeros(1,3);
model_pr_std = zeros(1,3);

for cc = 1:3
    %%% data
    load(datas{cc})
    opto_data = extractfield(Data, 'opto');
    dth_data = extractfield(Data, 'dth');
    stim_pos = find(opto_data>4);
    n_turns = zeros(1,length(stim_pos)); 
    n_pre = n_turns*1;
    beh_opt = dth_data(stim_pos);
    n_turns(abs(beh_opt)>turn_thre) = 1;
    beh_pre = dth_data(stim_pos-9);
    n_pre(abs(beh_pre)>turn_thre) = 1;
    data_pr_mean(cc) = (mean(n_turns) - mean(n_pre))*14/5;
    data_pr_std(cc) = (std(n_turns)+std(n_turns))/2/sqrt(length(temp));
    
    %%% model
    mlee = squeeze((mle_params_opto(cc,:,5)));
    Kodor = mlee(3:6)*cosBasis';
    Kopto = mlee(7:10)*cosBasis';
    odor_data = extractfield(Data, 'dc');
    opto_data = extractfield(Data, 'opto');
    stim_pos = find(opto_data>4);  % find impulse position
    filt_dc = conv_kernel(odor_data, Kodor);  % fix this
    filt_opto = conv_kernel(opto_data, Kopto); % looking at impulse... if this does not work we only look ar data...
    A_ = mlee(2); C_ = mlee(11);
    P = (A_-C_) ./ (1 + exp( -((filt_dc(stim_pos)) + filt_opto(stim_pos)) )) + C_;  %sigmoid(A_,B_,dc); 
    Pbase = (A_-C_) ./ (1 + exp( -((filt_dc(stim_pos)) + filt_opto(stim_pos-9)) )) + C_;
    dP = P - Pbase;
    model_pr_mean(cc) = mean(dP)*14/5;
    model_pr_std(cc) = std(dP);%/(sqrt(length(stim_pos)));
    
end
data_pr_mean
model_pr_mean

clear ctr ydt
figure;
hBar = bar([data_pr_mean; model_pr_mean]');

for k1 = 1:2
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');    hold on  
    ydt(k1,:) = hBar(k1).YData;  hold on
end
hold on 
errorbar(ctr, ydt,  [data_pr_std; model_pr_std], '.k')
xticklabels({'appetitive', 'naive', 'aversive'});
set(gcf,'color','w'); set(gca,'Fontsize',20);
