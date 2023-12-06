%% sample_odor_opto
% other than using posterior sampling, try to fit to sub-tracks and compare
% odor-opto kernels

%% load tracks from each condition
% load data files
datas = {'/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_app.mat',...
         '/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_nai.mat',...
         '/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_ave.mat'};
     
%% fitting loop
rep = 50;
samp = 100;
cond = length(datas);  % all exp conditions
n_params = 13;  % number of parameters in our model for now
nB = 4;  % nunber of basis
mle_params = zeros(rep, cond, n_params); % K x c x N
samp_std_dc = zeros(2,cond,rep);  % storing variance of dC for normalization; (opto x condition x repeats)
rec_logp = zeros(cond,rep);

for rr = 1:rep
    for cc = 1:cond
        %%% fit kernels
        load(datas{cc})
        drange = randperm(length(Data));%[1:length(Data)]; %
        Data_fit = Data(drange(1:samp));  %Data(1:100); %
        lfun = @(x)pop_nLL_opto(x, Data_fit);

        opts = optimset('display','iter');
        LB = [1e-0, 1e-1, ones(1,nB*2)*-inf, 0,    1., 0.1];%, -inf, -180];
        UB = [200, 1., ones(1,nB*2)*inf, 0.1    50, 1.];%, inf, 180];
        prs0 = [50, 0.2, randn(1,nB*2)*1, 0.01,     5, .1];%, 0, 10];
        prs0 = prs0 + prs0.*randn(1,length(UB))*0.01;
        try
            [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
        end
        mle_params(rr,cc,:) = x;
        
        %%% record input variances
        alldc = extractfield(Data_fit, 'dc');
        alltrials = extractfield(Data_fit, 'mask');
        allopto = extractfield(Data_fit, 'opto');
        trials_fit = ones(1,length(alltrials));
        trials_fit(find(alltrials==0)) = NaN;
        samp_std_dc(1,cc,rr) = nanstd(alldc.*trials_fit);
        samp_std_dc(2,cc,rr) = nanstd(allopto.*trials_fit);
        
        rec_logp(cc,rr) = -fval/length(alldc);  % record likelihood
        
    end
end

%% test analysis for odor-opto kernels
figure
cond_title = {'appetitive', 'naive', 'aversive'};
for ci = 1:3
% ci = 3;
k_norms = zeros(2, rep);
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);

for ii = 1:rep
    k_norms(1,ii) = norm(squeeze(mle_params(ii,ci,3:6))'*cosBasis');% / samp_std_dc(1,ci,ii);  %sign(sum(m_samps(3:6,ii)))*
    k_norms(2,ii) = norm(squeeze(mle_params(ii,ci,7:10))'*cosBasis');% / samp_std_dc(2,ci,ii);
   
%     k_norms(1,ii) = sum(squeeze(mle_params(ii,ci,3:4))) / samp_std_dc(1,ci,ii);
%     k_norms(2,ii) = sum(squeeze(mle_params(ii,ci,7:8)))/ samp_std_dc(2,ci,ii);

end

X = k_norms(1,:);
Y = k_norms(2,:);
mdl = fitlm(k_norms(1,:)', k_norms(2,:)')
coefficients = mdl.Coefficients.Estimate;
slope = coefficients(2);
intercept = coefficients(1);
Y_predicted = slope * X + intercept;
subplot(1,3,ci)
scatter(X, Y, 'filled', 'MarkerFaceColor', 'k');
hold on;
plot(X, Y_predicted, '--', 'LineWidth', 2);
hold off;
set(gca,'FontSize',20); set(gcf,'color','w');
xlabel('K_{odor}'); ylabel('K_{opto}'); title(cond_title{ci})

rsquare = mdl.Rsquared.Ordinary;
rsquareString = ['R^2 = ' num2str(rsquare, 4) ', p = ' num2str(mdl.Coefficients.pValue(2))]; % Format the number to 4 decimal places;
legendString = sprintf(rsquareString);
legend(legendString);
% mdl.Coefficients.pValue
end

%% show example
samp_id = 3;
figure
for ci = 1:cond
    x = squeeze(mle_params(samp_id,ci,:))';
    B_dc = x(3:6);
    B_opto = x(7:10);
    tt = [1:length(cosBasis)]*5/14;
    subplot(1,3,ci)
    yyaxis left; 
    plot(tt, B_dc * cosBasis'); ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
    ylabel('K_{dc}');
    yyaxis right;
    plot(tt, B_opto * cosBasis'); yliml = get(gca,'Ylim');
    ylabel('K_{opto}'); 
    if yliml(2)*ratio<yliml(1)
        set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
    else
        set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
    end

    xlabel('time (s)');
    set(gca,'FontSize',20); set(gcf,'color','w'); title(cond_title{ci})
end