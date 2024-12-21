% SI panels
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other scripts
%%% SI1
% numeric shown in S1_data.xlsx
%%% SI2
% find example in /odor_flow/
%%% SI3
% run /demo/mGLM.m script
%%% SI 12,13
% run /scripts/dPAW_gof_strains.m and mute_param_sim.m

%% S I4,5,7
%% load all data, with fully trained parameters
load(fullfile(datadir,'data4plots', 'classify_learn3_vars7.mat'));
datas = {fullfile(datadir,'data4plots', 'Data_app_test2.mat'),...
         fullfile(datadir,'data4plots', 'Data_nai_test2.mat'),...
         fullfile(datadir,'data4plots', 'Data_ave_test2.mat')};

%% empirical turning rate
model_Pturn = zeros(1,3);
exp_Pturn = zeros(1,3);
turn_thr = 150;
gamm = 0.20;
nbins = 50;
d2r = pi/180;

figure
for ii = 1:3
    %%% load data
    load(datas{ii});

    %%% test with the same CV testing data
    indices = ids{ii};  % indices for CV (pre-assigned!)
    test_set = (indices==kth);   %%% this matters!!
    train_set = ~test_set;
    Data_train = Data(train_set);
    %%%
    
    [xx_train, yy_train, mask_train] = data2xy(Data_train);%(Data(1:400));%(Data(idvec(1:400)));
    ddc_fit = xx_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
    ang_fit = yy_train;
    
    % emperical
    exp_Pturn(ii) = length(find(abs(ang_fit)>turn_thr)) / (length(ang_fit)-sum(isnan(trials_fit)));
    
    % model
    filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
    K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
    dc_dth = (filt_ddc + filt_dth + b_mle(ii)*0);
    Pturns = (A_mle(ii)-C_mle(ii))./ (1 + exp( -(dc_dth) )) + C_mle(ii);
    P_beta = nansum(Pturns) / (length(Pturns)-sum(isnan(Pturns)));
    [n, edges] = histcounts(yy_train, nbins, 'Normalization', 'probability');
    bb = edges(1:end-1)*pi/180;
    pos = find(abs(bb) > turn_thr*pi/180);
    scal = sum(1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)  + ... 
       ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta );
    p_thr_beta = (1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) * P_beta / scal;
    model_Pturn(ii) = sum(p_thr_beta(pos));
    exp_Pturn(ii) = sum(n(pos));

end

% figure
learn_bar_plot([exp_Pturn; model_Pturn], 'turn/s')

%% densities
nbins = 100;
ii = 2;
gamm = 0.2;
load(datas{ii});
figure
[xx_train, yy_train, mask_train] = data2xy(Data);
%%% turn P
filt_ddc = conv_kernel(ddc_fit.*trials_fit, Kc_mle(ii, :));
K_h_rec = Ah_mle(ii)*exp(-t_vec/th_mle(ii));
filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
dc_dth = (filt_ddc + filt_dth + b_mle(ii)*0);
Pturns = (A_mle(ii)-C_mle(ii))./ (1 + exp( -(dc_dth) )) + C_mle(ii);
P_beta = nansum(Pturns) / (length(Pturns)-sum(isnan(Pturns)));
% dth histogram
hh = histogram(yy_train, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
bb = hh.BinEdges(1:end-1)*pi/180;
scal = sum(1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)  + ... 
       ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta );
plot( bb*180/pi, 1/(2*pi*besseli(0,K1_mle(ii)^1)) * exp(K1_mle(ii)^1*cos( bb )) * (1-P_beta)*1/scal, 'b'); hold on
plot( bb*180/pi, ( 1/(2*pi*besseli(0,K2_mle(ii)^1)) * exp(K2_mle(ii)^1*cos( bb-pi ))*(gamm) + (1-gamm)/(2*pi) ) *P_beta*1/scal,'r');

%% concentration distribution
numBins = 30;
cols = ['b', 'k', 'r'];
figure
for ii = 1:3
    %%% load data
    load(datas{ii});
    dcs = zeros(1,length(Data));
    for jj = 1:length(Data)
        dcs(jj) = Data(jj).dc(end) - Data(jj).dc(1);
    end
    hh = histogram(dcs, numBins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% historgram of simulated tracks
%% low data
load(fullfile(datadir,'data4plots', 'Data_nai_test2.mat'));
mle = load(fullfile(datadir,'data4plots', 'Kfold_mle_param3.mat'));
mle = mle.mle_params;
Cmap = load(fullfile(datadir,'data4plots', 'Landscape_low_0623_2.mat'));
M = Cmap.vq1;
M = fliplr(flipud(M));

%% setup params
clear specs
specs = struct();

[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);

specs.M = M;
specs.fr = 14/5;
specs.cosBasis = cosBasis;
specs.T = floor(30*60*14/5);
specs.dt = 1;
specs.REP = 100;
cond_id = 1;
x = squeeze(mle(2,cond_id,:))';
x_temp= x;
[tracks, CI] = param2tracks(x_temp, specs, []);
x_temp(10) = 0;  % without history
[tracks_woh, CI] = param2tracks(x_temp, specs, []);

%% plots for tracks
thre = 50;
run_dur_model = [];
run_dur_data = [];
run_dur_woh = [];
for ii = 1:100
    dth_i = tracks(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_model = [run_dur_model   diff(pos)];
    dth_i = tracks_woh(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_woh = [run_dur_woh   diff(pos)];
    dth_i = Data(ii).dth;
    pos = find(abs(dth_i)>thre);
    run_dur_data = [run_dur_data   diff(pos)];
end

figure
subplot(121)
[xx_train, yy_train, mask_train] = data2xy(tracks); %tracks or Data
hist(yy_train,100)
subplot(122)
hist(run_dur_model,50);
figure
subplot(121)
[xx_train, yy_train, mask_train] = data2xy(Data);
hist(yy_train,100)
subplot(122)
hist(run_dur_data*14/5,50);

%% for duration density
figure
numBins = 40;
binEdges = linspace(0, 120, numBins + 1);
histogram(run_dur_data*14/5, binEdges, 'Normalization', 'probability',  'FaceAlpha', 0.5); hold on
histogram(run_dur_model, binEdges, 'Normalization', 'probability', 'FaceAlpha', 0.5);
histogram(run_dur_woh, binEdges, 'Normalization', 'probability',  'FaceAlpha', 0.5);

%% compute tc
[dur, ydata] = hist(run_dur_data,150);
[dur, ydata] = hist(run_dur_model,150);
y = ydata;
F = @(x,dur)x(1)*exp(-dur/x(2)) + x(3)*exp(-dur/x(4));
x0 = [100 10 1 50] ;
xunc = lsqcurvefit(F, x0, dur, y);
tlist = linspace(min(dur), max(dur));   % Plot Finer Resolution
tc = 1/(1/xunc(2)-1/xunc(4)) * log(xunc(1)/xunc(3)) / (5/14)

%% check concentration distributions
rng(1)
temp = load(fullfile(datadir,'data4plots', 'Kfold_mle_param7.mat'));
mle = temp.mle_params;
cond_id = 2;
specs.REP = 50;
x = squeeze(mle(4,cond_id,:))';
x_temp= x;
specs.T = floor(30*60*14/5);
[tracks, CI] = param2tracks(x_temp, specs, []);
c_sim = extractfield(tracks, 'dc');
dcp_sim = extractfield(tracks, 'dcp');
temp = load(fullfile(datadir,'data4plots', 'Data_nai_test2.mat'));
c_data = extractfield(Data, 'dc');
dcp_data = extractfield(Data, 'dcp');

%% compare simulated and data densities
nb = 20;
figure;
subplot(121)
mm = min([c_sim, c_data]); MM = max([c_sim, c_data]);
beds = linspace(mm,MM, nb);
histogram(c_data, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'black','FaceAlpha', 0.5);
hold on;
histogram(c_sim, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 0.5);
xlabel('ppm'); ylabel('probability'); set(gcf,'color','w'); set(gca,'Fontsize',20);
subplot(122)
mm = min([dcp_sim, dcp_data]); MM = max([dcp_sim, dcp_data]);
beds = linspace(mm,MM, nb);
histogram(dcp_data, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'black','FaceAlpha', 0.5);
hold on;
histogram(dcp_sim, 'BinEdges', beds,'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 0.5);
xlabel('ppm/mm'); ylabel('probability'); set(gcf,'color','w'); set(gca,'Fontsize',20);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SI10
%% load data
% redirect to load mutant tracks
load(fullfile(datadir,'track_data/mutants', 'rawTracks_AML580ave.mat'));
% load('/projects/LEIFER/Kevin/Publications/temp/track_data/mutants/rawTracks_AML580ave.mat')

%% plot tracks from mutatn chemotaxis
figure()
imagesc(M)
colormap()
hold on
for ii = 1:length(Tracks)
    xy = Tracks(ii).Path;
    plot(xy(:,1), xy(:,2),'k')
    hold on
    plot(xy(1,1), xy(1,2),'g.', 'MarkerSize',25)
    plot(xy(end,1), xy(end,2),'r.', 'MarkerSize',25)
end

%% local function
function learn_bar_plot(ap_na_av, y_name, hh)
    if nargin==2
        figure;
    else
        disp('hold on')
    end
    
    cols = ['b','k','r'];
    for ii = 1:3
        bar(ii, ap_na_av(:,ii), 'FaceColor',cols(ii));
        hold on
    end
    xticks(1:3);
    xticklabels( {'appetitive' 'naive', 'aversive'});
    % set(gca, 'XTickLabel', {'appetitive' 'naive', 'aversive'})
    ylabel(y_name)
    set(gcf,'color','w'); set(gca,'Fontsize',20);
end
