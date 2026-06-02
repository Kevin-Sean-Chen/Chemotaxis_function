%% Bearing analysis for Figure 5 (a-d)
%%% given the fitted gams (P(z)), compute bearing conditioned on states
clear
clc
% USER INSTRUCTION: download and unzip the demo data pack Chen_states_2024.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_states_2024');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_states_2024');
%% load data
load(fullfile(datadir,'data4plot/salt_chemotaxis', '20240520_023653_cv_staPWA.mat'));
% load(fullfile(datadir,'data4plot/salt_chemotaxis', 'Data_salt100_50_staPAW.mat'));

%% asign models
rr = 1;
dPAW_fit = all_record(rr,1,1).params;  % single state
staPAW_fit = all_record(rr,2,2).params;  % two-state model
temp_ws = dPAW_fit.wts;
temp_ws([2:5]) = zeros(1,4);
temp_ws([14:17]) = ones(1,4);
% null_fit = dPAW_fit;
% null_fit.wts = temp_ws;  % ad-hoc ablation!

staPAW_wo_fit = staPAW_fit;
staPAW_wo_fit.wts_state = ones(2,2,4)*.0;
% staPAW_wo_fit.A = ones(2,2); %staPAW_wo_fit.A*0.5;

%%% one at a time
model_choice = staPAW_fit; % for the full staPAW fitting
% model_choice = dPAW_fit;  %%% dPAW fit, loaded from file
% model_choice = staPAW_wo_fit; %%% staPAW without transitions

%% processing data
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data(:));
alldis = extractfield(Data(:), 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:160000];
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(model_choice,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end
xys = allxys(:,wind_test);

target = [3000, 1500];  % gradient vector for odor
% gams_ = ones(1,length(xx))+1;  %%% use this when it is the without states condition %%%%

%% if we want to sample from the fitted staPAW given states
% xy_sim = xys*0;
% dth_sim = zeros(1,length(xys));
% wind = 47 + 2;  % should be larger than the kernel size (47) along with prediction window (which should be larger than the pre-post predictions)
% xy_sim(:,wind) = xys(:,wind);
% 
% for tt = wind:length(xx)-1
%     [dth_pred, dr_pred] = model_predit_PAW(model_choice, xx(:,tt-wind+1:tt-0), yy(:,tt-wind+1:tt-0), gams_(:,tt-wind+1:tt-0));  % forward model
%     dr = dr_pred(end);
%     dth = dth_pred(end);
%     xy_sim(:,tt+1) = xys(:,tt) + ([dr*sin(dth*pi/180) dr*cos(dth*pi/180)])';  % make predictions about location (function of dr, dv) for one step
% %     xy_sim(:,tt+1) = xy_sim(:,tt) + ([dr*sin(dth*pi/180) dr*cos(dth*pi/180)])';  % without data
%     dth_sim(tt) = dth;
% end

%%
%%% replacing data back
xy_sim = xys;

%% defining states (change for dPAW comparison)
pre_t = 28; %8, 2, 20; in terms of 5/14 s units
post_t = 2;
state_vec = gams_(1,:)*0; %gams_*0;%
%%% staPAW states!
pos_state_stapaw = find(gams_(2,:)>0.5); %<gams_(2,:));
%%% mock states!
pos_state_turn = find(abs(yy(1,:))>50);  % use this for fare comparison across models
% pos_state_turn = find(abs(dth_sim)>50);  %50, 18

%%% logics (choose for different panels as model comparison)
pos_state = pos_state_stapaw; %%% for fitted staPAW
% pos_state = pos_state_turn;  %%% for dPAW along
% pos_state = randi(length(yy),1,length(pos_state_turn));  %%% for randomized control

% pos_state = setdiff(pos_state_turn, pos_state_stapaw);  % only turns not states
% pos_state = setdiff(pos_state_stapaw, pos_state_turn);  % only state-switches not turns

%%% test remove reversal
% pos = find(diff(pos_state)==1);
% pos_state(pos) = [];

%% test for classic pirouette from data
% windp = 18;  %20,28
% pos_state_turn = find(abs(yy(1,:))>50);
% % pos_state_turn = find(abs(dth_sim)>50); 
% state_vec = wind_test*0;
% state_vec(pos_state_turn) = ones(1,length(pos_state_turn));
% state_vec = conv(state_vec,ones(1,windp), 'same');
% state_vec(find(state_vec>1)) = 1;

%% iterations
% state_vec = zeros(2,floor(length(xx)/2)+1);  % making 2d for column difference
state_vec(pos_state) = ones(1,length(pos_state));  %%%% use for non-pirouette positions
% state_vec = state_vec(randperm(length(state_vec)));
trans_pos = diff(state_vec);

trans12 = find(trans_pos>0)-0;
trans21 = find(trans_pos<0)+0;
npairs = min([length(trans12), length(trans21)]);
Bpairs = zeros(2,npairs);
vecs = diff(xy_sim,1,2);
for bb = 20:npairs-20
    %%% pre vector
%     v1 = (target - xy_sim(:, trans12(bb))');  %target for odor
    v2 = -(xy_sim(:,trans12(bb)) - xy_sim(:,trans12(bb)-pre_t))';  % position-based 
    v1 = [3000 v2(2)+0]; %%% for linear gradients
    
    Bpre = angles(v1, v2);

    %%% post vector
%     v1 = (target - xy_sim(:, trans21(bb))'); %target for odor
    v2 = -(xy_sim(:,trans21(bb)) - xy_sim(:,trans21(bb)+post_t))';
%     v1 = -(target - xys(:, trans12(bb)+fix_t)'); %target
%     v2 = (xys(:,trans12(bb)+fix_t) - xys(:,trans12(bb)+fix_t+post_t))';  % another timed-control!
    v1 = [3000 v2(2)+0]; %%% for linear gradients
    
    Bpost = angles(v1, v2);
    
    %%% recording pairs
    Bpairs(:,bb) = [Bpre; Bpost];
end

%%
figure()
% hist(Bpairs(1,:),20, 'FaceColor', 'r'); hold on
% hist(Bpairs(2,:),20, 'FaceColor', 'b')
nbins = 12;
subplot(131)
histogram(Bpairs(1,:), nbins, 'FaceColor', 'r', 'EdgeColor', 'none', 'Normalization', 'probability');  % Use 20 bins for the histogram
hold on; xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20);

% subplot(132)
% histogram(Bpairs(2,:), nbins, 'FaceColor', 'r', 'EdgeColor', 'none', 'Normalization', 'probability');

subplot(133)
shuffle_rep = 100;
dB = -diff(Bpairs,1);
Bc  = [];
for ss = 1:shuffle_rep
    Bc = [Bc dB(randperm(length(dB))) + Bpairs(1,:)];
%     Bc = dB(randperm(length(dB))) + Bpairs(1,:);
end
Bc = wrapToPi(Bc*pi/180)*180/pi;
% histogram(Bc, nbins, 'FaceColor', 'k', 'EdgeColor', 'none', 'Normalization', 'probability'); hold on
H = histogram(Bc, nbins, 'Normalization', 'probability','Visible', 'off');
binEdges = H.BinEdges;  counts = H.BinCounts;
midb = (binEdges(2:end)+binEdges(1:end-1))/2;
plot(midb, counts/sum(counts)); hold on
histogram(Bpairs(2,:), nbins, 'FaceColor', 'b', 'EdgeColor', 'none', 'Normalization', 'probability'); hold on; xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20);

subplot(132)
dB = wrapToPi(dB*pi/180)*180/pi;
histogram(dB, nbins, 'FaceColor', 'g', 'EdgeColor', 'none', 'Normalization', 'probability'); xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20);

test = Bpairs(2,:);
length(find(abs(test)<90))/length(test)  % compute fraction that aligns

%% for proper distribution tests
% % Estimate PDFs using kernel density estimation
% [f1, x1] = ksdensity(Bpost_data); % PDF of data 1
% [f2, x2] = ksdensity(Bpost_dpaw); % PDF of data 2  %Bpost_stapaw
% % Compute KL divergence
% kl_divergence = sum(f2 .* log(f2 ./ f1))

