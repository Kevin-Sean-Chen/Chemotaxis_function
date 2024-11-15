% Bearing_turn
%%% given the fitted gams (P(z)), compute bearing conditioned on states
%%% here we focus on the analysis for bearing affected by turning events,
%%% rather than the states!

%% loading data
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_staPAW.mat')
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt100_50_staPAW.mat')

%% load trained models
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240319_110114_cv_staPWA.mat') %%% Data_salt0_50 %%%
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240324_114402_cv_staPWA.mat') %%% Data_salt100_50
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240520_023653_cv_staPWA.mat') %%% new Data_salt0_50_0513, for pirouettes demo!
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/20240527_060158_cv_staPWA.mat')
% rng(42) %37

%% asign models
rr = 1; %3%2
dPAW_fit = all_record(rr,1,1).params;  % single state
staPAW_fit = all_record(rr,2,1).params;  % two-state model
temp_ws = dPAW_fit.wts;
temp_ws([2:5]) = zeros(1,4);
temp_ws([14:17]) = ones(1,4);
% null_fit = dPAW_fit;
% null_fit.wts = temp_ws;  % ad-hoc ablation!

staPAW_wo_fit = staPAW_fit;
staPAW_wo_fit.wts_state = ones(2,2,4)*.0;
% staPAW_wo_fit.A = ones(2,2); %staPAW_wo_fit.A*0.5;

%%% one at a time
model_choice = staPAW_fit; %dPAW_fit;  %staPAW_wo_fit; %

%% processing data
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data(:));
alldis = extractfield(Data(:), 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:160000];  %150000, 250000, 298000  %150000:298500
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

target = [3000, 1500];  % gradient vector

%%
%%% replacing data back
xy_sim = xys;

%% defining states (change for dPAW comparison)
pre_t = 28; %4 %8, 2, 20
post_t = 2; %4
state_vec = gams_(1,:)*0; %gams_*0;%
%%% staPAW states!
pos_state_stapaw = find(gams_(1,:)>0.5); %<gams_(2,:));
%%% mock states!
pos_state_turn = find(abs(yy(1,:))>50);  % use this for fare comparison across models
% pos_state_turn = find(abs(dth_sim)>50);  %50, 18

%%% logics
% pos_state = setdiff(pos_state_turn, [pos_state_stapaw,pos_state_stapaw+1,pos_state_stapaw-1]);  % only turns not states
% pos_state = setdiff(pos_state_stapaw, pos_state_turn);  % only state-switches not turns

%%% emperical turn-based states
windp = 28;  %20,28
pos_purit = wind_test*0;
pos_purit(pos_state_turn) = ones(1,length(pos_state_turn));
pos_purit = conv(pos_purit,ones(1,windp), 'same');
pos_purit = find(pos_purit>1);

%%% more elaborate method to detect continued states, then remove late turn
%%% within state OR pirouette
diff_vec = diff(pos_purit); % pos_purit, pos_state_stapaw; % Find the differences between consecutive elements
chunk_starts = [1, find(diff_vec ~= 1) + 1]; % Identify the start of each new chunk (where difference is not 1)
chunk_ends = [chunk_starts(2:end) - 1, length(vec)];
chunks = arrayfun(@(s, e) pos_purit(s:e), chunk_starts, chunk_ends, 'UniformOutput', false); % Extract each chunk
%%
pos_turn_not_state = [];
for ii = 1:length(chunks)
    pos_i = chunks{ii};
    if length(pos_i)>1
        pos_turn_in_state = find(abs(yy(1,pos_i))>50);
        pos_turn_not_state = [pos_turn_not_state   pos_i(1)*1+pos_turn_in_state(2:end-1)];
    end
end
%%%

%%% assinging positions
% pos_state = pos_state_stapaw;  % staPAW
pos_state = pos_turn_not_state;  % turn within state, but not last
% pos_state = pos_state_turn;  % turns
% pos_state = randi(length(yy),1,length(pos_state_turn));  % random
fix_t = 10;

%%% test remove reversal
% pos = find(diff(pos_state)==1);
% pos_state(pos) = [];

%% test for classic pirouette! for data
% windp = 28;  %20,28
% pos_state_turn = find(abs(yy(1,:))>50);
% % pos_state_turn = find(abs(dth_sim)>50); 
% state_vec = wind_test*0;
% state_vec(pos_state_turn) = ones(1,length(pos_state_turn));
% state_vec = conv(state_vec,ones(1,windp), 'same');
% state_vec(find(state_vec>1)) = 1;

%% iterations
% state_vec = zeros(2,floor(length(xx)/2)+1);  % making 2d for column difference
% state_vec = state_vec(randperm(length(state_vec))); % randperm control
state_vec(pos_state) = ones(1,length(pos_state));  %%%% use for non-pirouette positions %%%
trans_pos = diff(state_vec);

trans12 = find(trans_pos>0)-0;
trans21 = find(trans_pos<0)+0;
npairs = min([length(trans12), length(trans21)]);
Bpairs = zeros(2,npairs);
vecs = diff(xy_sim,1,2);
for bb = 20:npairs-20
    %%% pre vector
    v1 = (target - xy_sim(:, trans12(bb))');  %target
%     v2 = vecs(:,trans12(bb)-pre_t)';   % vector based
    v2 = -(xy_sim(:,trans12(bb)) - xy_sim(:,trans12(bb)-pre_t))';  % position-based 
    v1 = [3000 v2(2)+0];
    
    Bpre = angles(v1, v2);

    %%% post vector
    v1 = (target - xy_sim(:, trans21(bb))'); %target
%     v2 = vecs(:,trans21(bb)+post_t)';
    v2 = -(xy_sim(:,trans21(bb)) - xy_sim(:,trans21(bb)+post_t))';

%     v1 = -(target - xys(:, trans12(bb)+fix_t)'); %target
%     v2 = (xys(:,trans12(bb)+fix_t) - xys(:,trans12(bb)+fix_t+post_t))';  % another timed-control!
    v1 = [3000 v2(2)+0];
    
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
%%
%%% NOTES
% should be slightly different from the old analysis
% the bearing pre-pirouette might be accounted by individual turns
% show that staPAW recap this

%%
% % Estimate PDFs using kernel density estimation
% [f1, x1] = ksdensity(Bpost_data); % PDF of data 1
% [f2, x2] = ksdensity(Bpost_dpaw); % PDF of data 2  %Bpost_stapaw
% % Compute KL divergence
% kl_divergence = sum(f2 .* log(f2 ./ f1))
