% Bear_staPAW_odor
%%%
% analyze bearing angle from turning, for odor landscape
%%%
%% load data
% load('/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/Data_nai_staPAW.mat') %%% odor naive
% load('/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/Data_ave_staPAW.mat') %%% odor aversive
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt100_50_staPAW2.mat') %%% 100-50 salt
% % % % load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_staPAW.mat') %%% 0-50 salt

% %% for mutants
load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIZ_ave_staPAW_vars.mat') %4500
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIB_app_staPAW_vars.mat')  %200000
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIY_app_staPAW_vars.mat') %130000

%%% new
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIB_nai_staPAW_vars2.mat') %% 300000
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIY_nai_staPAW_vars.mat') %% 130000

%% load landscape
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
% Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');
% 
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623.mat')
% Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat')
% 
% Fcon = Fcon.F;
% M = Cmap.vq1;
% M = fliplr(flipud(M));  %flipped camera

target = [1500  2500];  %[1250 2500]%for odor %[3000 1750];  % peak of odor landscape  %

%% for salt
% [rows, cols] = size(M);
% [x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
% % [x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
% gradient_x = x_ * 1;
% M = (y_*0+1) .* gradient_x;  figure; imagesc(M)

%% compile time seris
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

% testing
wind_test = [1:45000]; %[1:150000];%for nai, NaCl_100 %[1:200000];% [1:300000];%for ave
offset = min(wind_test)-1;
yy = yyf(:,wind_test);
xx = xxf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end
xys = allxys(:,wind_test);

%%
%%% replacing data back
xy_sim = xys;

%% defining states (change for dPAW comparison)
pre_t = 28; %8, 2, 12, 20
post_t = 20; %2 20
state_vec = gams_(1,:)*0; %gams_*0;%
%%% staPAW states!
pos_state_stapaw = find(gams_(1,:)>0.5); %<gams_(1,:));  %%% 2 for NaCl_100, 1 for others
%%% mock states!
pos_state_turn = find(abs(yy(1,:))>50);  % use this for fare comparison across models
% pos_state_turn = find(abs(dth_sim)>50);  %50, 18

pos_state = pos_state_stapaw; 
% pos_state = pos_state_turn;
% pos_state = randi(length(yy),1,length(pos_state_turn));
fix_t = 10;

%%% test remove reversal
% pos = find(diff(pos_state)==1);
% pos_state(pos) = [];

%% test for classic pirouette! for data
% windp = 8; %10
% pos_state_turn = find(abs(yy(1,:))>50);
% % pos_state_turn = find(abs(dth_sim)>50); 
% state_vec = wind_test*0;
% state_vec(pos_state_turn) = ones(1,length(pos_state_turn));
% state_vec = conv(state_vec,ones(1,windp), 'same');
% state_vec(find(state_vec>1)) = 1;

%% iterations
% state_vec = zeros(2,floor(length(xx)/2)+1);  % making 2d for column difference
state_vec(pos_state) = ones(1,length(pos_state));
% state_vec = state_vec(randperm(length(state_vec)));
trans_pos = diff(state_vec);

trans12 = find(trans_pos>0)-0;
trans21 = find(trans_pos<0)+0;
npairs = min([length(trans12), length(trans21)]);
Bpairs = zeros(2,npairs);
vecs = diff(xy_sim,1,2);
for bb = 8:npairs-8
    %%% pre vector
    v1 = (target - xy_sim(:, trans12(bb))');  %target
%     v2 = vecs(:,trans12(bb)-pre_t)';   % vector based
    v2 = -(xy_sim(:,trans12(bb)) - xy_sim(:,trans12(bb)-pre_t))';  % position-based 
%     v1 = [-3000 v2(2)+0]; %[-3000 v2(2)];
    
    Bpre = angles(v1, v2);

    %%% post vector
    v1 = (target - xy_sim(:, trans21(bb))'); %target
%     v2 = vecs(:,trans21(bb)+post_t)';
    v2 = -(xy_sim(:,trans21(bb)) - xy_sim(:,trans21(bb)+post_t))';

%     v1 = -(target - xys(:, trans12(bb)+fix_t)'); %target
%     v2 = (xys(:,trans12(bb)+fix_t) - xys(:,trans12(bb)+fix_t+post_t))';  % another timed-control!
%     v1 = [-3000 v2(2)+0]; %[-3000 v2(2)];
    
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
shuffle_rep = 250;
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
