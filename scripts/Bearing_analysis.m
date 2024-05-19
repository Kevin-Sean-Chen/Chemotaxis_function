% Bearing_analysis
%%% given the fitted gams (P(z)), compute bearing conditioned on states

%% loading
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_staPAW.mat')
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt100_50_staPAW.mat')

%% processing
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
yyf = [yyf; alldis];

wind_test = [1:200000];
offset = min(wind_test)-1;
xx = xxf(:,wind_test);
yy = yyf(:,wind_test);

maskf = alltrials;
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data); %Data_perm
allxys = [];
for ii = 1:length(Data);  allxys = [allxys  Data(ii).xy]; end
xys = allxys(:,wind_test);

target = [3000, 2500];  % gradient vector

%% 
pre_t = 8;
post_t = 8;
state_vec = gams_*0;
%%% staPAW states!
pos_state_stapaw = find(gams_(2,:)>0.5); %<gams_(2,:));
%%% mock states!
pos_state_turn = find(abs(yy(1,:))>18);

%%% logics
pos_state = setdiff(pos_state_turn, pos_state_stapaw);  % only turns not states
% pos_state = setdiff(pos_state_stapaw, pos_state_turn);  % only state-switches not turns

% pos_state = pos_state_stapaw; 
% pos_state = pos_state_turn;
fix_t = 10;

%% iterations
state_vec(pos_state) = ones(1,length(pos_state));
trans_pos = diff(state_vec);
trans12 = find(trans_pos>0)-0;
trans21 = find(trans_pos<0)+0;
npairs = min([length(trans12), length(trans21)]);
Bpairs = zeros(2,npairs);
vecs = diff(xys,1,2);
for bb = 2:npairs-1
    %%% pre vector
    v1 = -(target - xys(:, trans12(bb))');  %target
%     v2 = vecs(:,trans12(bb))';
    v2 = (xys(:,trans12(bb)) - xys(:,trans12(bb)-pre_t))';
%     v1 = -[3000 v2(2)];
    
    Bpre = angles(v1, v2);

    %%% post vector
    v1 = -(target - xys(:, trans21(bb))'); %target
%     v2 = -vecs(:,trans21(bb))';
    v2 = (xys(:,trans21(bb)) - xys(:,trans21(bb)+post_t))';

%     v1 = -(target - xys(:, trans12(bb)+fix_t)'); %target
%     v2 = (xys(:,trans12(bb)+fix_t) - xys(:,trans12(bb)+fix_t+post_t))';  % another timed-control!
%     v1 = -[3000 v2(2)];
    
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
dB = diff(Bpairs,1);
Bc = dB(randperm(length(dB))) + Bpairs(1,:);
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

%%
%%% NOTES
% should be slightly different from the old analysis
% the bearing pre-pirouette might be accounted by individual turns
% show that staPAW recap this

%%
function angle = angle_between(v1, v2)
    % Calculate the angle between two 2D vectors in the range of -pi to pi.
    if numel(v1) ~= 2 || numel(v2) ~= 2
        error('Vectors must be of length 2.');
    end
    
    dot_product = dot(v1, v2);
    magnitude_v1 = norm(v1);
    magnitude_v2 = norm(v2);
    
    if magnitude_v1 == 0 || magnitude_v2 == 0
        error('Vectors cannot have zero magnitude.');
    end
    
    angle = atan2(det([v1; v2]), dot_product);
end

