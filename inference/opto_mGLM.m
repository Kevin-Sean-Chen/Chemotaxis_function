%%% opto_GLM fitting with popluation tracks
clear
clc

%% load data
addpath('C:\Users\Kevin\Documents\GitHub\leifer-Behavior-Triggered-Averaging-Tracker_new\Experimental Analysis')
addpath('C:\Users\Kevin\Desktop\Chemotaxis_function')

%batch analysis
% fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothSpeed','AngSpeed','SmoothX','SmoothY','Behaviors', 'LEDVoltages'};
fields_to_load = {'Path','SmoothX','SmoothY', 'LEDVoltages','Time'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);
nn = length(Tracks);

%% load odor map
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat')

M = Cmap.vq1;
M = fliplr(flipud(M));  %flipped camera

%% retrieving data across tracks
clear Data
Data(1) = struct();  % all tracks and chemotaxis time series

%%% pre-processing parameters
min_t = 60*1;
min_x = 100;
poly_degree = 3;  %polynomial fit for the tracks
bin = 5;  %temporal binning  %~0.5 s
filt = 14*1;  %filtering tracks
l_window = 1;  %lag time
perp_dist = 1;  %perpendicular vectors for C^perp measurements
fr2sec = 1/14;
dis_thr = 5.5*fr2sec*30;  %max or min mm/s that is reasonable
%basis for filters
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);

%%% gather time series
Data(1).dth = []; %angle change
Data(1).opto = [];  %optogenetic input
Data(1).dis = [];  %displacements
Data(1).mask = [];  %mark for tracks
Data(1).xy = [];  %all track location in 2D
Data(1).theta = [];  %local angles
Data(1).time = [];  % real time in the assay

id = 1;
for ii = 1:nn
    
    %%% process path to angles and distances
    temp = Tracks(ii);  %for loading saved tracks
    temp1 = zeros(round(size(temp.Path,1)/1),2);
    temp1(:,1) = smooth(temp.SmoothX, filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp.SmoothY, filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    timev = temp.Time;
    timev = timev(1:bin:end);
    
    %%% checking criterion
    if temp1(end,1)>500 && temp1(end,1)<2500 && temp1(end,2)>500 && temp1(end,2)<2000 %&& ...
%        M(floor(temp1(1,2)), floor(temp1(1,1))) - M(floor(temp1(end,2)), floor(temp1(end,1))) < 0. % thresholding tracks
    if timev(end) - timev(1) > min_t
    displace = mean((Tracks(ii).Path(:,1)-mean(Tracks(ii).Path(:,1))).^2 + (Tracks(ii).Path(:,2)-mean(Tracks(ii).Path(:,2))).^2); %pixel displacement
    if displace > min_x^2 
        
    vecs = diff(subs);
%     vecs = [vecs(1,:); vecs];   % compensating for the size change/time shift after taking difference 
    angs = zeros(1,size(vecs,1));    
    ang_loc = zeros(1,size(vecs,1));
    opto_i = temp.LEDVoltages;
    opto_i = opto_i(1:bin:end);

    %%% for time step calculations
%     optot = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    trials = ones(1,size(vecs,1));
    xys = zeros(2,size(vecs,1));
    
    %%% iterate through worms
    for dd = (l_window+1):length(angs)
        %%% angle function
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        ang_loc(dd) = angles([1,0],vecs(dd,:)/norm(vecs(dd,:)));
        
        %%% concentration input
        dCs(dd) = M(floor(subs(dd-0,2)), floor(subs(dd-0,1)));
        
        %%% check displacement
        dds(dd) = norm(subs(dd,1)-subs(dd-l_window,1), subs(dd,2)-subs(dd-l_window,2));
        
        %%% record location in 2D
        xys(:,dd) = subs(dd,:)';
    end
    %remove zeros
    dds = dds((l_window+1):end);
    angs = angs((l_window+1):end);
    ang_loc = ang_loc(l_window+1:end);
    trials = trials((l_window+1):end);
    xys = xys(:, (l_window+1):end);
    timev = timev((l_window+2):end);  % one index off because of the derivative
    opto_i = opto_i((l_window+2):end);
    dCs = dCs((l_window+1):end);
    
    trials(1) = nan; trials(end) = nan;
    
    % store as structure
    Data(id).dth = angs; %angle change
    Data(id).opto = opto_i;  %perpendicular dC
    Data(id).dis = dds;  %displacements
    Data(id).mask = trials;  %mark for tracks
    Data(id).xy = xys;  %all track location in 2D
    Data(id).theta = ang_loc;  %local angles
    Data(id).Basis = cosBasis;
    Data(id).lambda = .1;
    Data(id).time = timev;
    Data(id).dc = dCs;
    id = id+1;
    end
    end
    end
end

%%
%     kappa_wv = THETA(1)^0.5;  % variance of von Mises
%     A_ = THETA(2);            % max turn probability
%     alpha_dc = THETA(3:6);    % kernel for dC transitional concentration difference (weights on opto kernel basis)
%     C_ = THETA(7);            % baseline turning rate
%     Amp_h = THETA(8);        % amplitude of turning history kernel
%     tau_h = THETA(9);        % time scale for dth history kernel
%     kappa_turn = THETA(10)^0.5;   % vairance of the sharp turn von Mises
%     gamma = THETA(11);        % weight for uniform angle in the turn

%% test with stats-model for kernels
drange = randperm(length(Data));%[1:length(Data)]; %
Data_fit = Data(drange(1:200));  %Data(1:100); %
lfun = @(x)pop_nLL_opto(x, Data_fit);

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [1e-0, 1e-1, ones(1,nB*2)*-inf, 0,    1., 0.1];
UB = [200, 1., ones(1,nB*2)*inf, 0.1    20, 1.];
prs0 = [50, 0.5, randn(1,nB*2)*1, 0.01,     5, .5];
prs0 = prs0 + prs0.*randn(1,length(UB))*0.2;
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

normedLL = fval / length(Data_fit);  % this is LL normalized by time and also by trails
x
normedLL

%% summary statistics
B_dc = x(3:6);
B_opto = x(7:10);
tt = [1:length(cosBasis)]*5/14;
figure;
plot(tt, B_dc * cosBasis'); hold on
plot(tt, B_opto * cosBasis')
xlabel('time (s)'); ylabel('weights');
set(gca,'FontSize',20); set(gcf,'color','w');

%% loop learned conditions
datas = {'/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_app.mat',...
         '/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_nai.mat',...
         '/projects/LEIFER/Kevin/Data_odor_opto/odor_opto_learn_data/Data_opto_ave.mat'};

rep = 5;
mle_params_opto = zeros(3,13, rep);
LB = [1e-0, 1e-1, ones(1,nB*2)*-inf, 0,    1., 0.1];%, -inf, -180];
UB = [200, 1., ones(1,nB*2)*inf, 0.1    20, 1.];%, inf, 180];
prs0 = [50, 0.5, randn(1,nB*2)*.1, 0.01,     5, .5];%, 0, 10];
    
for c = 1:3
    load(datas{c})
    for r = 1:rep
        drange = randperm(length(Data));
        Data_fit = Data(drange(1:200));
        lfun = @(x)pop_nLL_opto(x, Data_fit);
        opts = optimset('display','iter');
        prs0 = prs0 + prs0.*randn(1,length(UB))*0.;
        try
            [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
            mle_params_opto(c,:,r) = x;
        catch
            mle_params_opto(c,:,r) = zeros(1,13)*NaN;
        end
    end
    
end