%% mGLM fitting with population of tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load tracks and maps
% Odor map
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');
Fcon = Fcon.F;
M = Cmap.vq1;
M = fliplr(flipud(M));  %flipped camera
H = fspecial('average',700);
M = imfilter(M, H, 'replicate');

%batch analysis
fields_to_load = {'Path','Time','Runs','Pirouettes','SmoothX','SmoothY'};%,'Behaviors'};
folder_names = getfoldersGUI();
Tracks = loadtracks(folder_names,fields_to_load);

cand = 1:length(Tracks);  % if all

%% retrieving data across tracks
clear Data
Data(length(cand)) = struct();  % all tracks and chemotaxis time series

%%% pre-processing parameters
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
Data(1).dcp = [];  %perpendicular dC
Data(1).dc = [];  %tangential dC
Data(1).dis = [];  %displacements
Data(1).mask = [];  %mark for tracks
Data(1).xy = [];  %all track location in 2D
Data(1).theta = [];  %local angles
Data(1).time = [];  % real time in the assay

ntracks = length(cand);  %where cand is candidates checked somewhere else
for ii = 1:ntracks
    
    %%% process path to angles and distances
    temp = Tracks(cand(ii));  %for loading saved tracks
    temp1 = zeros(round(size(temp.Path,1)/1),2);
    temp1(:,1) = smooth(temp.SmoothX, filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp.SmoothY, filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    timev = temp.Time;
    timev = timev(1:bin:end);
    vecs = diff(subs);
%     vecs = [vecs(1,:); vecs];   % compensating for the size change/time shift after taking difference 
    angs = zeros(1,size(vecs,1));    
    ang_loc = zeros(1,size(vecs,1));

    %%% for time step calculations
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    trials = ones(1,size(vecs,1));
    xys = zeros(2,size(vecs,1));
    
    %%% iterate through worms
    for dd = (l_window+1):length(angs)
        %%% angle function
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        ang_loc(dd) = angles([1,0],vecs(dd,:)/norm(vecs(dd,:)));
        
        %%% perpendicular concentration change
        perp_dir = [-vecs(dd-0,2), vecs(dd-0,1)];
        perp_dir = perp_dir/norm(perp_dir)*perp_dist;
        [xil, yil] = plate_bound(M, subs(dd-0,1)+perp_dir(1), subs(dd-0,2)+perp_dir(2));
        [xir, yir] = plate_bound(M, subs(dd-0,1)-perp_dir(1), subs(dd-0,2)-perp_dir(2));
        dCp(dd) = (M(yil,xil) - M(yir,xir));
%         dCp(dd) = (M(floor(subs(dd-l_window,2)+perp_dir(2)*perp_dist), floor(subs(dd-l_window,1)+perp_dir(1)*perp_dist))...
%                  - M(floor(subs(dd-l_window,2)-perp_dir(2)*perp_dist), floor(subs(dd-l_window,1)-perp_dir(1)*perp_dist))) / 1;  %computing normal direction dcp

        %%% forward concentration
        dCs(dd) = M(floor(subs(dd-0,2)), floor(subs(dd-0,1)));
        
        %%% check displacement
        dds(dd) = norm([subs(dd,1)-subs(dd-l_window,1), subs(dd,2)-subs(dd-l_window,2)]);
        
        %%% record location in 2D
        xys(:,dd) = subs(dd,:)';

    end
    %remove zeros
    dCp = dCp((l_window+1):end);
    dCs = dCs((l_window+1):end);
    dds = dds((l_window+1):end);
    angs = angs((l_window+1):end);
    ang_loc = ang_loc(l_window+1:end);
    trials = trials((l_window+1):end);
    xys = xys(:, (l_window+1):end);
    timev = timev((l_window+2):end);  % one index off because of the derivative
    
%     pos = find(dds>dis_thr); %remove sharp jumps
%     trials(pos) = nan;
    trials(1:50) = nan; trials(end) = nan;
    
    % store as structure
    Data(ii).dth = angs; %angle change
    Data(ii).dcp = dCp;  %perpendicular dC
    Data(ii).dc = dCs;  %tangential dC
    Data(ii).dis = dds;  %displacements
    Data(ii).mask = trials;  %mark for tracks
    Data(ii).xy = xys;  %all track location in 2D
    Data(ii).theta = ang_loc;  %local angles
    Data(ii).Basis = cosBasis;
    Data(ii).lambda = .1;
    Data(ii).time = timev;
    
end

%%
%     K_ = THETA(1);     % variance of von Mises
%     A_ = THETA(2);     % max turn probability
%     B_ = THETA(3:8);   % kernel for dC transitional concentration difference (weights on kerenl basis)
%     C_ = THETA(9);     % baseline turning rate
%     Amp = THETA(10);   % amplitude for kernel for dcp normal concentration difference
%     tau = THETA(11);   % time scale for dcp kernel
%     K2_ = THETA(12);   % vairance of the sharp turn von Mises
%     gamma = THETA(13); % uniform angle weight

%% test with stats-model for kernels
drange = randperm(length(Data));%[1:length(Data)]; %
Data_fit = Data(drange(1:100));  %Data(1:100); %
lfun = @(x)pop_nLL(x, Data_fit);

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [1e-0, 1e-1, ones(1,nB)*-inf, 0 -inf, 1e-0, -inf, 1e-1, 1., 0.05];%, -inf, -180];
UB = [200, 1., ones(1,nB)*inf, 0.1, inf, 50, inf, 100, 20, 1.];%, inf, 180];
prs0 = [50, 0.1, randn(1,nB)*10, 0.01, -1, 25, 1, 25, 5, 1.];%, 0, 10];
prs0 = prs0 + prs0.*randn(1,length(UB))*0.2;
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

normedLL = fval / length(Data_fit);  % this is LL normalized by time and also by trails
x
normedLL

%% summary statistics
% inferred parameters
K_ = x(1); A_ = x(2); B_ = x(3:6); C_ = x(7); Amp = x(8); tau = x(9); Amp_h = x(10); tau_h = x(11); K2_ = x(12);  gamma = x(13); %base_dc = x(14); base_dcp = x(15);
base_dc = 0;  base_dcp = 0;

% unwrap data
id = 1;
dcp_fit = Data(id).dcp;
ang_fit = Data(id).dth;
ddc_fit = Data(id).dc;

figure
subplot(2,2,1)
xx = 0:length(cosBasis)-1;
yyaxis left; plot(Amp*exp(-xx/tau)); hold on
yyaxis right; plot(Amp_h*exp(-xx/tau_h))
title('\delta C^{\perp}, \delta \theta kernel')
subplot(2,2,2)
plot(B_ * cosBasis')
title('\delta C kernel')

subplot(2,2,3)
K_dcp_rec = Amp*exp(-xx/tau);
filt_dcp = conv_kernel(dcp_fit, K_dcp_rec);%conv(dcp_fit, fliplr(Amp*exp(-xx/tau)), 'same');
[aa,bb] = hist((ang_fit - filt_dcp - base_dcp)*pi/180 , 500);
bar( bb, 1/(2*pi*besseli(0,K_^1)) * exp(K_^1*cos( bb )) , 100); hold on
bar( bb, 1/(2*pi*besseli(0,K2_^1)) * exp(K2_^1*cos( bb-pi ))*(gamma) + (1-gamma)/(2*pi) , 100,'r');
title('von Mises for \delta C^{\perp}')
subplot(2,2,4)
K_dc_rec = B_*cosBasis';
filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
xx_h = 1:length(xx)*1;
K_h_rec = Amp_h*exp(-xx_h/tau_h);
filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
dc_dth = filt_ddc + 1*filt_dth;
Pturns = (A_-C_) ./ (1 + exp( -(dc_dth + base_dc))) + C_; %+sb
plot(dc_dth/length(K_dcp_rec) , Pturns,'o')
title('Logistic for \delta C')

%%% printing
disp(['K_=',num2str(K_),' K2_',num2str(K2_),'beta',num2str(sum(B_)),'alpha',num2str(Amp)])

%% color code with track LL
LLs = zeros(1,length(Data_fit));  % tracks used for fitting
for ii = 1:length(LLs)
    LLs(ii) = pop_nLL(x, Data_fit(ii))/length(Data_fit(ii).dth);
    %%% test with ration!
%     LLs(ii) = (pop_nLL(x_ave, Data_fit(ii))) / (pop_nLL(x, Data_fit(ii)));
end
normV = [(LLs-min(LLs))./(max(LLs)-min(LLs))]';
% blue to red. 
C = [normV normV normV];%[normV zeros(size(normV)) 1-normV];

figure();
imagesc(M); hold on;
delta_C = LLs*0;
for ii = 1:length(LLs) %10:60 %200:300
    plot(Data_fit(ii).xy(1,:),Data_fit(ii).xy(2,:), 'Color', C(ii,:)); 
    hold on
    plot(Data_fit(ii).xy(1,1),Data_fit(ii).xy(2,1),'g.', 'MarkerSize',15)
    plot(Data_fit(ii).xy(1,end),Data_fit(ii).xy(2,end),'r.', 'MarkerSize',15)
    delta_C(ii) = M(floor(Data_fit(ii).xy(2,end)),floor(Data_fit(ii).xy(1,end))) - M(floor(Data_fit(ii).xy(2,1)),floor(Data_fit(ii).xy(1,1)));
end

%% subset of tracks
cand_list = [37,29,33];%[94,65,87]; % 102722 droplet data as an example
figure();
imagesc(M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on;
for ll = 1:length(cand_list)
    i = cand_list(ll);
    plot(Data(i).xy(1,:)*pix2mm,Data(i).xy(2,:)*pix2mm, 'Color', C(i,:),'LineWidth',2); 
    hold on
    plot(Data(i).xy(1,1)*pix2mm,Data(i).xy(2,1)*pix2mm,'g.', 'MarkerSize',15)
    plot(Data(i).xy(1,end)*pix2mm,Data(i).xy(2,end)*pix2mm,'r.', 'MarkerSize',15)
    
end
axis([20,100 -10,100]); set(gca,'YDir','normal')

%% color code with strategy (BRW incidents)
id = 50;
dcp_fit = Data(id).dcp;
ang_fit = Data(id).dth;
ddc_fit = Data(id).dc;
track_fit = Data(id).xy;

filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
dc_dth = filt_ddc + filt_dth;
Pturns = A_ ./ (1 + exp( -(dc_dth + base_dc))) + C_;
turn_pos = find(Pturns>0.99);
figure();
imagesc(M); hold on;
plot(track_fit(1,:), track_fit(2,:),'k','LineWidth',1); 
hold on
plot(track_fit(1,1), track_fit(2,1),'g.', 'MarkerSize',15)
plot(track_fit(1,end), track_fit(2,end),'r.', 'MarkerSize',15)
plot(track_fit(1,turn_pos), track_fit(2,turn_pos),'b.')

%% subset of track
cand_list = [37,29,33];%[94,65,87]; % 102722 droplet data as an example
figure();
imagesc(M','XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); hold on;
for ll = 1:length(cand_list)
    %%% load track
    i = cand_list(ll);
    dcp_fit = Data(i).dcp;
    ang_fit = Data(i).dth;
    ddc_fit = Data(i).dc;
    track_fit = Data(i).xy;
    %%% identify turns
    filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
    dc_dth = filt_ddc + filt_dth;
    Pturns = A_ ./ (1 + exp( -(dc_dth + base_dc))) + C_;
    turn_pos = find(Pturns>0.99);
    %%% plot tracks with LL color
    plot(Data(i).xy(1,:)*pix2mm,Data(i).xy(2,:)*pix2mm, 'Color', C(i,:),'LineWidth',2); 
%     plot(Data(i).xy(1,:)*pix2mm,Data(i).xy(2,:)*pix2mm, 'Color', 'k','LineWidth',2); 
    hold on
    plot(Data(i).xy(1,1)*pix2mm,Data(i).xy(2,1)*pix2mm,'g.', 'MarkerSize',15)
    plot(Data(i).xy(1,end)*pix2mm,Data(i).xy(2,end)*pix2mm,'r.', 'MarkerSize',15)
    %%% plot turn incidents
    plot(Data(i).xy(1,turn_pos)*pix2mm, Data(i).xy(2,turn_pos)*pix2mm,'cyan.')
end

axis([20,100 -10,100]); set(gca,'YDir','normal')
